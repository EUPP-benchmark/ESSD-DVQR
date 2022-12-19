# select library (path may needs to be adapted!)
.libPaths("/home/jobst/miniconda3/envs/r_env/lib/R/library") 

# load packages
library(doMC)
library(lubridate)
library(vinereg)

# load data (paths may need to be adapted!)
load("/home/jobst/benchmark/benchmark_t2m_train_ext.Rdata")
load("/home/jobst/benchmark/benchmark_t2m_test_ext.Rdata")

# select parameters for cluster (parameters may need to be adapted!)
s <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # get array number from batch script 
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")) # get cores from batch script
registerDoMC(cores = cores)


# D-vine copula postprocessing function with weekly training period
dvine_pp <- function(tdata, 
                     vdata, 
                     tresponse, 
                     tpredictors,
                     vpredictors,
                     family_set = "all",
                     selcrit = "bic", 
                     order = NA,
                     location, 
                     leadtime,
                     ffd, 
                     lfd,
                     nweeks,
                     cores = 1){
  
  ###-----------------------------------------------------------------------------
  ### Input
  # tdata .............. training  (data.frame)
  # vdata .............. validation data (data.frame)
  # tresponse .......... response variable for the training data (character)
  # tpredictors ........ predictor variables for the training data (character)
  # vpredictors ........ predictor variables for the validation/test data (character)
  # family_set ......... copula families (character)
  # selcrit ............ predictor variable selection criterion (character)
  # order .............. predictor variable order (character, logical)
  # location ........... location/station-id (character)
  # leadtime ........... leadtime (integer)
  # ffd ................ first forecast date (character)
  # lfd ................ last forecast date (character)
  # nweeks ............. see repo documentation (integer)
  # cores .............. number of cores for multi-core usage (integer)
  ###-----------------------------------------------------------------------------
  ### Output (list)
  # tquantiles ......... predicted quantiles for the training data (data.frame)
  # vquantiles ......... predicted quantiles for the validation/test data (data.frame)
  # runtime.train ...... training time for the model (vector)
  # runtime.pp ......... time for producing the postprocessed ensemble forecasts (vector)
  # tcrps .............. CRPS in the training data (data.frame)
  ###-----------------------------------------------------------------------------
  
  # forecast days
  days <- seq(as.Date(ffd), as.Date(lfd), by = "day")
  # weeks belonging to the forecast days
  w <- week(days)
  # evaluate last day 365/2 = 52, rest 1 with week 52; set week 53 to 52 to have more data points
  w[w == 53] <- 52
  # weeks, for each a model is fitted
  unique_w <- unique(w) 
  
  # get valid_date, response, predictors and weeks for location and leadtime in training data 
  tloc <- tdata[tdata$id == location & tdata$leadtime == leadtime, c("valid_date", tresponse, tpredictors)]
  w_tloc <- week(tloc$valid_date)
  w_tloc[w_tloc == 53] <- 52
  tloc <- cbind(w_tloc, tloc)
  tloc.all <- tloc
  
  # get valid_date, response, predictors and weeks for location and leadtime in test = validation data 
  vloc <- vdata[vdata$id == location & vdata$leadtime == leadtime, c("valid_date", vpredictors)]
  w_vloc <- week(vloc$valid_date)
  w_vloc[w_vloc == 53] <- 52
  vloc <- cbind(w_vloc, vloc)
  
  # initialize output
  runtime.train <- runtime.pp <- c()
  tcrps <- data.frame()
  tquantiles <- data.frame()
  vquantiles <- data.frame()
  
  # formula for model
  formula <- as.formula(paste0(tresponse,  " ~ ."))
  
  # iterate through each week in unique_w
  for (k in 1:length(unique_w)) {
    
    # select the weeks for training
    w_choose <- (unique_w[k] + -nweeks:nweeks) %% 52
    w_choose[w_choose == 0] <- 52
    
    # determine final training data for the location
    loc <- tloc[tloc$w_tloc %in% w_choose, c(tresponse, tpredictors)]
    loc <- na.omit(loc)
    
    # only pp, if at least 2 observations are available  
    if(nrow(loc) >= 2) {
      
      start_tm <- Sys.time()
      # fit model
      fit <- vinereg(formula = formula,
                     data = loc, 
                     family_set = family_set,
                     selcrit = selcrit, 
                     order = order,
                     cores = cores)
      end_tm <- Sys.time()
      
      # add training time
      runtime.train <- c(runtime.train,
                         as.numeric(difftime(end_tm, start_tm, units = "mins")))
      
      # validate model by CRPS and quantiles in training data (in sample)
      t.loc.dat <-  tloc.all[tloc.all$w_tloc == unique_w[k], ]
      y <- t.loc.dat[, tresponse]
      x <- predict(fit, newdata = t.loc.dat, alpha = 1:51/52, cores = cores)
      colnames(x) <- paste0("quantile.", 1:51)
      crps.k <- crps(y = y, x = x)
      tcrps <- rbind(tcrps, data.frame(week = t.loc.dat$w_tloc,
                                       valid_date = t.loc.dat$valid_date,
                                       crps = crps.k))
      tquantiles <- rbind(tquantiles, data.frame(week = t.loc.dat$w_tloc,
                                                 valid_date = t.loc.dat$valid_date,
                                                 x))
      
      # predict quantiles in test data (out of sample)  
      v.loc.dat <- vloc[vloc$w_vloc == unique_w[k], ]    
      start_tm <- Sys.time()
      x <- predict(fit, newdata = v.loc.dat, alpha = 1:51/52, cores = cores)
      end_tm <- Sys.time()
      colnames(x) <- paste0("quantile.", 1:51)
      
      # add postprocessing time
      runtime.pp <- c(runtime.pp, 
                      as.numeric(difftime(end_tm, start_tm, units = "mins")))
      vquantiles <- rbind(vquantiles, data.frame(week = v.loc.dat$w_vloc,
                                                 valid_date = v.loc.dat$valid_date,
                                                 x))
    } else { # if only 1 observation is available per week, set values to NaN
      
      runtime.train <- c(runtime.train, NaN)
      runtime.pp <- c(runtime.pp, NaN)
      
      # training data evaluation
      t.loc.dat <-  tloc.all[tloc.all$w_tloc == unique_w[k], ]
      tcrps <- rbind(tcrps, data.frame(week = t.loc.dat$w_tloc,
                                       valid_date = t.loc.dat$valid_date,
                                       crps = NaN))
      tquantiles <- rbind(tquantiles, data.frame(week = t.loc.dat$w_tloc,
                                                 valid_date = t.loc.dat$valid_date,
                                                 matrix(NaN, nrow = nrow(t.loc.dat), ncol = 51,
                                                        dimnames = list(c(), paste0("quantile.", 1:51)))))
      # test data evaluation (validation data)
      v.loc.dat <- vloc[vloc$w_vloc == unique_w[k], ]
      vquantiles <- rbind(vquantiles, data.frame(week = v.loc.dat$w_vloc,
                                                 valid_date = v.loc.dat$valid_date,
                                                 matrix(NaN, nrow = nrow(v.loc.dat), ncol = 51,
                                                        dimnames = list(c(), paste0("quantile.", 1:51)))))
    }
    
  }
  
  # prepare output
  out <- list(tquantiles = tquantiles[order(tquantiles$valid_date), ],
              vquantiles = vquantiles[order(vquantiles$valid_date), ],
              runtime.train = runtime.train,
              runtime.pp = runtime.pp,
              tcrps =  tcrps[order(tcrps$valid_date), ])
  
  return(out)
  
}


# function for CRPS calculation
crps <- function(y, x) {
  
  ###-----------------------------------------------------------------------------
  ### Input
  # y .............. observations (n-dim. vector)
  # x .............. ensemble forecasts (nxm-dim. data.frame)
  ###-----------------------------------------------------------------------------
  ### Output 
  # CRPS ......... CRPS (n-dim. vector)
  ###-----------------------------------------------------------------------------
  
  # get necessary parameters
  m <- ncol(x)
  i1 <- rep(1:m, times = m)
  i2 <- rep(1:m, each = m)
  
  # CRPS calculation
  out <- rowSums(abs(x-y))/m - 1/(2*m^2) * rowSums(abs(x[, i1]-x[, i2]))
  
  return(out)
}

# D-vine copula postprocessing setting
tdata <- benchmark_t2m_train_ext
vdata <- benchmark_t2m_test_ext
tresponse <- vresponse <- "observation"
tpredictors <- vpredictors <- c("ctrl", "mean")   
family_set <- "all"
selcrit <- "aic" 
order <- NA
nweeks <- rep(seq(0, 10, 2), each = 21)
leadtime <- rep(seq(0, 120, 6), 6)
lt <- leadtime[s]
dates <- unique(vdata[vdata$leadtime == lt, "valid_date"])
nw <- nweeks[s]
ids <- unique(tdata$id)

# start evaluation for each station for a fixed leadtime
results.dvine <- foreach(l = 1:229) %dopar% {
  
  out <- dvine_pp(tdata, 
                  vdata,
                  tresponse, 
                  tpredictors, 
                  vpredictors, 
                  family_set = family_set,
                  selcrit = selcrit, 
                  order = order,                 
                  location = ids[l], 
                  leadtime = lt,
                  ffd = dates[1], 
                  lfd = dates[length(dates)],
                  nweeks = nw,
                  cores = 1)
  out
  
}

# save output (path may needs to be adapted!)
nam <- paste0("/home/jobst/benchmark/results/results.dvine.lt", lt, ".nw", nw, ".Rdata")
save(results.dvine, file = nam)
