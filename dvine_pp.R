# select library
.libPaths("/home/jobst/miniconda3/envs/r_env/lib/R/library") 

# load packages
library(doMC)
library(lubridate)
library(vinereg)
library(rvinecopulib)

# load data
load("/home/jobst/benchmark/benchmark_t2m_train_ext.Rdata")
load("/home/jobst/benchmark/benchmark_t2m_test_ext.Rdata")

# select parameters for cluster
s <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # get array number from batch script 
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")) # get cores from batch script
registerDoMC(cores = cores)


# D-vine copula postprocessing function with weekly training period
dvine_pp <- function(tdata, vdata, 
                     tresponse, tpredictors,
                     vresponse, vpredictors,
                     family_set = "all",
                     selcrit = "bic", 
                     order = NA,
                     location, 
                     leadtime,
                     ffd, 
                     lfd,
                     nweeks,
                     cores = 1){
  
  days <- seq(as.Date(ffd), as.Date(lfd), by = "day")
  w <- week(days)
  # evaluate last day 365/2 = 52, rest 1 with week 52
  w[w == 53] <- 52
  unique_w <- unique(w) # weeks, for each a model is fitted
  
  # get valid_date, response, predictors and weeks for location and leadtime in training data 
  tloc <- tdata[tdata$id == location & tdata$leadtime == leadtime, c("valid_date", tresponse, tpredictors)]
  w_tloc <- week(tloc$valid_date)
  w_tloc[w_tloc == 53] <- 52
  tloc <- cbind(w_tloc, tloc)
  tloc.all <- tloc
  
  # get valid_date, response, predictors and weeks for location and leadtime in test = validation data 
  vloc <- vdata[vdata$id == location & vdata$leadtime == leadtime, c("valid_date", vresponse, vpredictors)]
  w_vloc <- week(vloc$valid_date)
  w_vloc[w_vloc == 53] <- 52
  vloc <- cbind(w_vloc, vloc)
  
  # initialize output
  dvine.order <- list()
  dvine.copulas <- list()
  runtime.train <- runtime.pp <- c()
  tcrps <- data.frame()
  tquantiles <- data.frame()
  vquantiles <- data.frame()
  
  formula <- as.formula(paste0(tresponse,  " ~ ."))
  
  for (k in 1:length(unique_w)) {
    
    w_choose <- (unique_w[k] + -nweeks:nweeks) %% 52
    w_choose[w_choose == 0] <- 52
    loc <- tloc[tloc$w_tloc %in% w_choose, c(tresponse, tpredictors)]
    loc <- na.omit(loc)
    
    # only pp, if at least 2 observations are available  
    if(nrow(loc) >= 2) {
      
      start_tm <- Sys.time()
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
      
      # get info of D-vine copula
      dvine.order[[k]] <- fit$order
      dvine.copulas[[k]] <- unlist(get_all_families(fit$vine))
      
      
      # validate model by crps and quantiles in training data (in sample)
      t.loc.dat <-  tloc.all[tloc.all$w_tloc == unique_w[k], ]
      y <- t.loc.dat[, tresponse]
      x <- predict(fit, newdata = t.loc.dat, alpha = 1:51/52, cores = cores)
      colnames(x) <- paste0("quantile.", 1:51)
      crps.k <- crps(y = y, x = x)
      tcrps <- rbind(tcrps, data.frame(week = t.loc.dat$w_tloc,
                                       valid_date = t.loc.dat$valid_date,
                                       nobs = nrow(loc),
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
      dvine.order[[k]] <- NaN
      dvine.copulas[[k]] <- NaN
      
      # training data evaluation
      t.loc.dat <-  tloc.all[tloc.all$w_tloc == unique_w[k], ]
      tcrps <- rbind(tcrps, data.frame(week = t.loc.dat$w_tloc,
                                       valid_date = t.loc.dat$valid_date,
                                       nobs = nrow(loc),
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
  
  out <- list(tquantiles = tquantiles[order(tquantiles$valid_date), ],
              vquantiles = vquantiles[order(vquantiles$valid_date), ],
              order = dvine.order,
              copulas = dvine.copulas,
              runtime.train = runtime.train,
              runtime.pp = runtime.pp,
              tcrps =  tcrps[order(tcrps$valid_date), ])
  
  return(out)
  
}

crps <- function(y, x) {
  
  m <- ncol(x)
  i1 <- rep(1:m, times = m)
  i2 <- rep(1:m, each = m)
  
  out <- rowSums(abs(x-y))/m - 1/(2*m^2) * rowSums(abs(x[, i1]-x[, i2]))
  return(out)
}

# D-vine copula postprocessing settings
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


results.dvine <- foreach(l = 1:229) %dopar% {
  
  out <- dvine_pp(tdata, vdata,
                  tresponse, tpredictors, 
                  vresponse, vpredictors, 
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

nam <- paste0("/home/jobst/benchmark/results/results.dvine.lt", lt, ".nw", nw, ".Rdata")
save(results.dvine, file = nam)



