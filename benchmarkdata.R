# load packages
library(ncdf4)
library(doMC)
registerDoMC(cores = detectCores())

# set working directory
setwd("./data/")

################################################################################

# PREPARE TRAINING DATA

# read data
data_fc <- nc_open("ESSD_benchmark_training_data_forecasts.nc")
data_obs <- nc_open("ESSD_benchmark_training_data_observations.nc")

# observations and forecasts of t2m
t2m_obs <- ncvar_get(data_obs, "t2m")
t2m_fc <- ncvar_get(data_fc, "t2m")

# station information
id <- as.character(data_fc[["dim"]][["station_id"]][["vals"]])
station_name <- ncvar_get(data_fc, "station_name")
station_altitude <- ncvar_get(data_fc, "station_altitude")
station_latitude <- ncvar_get(data_fc, "station_latitude")
station_longitude <- ncvar_get(data_fc, "station_longitude")
station_land_usage <- ncvar_get(data_fc, "station_land_usage")

# model information
model_altitude <- ncvar_get(data_fc, "model_altitude")
model_latitude <- ncvar_get(data_fc, "model_latitude")
model_longitude <- ncvar_get(data_fc, "model_longitude")
model_land_usage <- ncvar_get(data_fc, "model_land_usage")
ens_names <- c("ctrl", paste("ens", 1:10, sep = ""))

# time information
lead_times <- data_fc[["dim"]][["step"]][["vals"]] 
steps <- data_fc[["dim"]][["time"]][["vals"]]

# steps + 1, because of variable steps starts at 0
(reforecast_dates <- seq(as.Date("2017-01-02"), as.Date("2018-12-31"), by = "day")[steps+1])
n_dates <- length(reforecast_dates)
# check if reforecast_dates are on mondays and tuesdays; reforecasts dates can be seen as reference dates
weekdays(reforecast_dates)

# create date matrix (not elegant, but it works); 20 = years back for reforecasts
date_mat <- matrix(data = NA, nrow = n_dates, ncol = 20)
# dates are not necessarily mondays and tuesdays anymore; [-1] to kick out the first reforecast day, i.e. the reference day; rev to get date structure as in netcdf
for (i in 1:n_dates) {
  date_mat[i, ] <- rev(as.character(seq(reforecast_dates[i], 
                                        by = "-1 year", 
                                        length.out = 21))[-1])
}


# with t(t2m_obs[k, , , i]) the observation matrix is of the same structure as the matrix date_mat; the same applies to the forecasts
benchmark_t2m_train <- foreach(k = 1:length(lead_times), .combine = rbind) %dopar% {
  
  # control dataset construction progress in log file report.txt
  sink("report.txt", append = TRUE) 
  cat("leadtime", k, "started!\n")
  
  # initialize data frame for all stations
  df_ids <- data.frame()
  
  for (i in 1:length(id)) {
    
    # first member is control forecast; remaining 10 forecasts are perturbed; data structure matches to matrix date_mat
    ens_fc <- sapply(1:length(ens_names), function(j) as.vector(t(t2m_fc[j, k, , , i])))
    colnames(ens_fc) <- ens_names
    
    df_id <- data.frame(init_date = as.POSIXct(as.vector(date_mat), tz = "UTC"),
                        valid_date = as.POSIXct(as.vector(date_mat), tz = "UTC") + 3600*lead_times[k],
                        leadtime = lead_times[k],
                        id = id[i],
                        name = station_name[i],
                        station_altitude = station_altitude[i],
                        station_longitude = station_longitude[i], 
                        station_latitude = station_latitude[i],
                        station_land_usage = station_land_usage[i],
                        model_altitude = model_altitude[i],
                        model_longitude = model_longitude[i], 
                        model_latitude = model_latitude[i],
                        model_land_usage = model_land_usage[i],
                        observation = as.vector(t(t2m_obs[k, , , i])),
                        ens_fc
                        )
    
    # add the data of all stations
    df_ids <- rbind(df_ids, df_id)
  }
  
  # add the data of all stations for one lead time
  df_ids

}

# set NA values to NaN
benchmark_t2m_train[is.na(benchmark_t2m_train)] <- NaN
# save data
save(benchmark_t2m_train, file = "benchmark_t2m_train.Rdata")

################################################################################

# PREPARE TEST DATA

# read data
data_fc <- nc_open("ESSD_benchmark_test_data_forecasts.nc")
data_obs <- nc_open("ESSD_benchmark_test_data_observations.nc")

# observations and forecasts of t2m
t2m_obs <- ncvar_get(data_obs, "t2m")
t2m_fc <- ncvar_get(data_fc, "t2m")

# station information
id <- as.character(data_fc[["dim"]][["station_id"]][["vals"]])
station_name <- ncvar_get(data_fc, "station_name")
station_altitude <- ncvar_get(data_fc, "station_altitude")
station_latitude <- ncvar_get(data_fc, "station_latitude")
station_longitude <- ncvar_get(data_fc, "station_longitude")
station_land_usage <- ncvar_get(data_fc, "station_land_usage")

# model information
model_altitude <- ncvar_get(data_fc, "model_altitude")
model_latitude <- ncvar_get(data_fc, "model_latitude")
model_longitude <- ncvar_get(data_fc, "model_longitude")
model_land_usage <- ncvar_get(data_fc, "model_land_usage")
ens_names <- c("ctrl", paste("ens", 1:50, sep = ""))

# time information
lead_times <- data_fc[["dim"]][["step"]][["vals"]] 
dates <- as.POSIXct(data_obs[["dim"]][["time"]][["vals"]], origin = "1970-01-01", tz = "UTC")


benchmark_t2m_test <- foreach(k = 1:length(lead_times), .combine = rbind) %dopar% {
  
  # control dataset construction progress in log file report.txt
  sink("report.txt", append = TRUE) 
  cat("leadtime", k, "started!\n")
  
  # initialize data frame for all stations
  df_ids <- data.frame()
  
  for (i in 1:length(id)) {
    
    # first member is control forecast; remaining 50 forecasts are perturbed
    ens_fc <- t(t2m_fc[, k, , i])
    colnames(ens_fc) <- ens_names
    
    df_id <- data.frame(init_date = dates,
                        valid_date = dates + 3600*lead_times[k],
                        leadtime = lead_times[k],
                        id = id[i],
                        name = station_name[i],
                        station_altitude = station_altitude[i],
                        station_longitude = station_longitude[i], 
                        station_latitude = station_latitude[i],
                        station_land_usage = station_land_usage[i],
                        model_altitude = model_altitude[i],
                        model_longitude = model_longitude[i], 
                        model_latitude = model_latitude[i],
                        model_land_usage = model_land_usage[i],
                        observation = t2m_obs[k, , i],
                        ens_fc
    )
    
    # add the data of all stations
    df_ids <- rbind(df_ids, df_id)
  }
  
  # add the data of all stations for one lead time
  df_ids
  
}

# set NA values to NaN
benchmark_t2m_test[is.na(benchmark_t2m_test)] <- NaN
# save data
save(benchmark_t2m_test, file = "benchmark_t2m_test.Rdata")


################################################################################

# EXTEND TRAINING DATASET WITH PREDICTORS

load("./benchmark_t2m_train.Rdata")
lead_times <- unique(benchmark_t2m_train$leadtime)
id <- unique(benchmark_t2m_train$id)
ens_names <- c("ctrl", paste("ens", 1:10, sep = ""))

benchmark_t2m_train_ext <- foreach(k = 1:length(lead_times), .combine = rbind) %dopar% {
  
  # control dataset construction progress in log file report.txt
  sink("report.txt", append = TRUE) 
  cat("leadtime", k, "started!\n")
  
  # initialize data frame for all stations
  df_ids <- data.frame()
  
  for (i in 1:length(id)) {
    
    loc <- benchmark_t2m_train[benchmark_t2m_train$leadtime == lead_times[k] & benchmark_t2m_train$id == id[i], ]
    # calculate additional predictors only from perturbed forecasts
    ens <- loc[, ens_names]
    m <- apply(ens[, -1], 1, mean, na.rm = TRUE)
    q <- t(apply(ens, 1, quantile, prob = c(0.1, 0.5, 0.9), type = 1, na.rm = TRUE))
    colnames(q) <- c("q10", "median", "q90")
    s <-  apply(ens, 1, sd, na.rm = TRUE)
    iqr <-  apply(ens, 1, IQR, type = 1, na.rm = TRUE)
    rand <-  apply(ens, 1, sample, size = 1)
    df_id <- data.frame(loc, 
                        mean = m,
                        q,
                        sd = s,
                        iqr = iqr, 
                        rand = rand)
    
    # add the data of all stations
    df_ids <- rbind(df_ids, df_id)
  }
  
  # add the data of all stations for one lead time
  df_ids
  
}

# set NA values to NaN
benchmark_t2m_train_ext[is.na(benchmark_t2m_train_ext)] <- NaN
# save data
save(benchmark_t2m_train_ext, file = "benchmark_t2m_train_ext.Rdata")


################################################################################

# EXTEND TEST DATASET WITH PREDICTORS

load("./benchmark_t2m_test.Rdata")
lead_times <- unique(benchmark_t2m_test$leadtime)
id <- unique(benchmark_t2m_test$id)
ens_names <- c("ctrl", paste("ens", 1:50, sep = ""))

benchmark_t2m_test_ext <- foreach(k = 1:length(lead_times), .combine = rbind) %dopar% {
  
  # control dataset construction progress in log file report.txt
  sink("report.txt", append = TRUE) 
  cat("leadtime", k, "started!\n")
  
  # initialize data frame for all stations
  df_ids <- data.frame()
  
  for (i in 1:length(id)) {
    
    loc <- benchmark_t2m_test[benchmark_t2m_test$leadtime == lead_times[k] & benchmark_t2m_test$id == id[i], ]
    # calculate additional predictors only from perturbed forecasts
    ens <- loc[, ens_names]
    m <- apply(ens[, -1], 1, mean, na.rm = TRUE)
    q <- t(apply(ens, 1, quantile, prob = c(0.1, 0.5, 0.9), type = 1, na.rm = TRUE))
    colnames(q) <- c("q10", "median", "q90")
    s <-  apply(ens, 1, sd, na.rm = TRUE)
    iqr <-  apply(ens, 1, IQR, type = 1, na.rm = TRUE)
    rand <-  apply(ens, 1, sample, size = 1)
    df_id <- data.frame(loc, 
                        mean = m,
                        q,
                        sd = s,
                        iqr = iqr, 
                        rand = rand)
    
    # add the data of all stations
    df_ids <- rbind(df_ids, df_id)
  }
  
  # add the data of all stations for one lead time
  df_ids
  
}

# set NA values to NaN
benchmark_t2m_test_ext[is.na(benchmark_t2m_test_ext)] <- NaN
# save data
save(benchmark_t2m_test_ext, file = "benchmark_t2m_test_ext.Rdata")
