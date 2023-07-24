library(ncdf4)
# for faster construction! 
library(doMC)
registerDoMC(cores = detectCores()-1)

# create empty array for saving data
dvine.pp <- array(NaN, dim = c(51, 21, 730, 229))
lt <- seq(0, 120, by = 6)

# get quantiles for each lead time
out <- foreach(l = 1:21) %dopar% {
  load(paste0("results.dvine.lt", lt[l], ".nw2.Rdata"))
  stat.quantiles <- array(NaN, dim = c(51, 730, 229))
  for (s in 1:229) {
    stat.quantiles[, , s] <- t(as.matrix(results.dvine[[s]]$vquantiles[, -c(1:2)]))
  }
  stat.quantiles
}

# save quantiles for each lead time
for (l in 1:21) {
  dvine.pp[, l, , ] <- out[[l]]
}


# write TRUE is important, as otherwise file is not writeable!
data <- nc_open("1_ESSD-benchmark_University-of-Hildesheim_D-Vine-Copula_v1.0.nc", write = TRUE)

# get global attributes of t2m
ncatt_get(data, "t2m")

# redefine global attributes of t2m
ncatt_put(data, 
          varid = "t2m", 
          attname = "institution",
          attval = "University_of_Hildesheim")
ncatt_put(data, 
          varid = "t2m", 
          attname = "model",
          attval = "D-Vine-Copula")
ncatt_put(data, 
          varid = "t2m", 
          attname = "version",
          attval = "v.1.0")
ncatt_put(data, 
          varid = "t2m", 
          attname = "output",
          attval = "quantiles")
ncatt_get(data, "t2m")

# provide quantiles in the same structure, as the netcdf file, i.e.
# 50 x 21 x 730 x 229
x <- as.vector(dvine.pp)

# overwrite value as vector (order automatically correct, if structure is the same as in original file)
ncvar_put(data, varid = "t2m", vals = x)

# saves changes on the netcdf file
nc_sync(data)

