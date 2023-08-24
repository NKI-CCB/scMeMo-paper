df <- read.csv("imagej_quantification.csv", stringsAsFactors = F, check.names = F)

time <- as.numeric(df[-1,1]) * 60 * 60
e2f1 <- as.numeric(df[-1,2]) / max(as.numeric(df[-1,2]))
e2f7 <- as.numeric(df[-1,3]) / max(as.numeric(df[-1,3]))

plot(time, e2f1, pch=19)
points(time, e2f7, col='red', pch=19)

diff(as.numeric(quantile(e2f1, c(0.25, 0.75), na.rm=T)))
diff(as.numeric(quantile(e2f7, c(0.25, 0.75), na.rm=T)))

library(ncdf4)

dim <- ncdim_def("westendorp_fig1a/time", "seconds", time)
vars <- list(ncvar_def("westendorp_fig1a/E2F1", "signal", dim, NA),
             ncvar_def("westendorp_fig1a/E2F7", "signal", dim, NA))

ncnew <- nc_create("Westendorp2012.nc", vars, force_v4=T)
ncvar_put(ncnew, vars[[1]], e2f1)
ncvar_put(ncnew, vars[[2]], e2f7)
nc_close(ncnew)
