source("../extract_data_functions.r")

timings <- extract_data_from_SVG("Figure1E.svg", c(0, 200, 400), c(0, 1), point_layer_name = "timing")

barplot(timings$x, horiz = T, xlim=c(0, 400))
barplot(timings$x)


library(ncdf4)
dim2 <- ncdim_def("Lu_Fig5E/cells", "id", 1:nrow(timings))
var3 <- ncvar_def("Lu_Fig5E/NEBD_to_AO_duration", "seconds", dim2, NA)
ncnew <- nc_create("Lu2013.nc", list(var3), force_v4=T)
ncvar_put(ncnew, var3, as.numeric(timings$x)*60, start=1, count=nrow(timings))
nc_close(ncnew)

use_ix <- sample(nrow(timings), 20)
barplot(sort(timings$x[use_ix]))

dim2 <- ncdim_def("Lu_Fig5E/cells", "id", 1:20)
var3 <- ncvar_def("Lu_Fig5E/NEBD_to_AO_duration", "seconds", dim2, NA)
ncnew <- nc_create("Lu2013_20cells.nc", list(var3), force_v4=T)
ncvar_put(ncnew, var3, sort(as.numeric(timings$x[use_ix]))*60, start=1, count=20)
nc_close(ncnew)


diff(as.numeric(quantile(timings$x, c(0.25, 0.75))))*60
log10(0.5 * 300 / (diff(as.numeric(quantile(timings$x, c(0.25, 0.75))))*60))
