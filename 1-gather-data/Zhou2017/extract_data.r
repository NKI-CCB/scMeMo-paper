source("../extract_data_functions.r")

timings <- extract_data_from_SVG("FigureS1E.svg", seq(120, 0, -20), c(0, 1), point_layer_name = "timing")

barplot(timings$x)

library(ncdf4)
dim2 <- ncdim_def("Zhou_FigS1E/cells", "id", 1:nrow(timings))
var3 <- ncvar_def("Zhou_FigS1E/NEBD_to_AO_duration", "seconds", dim2, NA)
ncnew <- nc_create("Zhou2017.nc", list(var3), force_v4=T)
ncvar_put(ncnew, var3, as.numeric(timings$x)*60, start=1, count=nrow(timings))
nc_close(ncnew)

use_ix <- sample(nrow(timings), 20)
barplot(sort(timings$x[use_ix]))

dim2 <- ncdim_def("Zhou_FigS1E/cells", "id", 1:20)
var3 <- ncvar_def("Zhou_FigS1E/NEBD_to_AO_duration", "seconds", dim2, NA)
ncnew <- nc_create("Zhou2017_20cells.nc", list(var3), force_v4=T)
ncvar_put(ncnew, var3, sort(as.numeric(timings$x[use_ix]))*60, start=1, count=20)
nc_close(ncnew)

diff(as.numeric(quantile(timings$x, c(0.25, 0.75))))*60
