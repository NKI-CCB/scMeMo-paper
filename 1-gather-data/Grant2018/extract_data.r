source("../extract_data_functions.r")

df <- extract_data_from_SVG("Figure4C.svg", seq(1, 4, by=3), seq(25, 0, by=-5), point_layer_name = "u2os")

plot(df$y)

timings <- sample(df$y, 15)
points(timings, pch=19, col='red')


diff(as.numeric(quantile(df$y, c(0.25, 0.75))))*60*60


library(ncdf4)
dim2 <- ncdim_def("Grant2018_Fig4C/cells", "id", 1:length(timings))
var3 <- ncvar_def("Grant2018_Fig4C/Sphase_duration_U2OS", "seconds", dim2, NA)
ncnew <- nc_create("Grant2018_15cells.nc", list(var3), force_v4=T)
ncvar_put(ncnew, var3, as.numeric(timings)*60*60, start=1, count=length(timings))
nc_close(ncnew)

# use_ix <- sample(nrow(timings), 20)
# barplot(sort(timings$x[use_ix]))
# 
# dim2 <- ncdim_def("Lu_Fig5E/cells", "id", 1:20)
# var3 <- ncvar_def("Lu_Fig5E/NEBD_to_AO_duration", "seconds", dim2, NA)
# ncnew <- nc_create("Lu2013_20cells.nc", list(var3), force_v4=T)
# ncvar_put(ncnew, var3, sort(as.numeric(timings$x[use_ix]))*60, start=1, count=20)
# nc_close(ncnew)


