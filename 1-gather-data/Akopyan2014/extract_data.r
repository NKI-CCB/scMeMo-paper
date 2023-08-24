source("../extract_data_functions.r")

cycA <- extract_data_from_SVG("Figure3D.svg", seq(-20, -5, by=5), seq(5, 0, by=-1), point_layer_name="CycA")
cycB <- extract_data_from_SVG("Figure3D_cycB.svg", seq(-20, 5, by=5), seq(8, 0, by=-2), point_layer_name="CycB")

#add_t <- seq(5,15, len=10)
#add_y <- rep(NA, length(add_t))
#cycA <- rbind(cycA, data.frame(x=add_t, y=add_y))
#cycB <- rbind(cycB, data.frame(x=add_t, y=add_y))

plot(cycA$x, cycA$y)
plot(cycB$x, cycB$y)

library(ncdf4)

# dims <- list(ncdim_def("Akopyan_fig3d_cycA/time", "seconds", cycA$x * 3600),
#              ncdim_def("Akopyan_fig3d_cycB/time", "seconds", cycB$x * 3600))
# vars <- list(ncvar_def("Akopyan_fig3d_cycA/cycA", "fluorescence", dims[[1]], NA),
#              ncvar_def("Akopyan_fig3d_cycB/cycB", "fluorescence", dims[[2]], NA))
# 
# ncnew <- nc_create("Akopyan2014.nc", vars, force_v4=T)
# ncvar_put(ncnew, vars[[1]], cycA$y / 5)
# ncvar_put(ncnew, vars[[2]], cycB$y / 8)
# nc_close(ncnew)

merged_t <- sort(union(cycA$x, cycB$x))
merged_cycA <- rep(NA, length(merged_t))
merged_cycA[match(cycA$x, merged_t)] <- cycA$y
merged_cycB <- rep(NA, length(merged_t))
merged_cycB[match(cycB$x, merged_t)] <- cycB$y

min_cycA <- min(merged_cycA,na.rm=T)
min_cycB <- min(merged_cycB,na.rm=T)
scale_cycA <- max(merged_cycA,na.rm=T)
scale_cycB <- max(merged_cycB,na.rm=T)

scaled_cycA <- 2*(merged_cycA-min_cycA)/(scale_cycA-min_cycA)
scaled_cycB <- (merged_cycB-min_cycB)/(scale_cycB-min_cycB)

plot(merged_t, scaled_cycA)
points(merged_t, scaled_cycB,col='red')

dims <- list(ncdim_def("Akopyan_fig3d/time", "seconds", merged_t * 3600),
             ncdim_def("Akopyan_fig3d/cells_timing", "id", 1:32))
vars <- list(ncvar_def("Akopyan_fig3d/cycA", "fluorescence", dims[[1]], NA),
             ncvar_def("Akopyan_fig3d/cycB", "fluorescence", dims[[1]], NA),
             ncvar_def("Akopyan_fig3d/NEBD_to_AO_duration", "seconds", dims[[2]], NA))

ncnew <- nc_create("Akopyan2014.nc", vars, force_v4=T)
ncvar_put(ncnew, vars[[1]], scaled_cycA)
ncvar_put(ncnew, vars[[2]], scaled_cycB)
nc_close(ncnew)

diff(as.numeric(quantile(scaled_cycA, c(0.25, 0.75), na.rm=T)))/sqrt(13)
diff(as.numeric(quantile(scaled_cycB, c(0.25, 0.75), na.rm=T)))/sqrt(17)
