source("../extract_data_functions.r")

entry <- extract_data_from_SVG("Figure1A.svg", seq(0, 105, by=15), seq(125, 85, by=-5), point_layer_name = "entry")
exit <- extract_data_from_SVG("Figure1A.svg", seq(0, 105, by=15), seq(125, 85, by=-5), point_layer_name = "exit")
nebd_time <- extract_data_from_SVG("Figure1A.svg", seq(0, 105, by=15), seq(125, 85, by=-5), point_layer_name = "NEBD")
ana_time <- extract_data_from_SVG("Figure1A.svg", seq(0, 105, by=15), seq(125, 85, by=-5), point_layer_name = "ana")
mutant <- extract_data_from_SVG("Figure1B.svg", seq(0, 320, by=20), seq(125, 85, by=-5), point_layer_name = "points")
rpe_entry <- extract_data_from_SVG("Figure1E.svg", seq(0, 225, by=15), seq(120, 85, by=-5), point_layer_name = "entry")
rpe_exit <- extract_data_from_SVG("Figure1E.svg", seq(0, 225, by=15), seq(120, 85, by=-5), point_layer_name = "exit")
rpe_nebd_time <- extract_data_from_SVG("Figure1E.svg", seq(0, 225, by=15), seq(120, 85, by=-5), point_layer_name = "NEBD")
ana_nebd_time <- extract_data_from_SVG("Figure1E.svg", seq(0, 225, by=15), seq(120, 85, by=-5), point_layer_name = "ana")

plot(entry$x, entry$y, xlim=c(0, 195), ylim=c(85,125))
points(exit$x, exit$y)

normalized_entry <- entry$y - 100
normalized_exit <- exit$y - 100
max_signal <- max(c(normalized_entry, normalized_exit))
normalized_entry <- normalized_entry / max_signal
normalized_exit <- normalized_exit / max_signal

#mutant_normalized <- (mutant$y - 100) / max_signal
#plot(mutant$x, mutant_normalized)

normalized_rpe_entry <- rpe_entry$y - 100
normalized_rpe_exit <- rpe_exit$y - 100
normalized_rpe_entry <- normalized_rpe_entry / max_signal
normalized_rpe_exit <- normalized_rpe_exit / max_signal


plot(entry$x, normalized_entry, xlim=c(0, 350), ylim=c(-0.5, 1))
points(exit$x, normalized_exit)
points(mutant$x, mutant_normalized, col='red')
points(rpe_entry$x, normalized_rpe_entry, col='green')
points(rpe_exit$x, normalized_rpe_exit, col='green')


hela_individual_1 <- extract_data_from_SVG("FigureS6.svg", seq(0, 225, by=15), seq(120, 70, by=-10), point_layer_name = "cellAsensor")
hela_individual_2 <- extract_data_from_SVG("FigureS6.svg", seq(0, 225, by=15), seq(120, 70, by=-10), point_layer_name = "cellBsensor")

# Align by NEBD instead of metaphase
hela_individual_1$x <- hela_individual_1$x - 105
hela_individual_2$x <- hela_individual_2$x - 60

timepoints <- unique(c(hela_individual_1$x, hela_individual_2$x))
df <- data.frame(time=timepoints, cellA=NA, cellB=NA)
df$cellA[match(hela_individual_1$x, df$time)] <- (hela_individual_1$y - 100) / 12
df$cellB[match(hela_individual_2$x, df$time)] <- (hela_individual_2$y - 100) / 12
df <- df[order(df$time),]

plot(df$time, df$cellA, pch=19, col='black', ylim=c(-0.5, 1.2), xlim=c(-40, 250))
points(df$time, df$cellB, pch=19, col='grey')


library(ncdf4)
dim1 <- ncdim_def("GavetPines_HeLa/time_single_cells", "seconds", (df$time) * 60)
dim2 <- ncdim_def("GavetPines_HeLa/time_avg_entry", "seconds", (entry$x-nebd_time$x) * 60)
dim3 <- ncdim_def("GavetPines_HeLa/time_avg_exit", "seconds", (exit$x-ana_time$x) * 60)
dim4 <- ncdim_def("GavetPines_HeLa/cells", "id", 1:2)
dim5 <- ncdim_def("GavetPines_HeLa/cells_timing", "id", 1:32)
var1 <- ncvar_def("GavetPines_HeLa/CDK1_sensor_mean_entry", "signal", dim2, NA)
var2 <- ncvar_def("GavetPines_HeLa/CDK1_sensor_mean_exit", "signal", dim3, NA)
var3 <- ncvar_def("GavetPines_HeLa/CDK1_sensor", "signal", list(dim4, dim1), NA)
var4 <- ncvar_def("GavetPines_HeLa/NEBD_to_AO_duration", "seconds", dim5, NA)
ncnew <- nc_create("GavetPines2010.nc", list(var1, var2, var3, var4), force_v4=T)
ncvar_put(ncnew, var1, normalized_entry)
ncvar_put(ncnew, var2, normalized_exit)
ncvar_put(ncnew, var3, df$cellA, start=c(1,1), count=c(1, nrow(df)))
ncvar_put(ncnew, var3, df$cellB, start=c(2,1), count=c(1, nrow(df)))
nc_close(ncnew)

#diff(as.numeric(quantile(c(normalized_entry), c(0.25, 0.75))))/sqrt(6)
#diff(as.numeric(quantile(c(normalized_exit), c(0.25, 0.75))))/sqrt(6)
diff(as.numeric(quantile(c(normalized_entry, normalized_exit), c(0.25, 0.75))))/sqrt(6)
diff(as.numeric(quantile(c(df$cellA, df$cellB), c(0.25, 0.75), na.rm=T)))
