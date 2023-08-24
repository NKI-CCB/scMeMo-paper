source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))
library(ncdf4)

dir.create("figures", showWarnings = F)

# We want to load the simulation results of the fitted model, but with a different likelihood file
# that has a more dense time sampling during mitosis in the standard simulation output
model <- bcm3.load.results(".", "../2-mitosis-model/output_gavetpines", prior_file="../2-mitosis-model/prior_gavetpines.xml", likelihood_file="likelihood_gavetpines_dense.xml", load_sampler_adaptation = F)

nthreads <- 4
temperature_ix <- dim(model$posterior$samples)[2]
map_sample_ix <- which.max(model$posterior$llikelihood[temperature_ix,] + model$posterior$lprior[temperature_ix,])
sample <- model$posterior$samples[,temperature_ix,map_sample_ix]
names(sample) <- model$variables

# Retrieve trajectories from the MAP estimate of the Gavet&Pines data fit
model <- bcm3.init.cpp(model, "", nthreads)
traj <- bcm3.cellpop.get.simulated.trajectories(model, model$likelihood$experiments[[1]]$name, sample)
species_names <- rep(NA, dim(traj$cells)[1])
for (i in 1:length(species_names)) {
  species_names[i] <- bcm3.cellpop.get.species.name(model, model$likelihood$experiments[[1]]$nam, i)
}
model <- bcm3.release.cpp(model)

plot(traj$time, traj$cells[2,,1], type='l', xlim=c(0, 2e4))
for (i in 2:32) {
  lines(traj$time, traj$cells[2,,i], type='l')
}

# Add some measurement noise
time_ix <- seq(100,400, by=5)
set.seed(100)
noise_sd <- 0.1
generated_data <- list()
generated_data[["CDK1"]] <- traj$cells[2,time_ix,]
generated_data[["CDK1"]] <- generated_data[["CDK1"]] + rnorm(length(generated_data[["CDK1"]]), 0, noise_sd)
generated_data[["CycB"]] <- traj$cells[1,time_ix,] + traj$cells[11,time_ix,]
generated_data[["CycB"]] <- generated_data[["CycB"]] + rnorm(length(generated_data[["CycB"]]), 0, noise_sd)

pdf("figures/generated_data.pdf", width=60/25.4, height=50/25.4)
par(cex=0.5)
plot(traj$time[time_ix]/60, generated_data[["CDK1"]][,1], type='l', xlim=c(0, 2e4/60), ylim=c(-0.5, 1.5), xlab="Time (minutes)", ylab="CDK1 sensor level")
for (i in 2:ncol(generated_data[["CDK1"]])) {
  lines(traj$time[time_ix]/60, generated_data[["CDK1"]][,i], type='l')
}
res <- dev.off()

dir.create("data")
dim1 <- ncdim_def("CDK1_sensor_sim/time_single_cells", "seconds", traj$time[time_ix])
dim2 <- ncdim_def("CDK1_sensor_sim/cells", "id", 1:ncol(generated_data[["CDK1"]]))
dim3 <- ncdim_def("CDK1_sensor_sim/cells_timing", "id", 1:32)
var1 <- ncvar_def("CDK1_sensor_sim/CDK1_sensor_sim", "signal", list(dim2, dim1), NA)
var2 <- ncvar_def("CDK1_sensor_sim/CycB_sim", "signal", list(dim2, dim1), NA)
var3 <- ncvar_def("CDK1_sensor_sim/NEBD_to_AO_duration", "seconds", dim3, NA)
ncnew <- nc_create("data/CDK1_sensor_simulated_data.nc", list(var1, var2, var3), force_v4=T)
for (i in 1:ncol(generated_data[["CDK1"]])) {
  ncvar_put(ncnew, var1, generated_data[["CDK1"]][,i], start=c(i,1), count=c(1, nrow(generated_data[["CDK1"]])))
}
for (i in 1:ncol(generated_data[["CycB"]])) {
  ncvar_put(ncnew, var2, generated_data[["CycB"]][,i], start=c(i,1), count=c(1, nrow(generated_data[["CycB"]])))
}
nc_close(ncnew)
diff(as.numeric(quantile(c(generated_data[["CDK1"]]), c(0.25, 0.75), na.rm=T)))
diff(as.numeric(quantile(c(generated_data[["CycB"]]), c(0.25, 0.75), na.rm=T)))




# Add timings for merged fit
timings20 <- list()

meraldi <- H5File$new("../2-mitosis-model/data/Meraldi2004.nc", 'r')
timings20[["meraldi"]] <- meraldi[["meraldi_fig1b/NEBD_to_AO_duration"]][]
meraldi$close_all()

zhou <- H5File$new("../2-mitosis-model/data/Zhou2017_20cells.nc", 'r')
timings20[["zhou"]] <- zhou[["Zhou_FigS1E/NEBD_to_AO_duration"]][]
zhou$close_all()

lu <- H5File$new("../2-mitosis-model/data/Lu2013_20cells.nc", 'r')
timings20[["lu"]] <- lu[["Lu_Fig5E/NEBD_to_AO_duration"]][]
lu$close_all()

liu <- H5File$new("../2-mitosis-model/data/Liu.nc", 'r')
timings20[["liu"]] <- liu[["Liu_fig4e/NEBD_to_AO_duration"]][]
liu$close_all()

dim4 <- ncdim_def("CDK1_sensor_sim/cells_timing_liu", "id", 1:length(timings20[["liu"]]))
dim5 <- ncdim_def("CDK1_sensor_sim/cells_timing_lu", "id", 1:length(timings20[["lu"]]))
dim6 <- ncdim_def("CDK1_sensor_sim/cells_timing_meraldi", "id", 1:length(timings20[["meraldi"]]))
dim7 <- ncdim_def("CDK1_sensor_sim/cells_timing_zhou", "id", 1:length(timings20[["zhou"]]))
var4 <- ncvar_def("CDK1_sensor_sim/NEBD_to_AO_duration_liu", "seconds", dim4, NA)
var5 <- ncvar_def("CDK1_sensor_sim/NEBD_to_AO_duration_lu", "seconds", dim5, NA)
var6 <- ncvar_def("CDK1_sensor_sim/NEBD_to_AO_duration_meraldi", "seconds", dim6, NA)
var7 <- ncvar_def("CDK1_sensor_sim/NEBD_to_AO_duration_zhou", "seconds", dim7, NA)
ncnew <- nc_create("data/CDK1_sensor_simulated_data_withtimings.nc", list(var1, var2, var4, var5, var6, var7), force_v4=T)
for (i in 1:ncol(generated_data[["CDK1"]])) {
  ncvar_put(ncnew, var1, generated_data[["CDK1"]][,i], start=c(i,1), count=c(1, nrow(generated_data[["CDK1"]])))
}
for (i in 1:ncol(generated_data[["CycB"]])) {
  ncvar_put(ncnew, var2, generated_data[["CycB"]][,i], start=c(i,1), count=c(1, nrow(generated_data[["CycB"]])))
}
ncvar_put(ncnew, var4, timings20[["liu"]])
ncvar_put(ncnew, var5, timings20[["lu"]])
ncvar_put(ncnew, var6, timings20[["meraldi"]])
ncvar_put(ncnew, var7, timings20[["zhou"]])
nc_close(ncnew)
