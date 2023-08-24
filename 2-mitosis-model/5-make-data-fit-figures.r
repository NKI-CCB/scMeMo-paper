source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))
library(pals)
source("../scripts/cellpop_functions.r")

gavetpines <- bcm3.load.results(".", "output_gavetpines", prior_file="prior_gavetpines.xml", likelihood_file="likelihood_gavetpines.xml")

temperature_ix <- dim(gavetpines$posterior$samples)[2]

gavetpines <- bcm3.init.cpp(gavetpines, threads=4)

pp <- get_posterior_predictive(gavetpines, 1, temperature_ix, numppdsamples=500)

pdf("figures/fit_gavetpines.pdf", width=100/25.4, height=100/25.4)
par(mfrow=c(2,2))
par(cex=0.5, mgp=c(2.0,0.8,0))
posterior_predictive_plot(gavetpines, pp, experiment_ix=1, data_ix=3, cell_ixs=1, use_time="minutes", xlim=c(-100,150), ylim=c(-0.7,1.2), ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(gavetpines, pp, experiment_ix=1, data_ix=3, cell_ixs=2, use_time="minutes", xlim=c(-100,150), ylim=c(-0.7,1.2), ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(gavetpines, pp, experiment_ix=1, data_ix=1, cell_ixs=1, use_time="minutes", xlim=c(-100,150), ylim=c(-0.7,1.2), ppd_point_symbol=20, ppd_color = "#00beff");
posterior_predictive_plot(gavetpines, pp, experiment_ix=1, data_ix=2, cell_ixs=1, use_time="minutes", xlim=c(-100,150), ylim=c(-0.7,1.2), ppd_point_symbol=20, ppd_color = "#00beff");
par(mfrow=c(1,1))
res <- dev.off()

gavetpines <- bcm3.release.cpp(gavetpines)

akopyan <- bcm3.load.results(".", "output_akopyan", prior_file="prior_akopyan.xml", likelihood_file="likelihood_akopyan.xml")

akopyan <- bcm3.init.cpp(akopyan, threads=4)

temperature_ix <- dim(akopyan$posterior$samples)[2]

pp <- get_posterior_predictive(akopyan, 1, temperature_ix, numppdsamples=500)

pdf("figures/fit_akopyan.pdf", width=70/25.4, height=50/25.4)
par(cex=0.5, mgp=c(2.0,0.8,0))
posterior_predictive_plot(akopyan, pp, experiment_ix=1, data_ix=1, cell_ixs=1, use_time="minutes", xlim=c(-1000, 320), ylim=c(-0.2,1.2), ppd_point_symbol=20, ppd_color = "#00beff");
res <- dev.off()


experiment_ix <- 1
plot_sample_ix <- which.max(akopyan$posterior$lposterior[temperature_ix,])
sample <- akopyan$posterior$samples[,temperature_ix,plot_sample_ix]

traj <- bcm3.cellpop.get.simulated.trajectories(akopyan, akopyan$likelihood$experiments[[experiment_ix]]$name, sample)

ix <- which(dimnames(traj$cells)[[1]] == "CyclinB" | dimnames(traj$cells)[[1]] == "active_CyclinB_CDK1")

pdf("figures/cell_trajectories.pdf", width=100/25.4, height=50/25.4)
par(mfrow=c(1,2))
par(cex=0.5, mgp=c(2.0,0.8,0))
plot(traj$time/60, apply(traj$cells[ix,,1], 2, sum), type='n', ylim=c(0, 0.9), xlab='Time (minutes)', ylab='Cyclin B', xlim=c(0,1600))
for (j in 1:dim(traj$cells)[3]) {
  lines(traj$time/60, apply(traj$cells[ix,,j], 2, sum))
}

nuclear_envelope_species_x <- which(dimnames(traj$cells)[[1]] == "nuclear_envelope")

plot(traj$time/60, apply(traj$cells[ix,,1], 2, sum), type='n', ylim=c(0, 0.9), xlab='Time relative to NEBD (minutes)', ylab='Cyclin B', xlim=c(-1000,200))
for (j in 1:dim(traj$cells)[3]) {
  NEBD_time <- traj$time[head(which(traj$cells[nuclear_envelope_species_x,,j] < 0.5), n=1)]
  if (length(NEBD_time) > 0) {
    lines((traj$time-NEBD_time)/60, apply(traj$cells[ix,,j], 2, sum))
  } 
}
res <- dev.off()

akopyan <- bcm3.release.cpp(akopyan)


pdf("figures/posterior_densities.pdf", width=230/25.4, height=100/25.4)
par(mfcol=c(2,6))
par(cex=0.5, mgp=c(2.0,0.8,0))
plot_variable_distribution(gavetpines, var_name="spindle_assembly")
plot_variable_distribution(gavetpines, var_name="spindle_assembly_variability")
plot_variable_distribution(gavetpines, var_name="FOXM1_synthesis")
plot_variable_distribution(gavetpines, var_name="FOXM1_synthesis_variability", ylim=c(0,1))
plot_variable_distribution(gavetpines, var_name="mitotic_entry")
plot_variable_distribution(gavetpines, var_name="mitotic_entry_variability")
plot_variable_distribution(akopyan, var_name="spindle_assembly")
plot_variable_distribution(akopyan, var_name="spindle_assembly_variability", ylim=c(0,1))
plot_variable_distribution(akopyan, var_name="FOXM1_synthesis")
plot_variable_distribution(akopyan, var_name="FOXM1_synthesis_variability")
plot_variable_distribution(akopyan, var_name="mitotic_entry", ylim=c(0,1))
plot_variable_distribution(akopyan, var_name="mitotic_entry_variability", ylim=c(0,1))
res <- dev.off()
