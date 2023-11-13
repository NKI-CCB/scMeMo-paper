source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))
library(pals)
library(ComplexHeatmap)
library(circlize)
source("../scripts/cellpop_functions.r")

models <- list()
models[["akopyan"]]                 <- bcm3.load.results(".", "output_akopyan_t24n50e2_j12k16",                prior_file="prior_akopyan.xml",    likelihood_file="likelihood_akopyan.xml",                load_sampler_adaptation=F)
models[["barr2016"]]                <- bcm3.load.results(".", "output_barr2016_t24n50e2_j12k16",               prior_file="prior_barr2016.xml",   likelihood_file="likelihood_barr2016.xml",               load_sampler_adaptation=F)
models[["barr2017"]]                <- bcm3.load.results(".", "output_barr2017_t24n20e2_j12k16",               prior_file="prior_barr2017.xml",   likelihood_file="likelihood_barr2017.xml",               load_sampler_adaptation=F)
models[["cappell"]]                 <- bcm3.load.results(".", "output_cappell_t24n20e2_j12k16",                prior_file="prior_cappell.xml",    likelihood_file="likelihood_cappell.xml",                load_sampler_adaptation=F)
models[["eward"]]                   <- bcm3.load.results(".", "output_eward_t24n20e2_j12k16",                  prior_file="prior_eward.xml",      likelihood_file="likelihood_eward.xml",                  load_sampler_adaptation=F)
models[["eward_noCycArepression"]]  <- bcm3.load.results(".", "output_eward_noCycArepression_t24n20e2_j12k16", prior_file="prior_eward.xml",      likelihood_file="likelihood_eward_noCycArepression.xml", load_sampler_adaptation=F)
models[["eward_noE2F7"]]            <- bcm3.load.results(".", "output_eward_noE2F7_t24n20e2_j12k16",           prior_file="prior_eward.xml",      likelihood_file="likelihood_eward_noE2F7.xml",           load_sampler_adaptation=F)
models[["westendorp"]]              <- bcm3.load.results(".", "output_westendorp_t24n20e2_j12k16",             prior_file="prior_westendorp.xml", likelihood_file="likelihood_westendorp.xml",             load_sampler_adaptation=F)

numthreads <- 8

marginal_likelihood(models[["eward"]])
marginal_likelihood(models[["eward_noCycArepression"]])
marginal_likelihood(models[["eward_noE2F7"]])

exp(marginal_likelihood(models[["eward"]]) - marginal_likelihood(models[["eward_noCycArepression"]]))
exp(marginal_likelihood(models[["eward"]]) - marginal_likelihood(models[["eward_noE2F7"]]))

pps <- readRDS("pps.rds")

dir.create("figures")

pdf("figures/posterior_predictive_akopyan.pdf", width=80/25.4, height=40/25.4)
par(mfrow=c(1,2))
par(cex=0.5, mgp=c(2.0,0.8,0))
posterior_predictive_plot(models[["akopyan"]], pps[["akopyan"]], 1, data_ix=1, cell_ixs=1, ylim=c(-0.2,2.2), ppd_point_symbol=20, ppd_color = "#00beff");
posterior_predictive_plot(models[["akopyan"]], pps[["akopyan"]], 1, data_ix=2, cell_ixs=1, ylim=c(-0.2,1.2), ppd_point_symbol=20, ppd_color = "#00beff");
res <- dev.off()

pdf("figures/posterior_predictive_barr2016.pdf", width=120/25.4, height=120/25.4)
par(mfcol=c(3,3))
par(cex=0.5, mgp=c(2.0,0.8,0))
posterior_predictive_plot(models[["barr2016"]], pps[["barr2016"]], 1, data_ix=1, cell_ixs=c(1), ylim=c(-0.4,2.6), ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(models[["barr2016"]], pps[["barr2016"]], 1, data_ix=1, cell_ixs=c(2), ylim=c(-0.4,2.6), ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(models[["barr2016"]], pps[["barr2016"]], 1, data_ix=1, cell_ixs=c(3), ylim=c(-0.4,2.6), ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(models[["barr2016"]], pps[["barr2016"]], 1, data_ix=2, cell_ixs=c(1), ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(models[["barr2016"]], pps[["barr2016"]], 1, data_ix=2, cell_ixs=c(2), ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(models[["barr2016"]], pps[["barr2016"]], 1, data_ix=2, cell_ixs=c(3), ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(models[["barr2016"]], pps[["barr2016"]], 1, data_ix=3, cell_ixs=c(1), ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(models[["barr2016"]], pps[["barr2016"]], 1, data_ix=3, cell_ixs=c(2), ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(models[["barr2016"]], pps[["barr2016"]], 1, data_ix=3, cell_ixs=c(3), ppd_point_symbol=4, ppd_color = "#00b25d");
res <- dev.off()

pdf("figures/posterior_predictive_barr2017.pdf", width=80/25.4, height=80/25.4)
par(mfrow=c(2,2))
par(cex=0.5, mgp=c(2.0,0.8,0))
posterior_predictive_plot(models[["barr2017"]], pps[["barr2017"]], 1, data_ix=1, cell_ixs=c(4) , ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(models[["barr2017"]], pps[["barr2017"]], 1, data_ix=1, cell_ixs=c(16), ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(models[["barr2017"]], pps[["barr2017"]], 1, data_ix=1, cell_ixs=c(13), ppd_point_symbol=4, ppd_color = "#00b25d");
posterior_predictive_plot(models[["barr2017"]], pps[["barr2017"]], 1, data_ix=1, cell_ixs=c(1) , ppd_point_symbol=4, ppd_color = "#00b25d");
res <- dev.off()

pdf("figures/posterior_predictive_cappell.pdf", width=80/25.4, height=80/25.4)
par(mfrow=c(2,2))
par(cex=0.5, mgp=c(2.0,0.8,0))
posterior_predictive_plot_replicate_data(models[["cappell"]], pps[["cappell"]], 1, data_ix=1, cell_ix=1, ppd_point_symbol=20, ppd_color = "#00beff");
posterior_predictive_plot_replicate_data(models[["cappell"]], pps[["cappell"]], 1, data_ix=3, cell_ix=1, ppd_point_symbol=20, ppd_color = "#00beff");
posterior_predictive_plot_replicate_data(models[["cappell"]], pps[["cappell"]], 1, data_ix=5, cell_ix=1, ppd_point_symbol=20, ppd_color = "#00beff");
posterior_predictive_plot_replicate_data(models[["cappell"]], pps[["cappell"]], 1, data_ix=7, cell_ix=1, ylim=c(0.8,2.2), ppd_point_symbol=20, ppd_color = "#00beff");
res <- dev.off()

pdf("figures/posterior_predictive_eward.pdf", width=140/25.4, height=40/25.4)
par(mfrow=c(1,3))
par(cex=0.5, mgp=c(2.0,0.8,0))
posterior_predictive_plot_replicate_data(models[["eward"]], pps[["eward"]], 1, 1, 1, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward"]], pps[["eward"]], 1, 3, 1, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward"]], pps[["eward"]], 1, 5, 1, ppd_point_symbol=20, ppd_color = "#00beff")
res <- dev.off()

pdf("figures/posterior_predictive_eward_noCycArepression.pdf", width=140/25.4, height=40/25.4)
par(mfrow=c(1,3))
par(cex=0.5, mgp=c(2.0,0.8,0))
posterior_predictive_plot_replicate_data(models[["eward_noCycArepression"]], pps[["eward_noCycArepression"]], 1, 1, 1, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward_noCycArepression"]], pps[["eward_noCycArepression"]], 1, 3, 1, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward_noCycArepression"]], pps[["eward_noCycArepression"]], 1, 5, 1, ppd_point_symbol=20, ppd_color = "#00beff")
res <- dev.off()

pdf("figures/posterior_predictive_eward_noE2F7.pdf", width=140/25.4, height=40/25.4)
par(mfrow=c(1,3))
par(cex=0.5, mgp=c(2.0,0.8,0))
posterior_predictive_plot_replicate_data(models[["eward_noE2F7"]], pps[["eward_noE2F7"]], 1, 1, 1, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward_noE2F7"]], pps[["eward_noE2F7"]], 1, 3, 1, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward_noE2F7"]], pps[["eward_noE2F7"]], 1, 5, 1, ppd_point_symbol=20, ppd_color = "#00beff")
res <- dev.off()

pdf("figures/posterior_predictive_westendorp.pdf", width=80/25.4, height=40/25.4)
par(mfrow=c(1,2))
par(cex=0.5, mgp=c(2.0,0.8,0))
posterior_predictive_plot(models[["westendorp"]], pps[["westendorp"]], 1, data_ix=1, cell_ixs=1, ppd_point_symbol=20, ppd_color = "#00beff");
posterior_predictive_plot(models[["westendorp"]], pps[["westendorp"]], 1, data_ix=2, cell_ixs=1, ppd_point_symbol=20, ppd_color = "#00beff");
res <- dev.off()

models[["akopyan"]] <- bcm3.init.cpp(models[["akopyan"]], "", numthreads)
pdf("figures/trajectory_heatmaps_akopyan.pdf", width=5, height=5)
trajectory_heatmap(models[["akopyan"]], c("CyclinA", "active_CyclinA_CDK2", "active_CyclinA_CDK1"), range=c(0,2))
trajectory_heatmap(models[["akopyan"]], c("CyclinB", "active_CyclinB_CDK1"))
res <- dev.off()
models[["akopyan"]] <- bcm3.release.cpp(models[["akopyan"]])

models[["barr2016"]] <- bcm3.init.cpp(models[["barr2016"]], "", numthreads)
pdf("figures/trajectory_heatmaps_barr2016.pdf", width=5, height=5)
trajectory_heatmap(models[["barr2016"]], c("active_p27", "inactive_p27"))
trajectory_heatmap(models[["barr2016"]], c("CyclinE", "active_CyclinE_CDK2"))
trajectory_heatmap(models[["barr2016"]], c("CyclinA", "active_CyclinA_CDK2", "active_CyclinA_CDK1"), range=c(0,2))
trajectory_heatmap(models[["barr2016"]], c("licensed_DNA"), range=c(0,2))
trajectory_heatmap(models[["barr2016"]], c("replicating_DNA", "replicated_DNA"), range=c(0,2))
res <- dev.off()
models[["barr2016"]] <- bcm3.release.cpp(models[["barr2016"]])

models[["barr2017"]] <- bcm3.init.cpp(models[["barr2017"]], "", numthreads)
pdf("figures/trajectory_heatmaps_barr2017.pdf", width=5, height=5)
trajectory_heatmap(models[["barr2017"]], c("p21"))
trajectory_heatmap(models[["barr2017"]], c("replication_induced_DNA_damage"))
res <- dev.off()
models[["barr2017"]] <- bcm3.release.cpp(models[["barr2017"]])

models[["cappell"]] <- bcm3.init.cpp(models[["cappell"]], "", numthreads)
pdf("figures/trajectory_heatmaps_cappell.pdf", width=5, height=5)
trajectory_heatmap(models[["cappell"]], c("CyclinE", "active_CyclinE_CDK2"))
trajectory_heatmap(models[["cappell"]], c("Emi1"))
trajectory_heatmap(models[["cappell"]], c("Geminin"))
trajectory_heatmap(models[["cappell"]], c("licensed_DNA", "DNA", "replicating_DNA", "replicated_DNA"), range=c(0,2))
res <- dev.off()
models[["cappell"]] <- bcm3.release.cpp(models[["cappell"]])

models[["westendorp"]] <- bcm3.init.cpp(models[["westendorp"]], "", numthreads)
pdf("figures/trajectory_heatmaps_westendorp.pdf", width=5, height=5)
trajectory_heatmap(models[["westendorp"]], c("active_E2F1", "inactive_E2F1"))
trajectory_heatmap(models[["westendorp"]], c("active_E2F7"))
trajectory_heatmap(models[["westendorp"]], c("replicating_DNA", "replicated_DNA"), range=c(0,2))
res <- dev.off()
models[["westendorp"]] <- bcm3.release.cpp(models[["westendorp"]])


pps[["eward"]] <- calculate_Rsq_replicated_data(models[["eward"]], pps[["eward"]], 1)
pps[["eward_noCycArepression"]] <- calculate_Rsq_replicated_data(models[["eward_noCycArepression"]], pps[["eward_noCycArepression"]], 1)
pps[["eward_noE2F7"]] <- calculate_Rsq_replicated_data(models[["eward_noE2F7"]], pps[["eward_noE2F7"]], 1)

sapply(pps[["eward"]]$Rsq, mean)
sapply(pps[["eward"]]$Rsq, quantile, probs=c(0.05, 0.95))
eward_title1 <- paste("R2=", format(mean(pps[["eward"]]$Rsq[[1]]), digits=2), " [", format(quantile(pps[["eward"]]$Rsq[[1]], 0.05), digits=2), "-", format(quantile(pps[["eward"]]$Rsq[[1]], 0.95), digits=2), "]", sep="")
eward_title2 <- paste("R2=", format(mean(pps[["eward"]]$Rsq[[2]]), digits=2), " [", format(quantile(pps[["eward"]]$Rsq[[2]], 0.05), digits=2), "-", format(quantile(pps[["eward"]]$Rsq[[2]], 0.95), digits=2), "]", sep="")
eward_title3 <- paste("R2=", format(mean(pps[["eward"]]$Rsq[[3]]), digits=2), " [", format(quantile(pps[["eward"]]$Rsq[[3]], 0.05), digits=2), "-", format(quantile(pps[["eward"]]$Rsq[[3]], 0.95), digits=2), "]", sep="")

sapply(pps[["eward_noCycArepression"]]$Rsq, mean)
sapply(pps[["eward_noCycArepression"]]$Rsq, quantile, probs=c(0.05, 0.95))
eward_noCycArepression_title1 <- paste("R2=", format(mean(pps[["eward_noCycArepression"]]$Rsq[[1]]), digits=2), " [", format(quantile(pps[["eward_noCycArepression"]]$Rsq[[1]], 0.05), digits=2), "-", format(quantile(pps[["eward_noCycArepression"]]$Rsq[[1]], 0.95), digits=2), "]", sep="")
eward_noCycArepression_title2 <- paste("R2=", format(mean(pps[["eward_noCycArepression"]]$Rsq[[2]]), digits=2), " [", format(quantile(pps[["eward_noCycArepression"]]$Rsq[[2]], 0.05), digits=2), "-", format(quantile(pps[["eward_noCycArepression"]]$Rsq[[2]], 0.95), digits=2), "]", sep="")
eward_noCycArepression_title3 <- paste("R2=", format(mean(pps[["eward_noCycArepression"]]$Rsq[[3]]), digits=2), " [", format(quantile(pps[["eward_noCycArepression"]]$Rsq[[3]], 0.05), digits=2), "-", format(quantile(pps[["eward_noCycArepression"]]$Rsq[[3]], 0.95), digits=2), "]", sep="")

sapply(pps[["eward_noE2F7"]]$Rsq, mean)
sapply(pps[["eward_noE2F7"]]$Rsq, quantile, probs=c(0.05, 0.95))
eward_noE2F7_title1 <- paste("R2=", format(mean(pps[["eward_noE2F7"]]$Rsq[[1]]), digits=2), " [", format(quantile(pps[["eward_noE2F7"]]$Rsq[[1]], 0.05), digits=2), "-", format(quantile(pps[["eward_noE2F7"]]$Rsq[[1]], 0.95), digits=2), "]", sep="")
eward_noE2F7_title2 <- paste("R2=", format(mean(pps[["eward_noE2F7"]]$Rsq[[2]]), digits=2), " [", format(quantile(pps[["eward_noE2F7"]]$Rsq[[2]], 0.05), digits=2), "-", format(quantile(pps[["eward_noE2F7"]]$Rsq[[2]], 0.95), digits=2), "]", sep="")
eward_noE2F7_title3 <- paste("R2=", format(mean(pps[["eward_noE2F7"]]$Rsq[[3]]), digits=2), " [", format(quantile(pps[["eward_noE2F7"]]$Rsq[[3]], 0.05), digits=2), "-", format(quantile(pps[["eward_noE2F7"]]$Rsq[[3]], 0.95), digits=2), "]", sep="")

pdf("figures/posterior_predictive_eward_model_comparison.pdf", width=210/25.4, height=210/25.4)
par(mfcol=c(3,3))
posterior_predictive_plot_replicate_data(models[["eward"]], pps[["eward"]], 1, 1, 1, main=eward_title1, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward"]], pps[["eward"]], 1, 3, 1, main=eward_title2, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward"]], pps[["eward"]], 1, 5, 1, main=eward_title3, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward_noCycArepression"]], pps[["eward_noCycArepression"]], 1, 1, 1, main=eward_noCycArepression_title1, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward_noCycArepression"]], pps[["eward_noCycArepression"]], 1, 3, 1, main=eward_noCycArepression_title2, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward_noCycArepression"]], pps[["eward_noCycArepression"]], 1, 5, 1, main=eward_noCycArepression_title3, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward_noE2F7"]], pps[["eward_noE2F7"]], 1, 1, 1, main=eward_noE2F7_title1, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward_noE2F7"]], pps[["eward_noE2F7"]], 1, 3, 1, main=eward_noE2F7_title2, ppd_point_symbol=20, ppd_color = "#00beff")
posterior_predictive_plot_replicate_data(models[["eward_noE2F7"]], pps[["eward_noE2F7"]], 1, 5, 1, main=eward_noE2F7_title3, ppd_point_symbol=20, ppd_color = "#00beff")
res <- dev.off()

temperature_ix <- dim(models[["eward"]]$posterior$samples)[2]
MAP_sample_ix <- which.max(models[["eward"]]$posterior$lposterior[temperature_ix,])
MAP_sample <- models[["eward"]]$posterior$samples[,temperature_ix,MAP_sample_ix]

models[["eward"]] <- bcm3.init.cpp(models[["eward"]], "", numthreads)
pdf("figures/trajectory_heatmaps_eward.pdf", width=5, height=5)
trajectory_heatmap(models[["eward"]], c("CyclinE_mRNA"))
trajectory_heatmap(models[["eward"]], c("CyclinA_mRNA"))
trajectory_heatmap(models[["eward"]], c("CyclinB_mRNA"))
dev.off()

png("figures/trajectory_heatmaps_eward_cycE.png", width=710, height=710)
trajectory_heatmap(models[["eward"]], c("CyclinE_mRNA"), draw_whole_population = T, draw_population_average = F, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_cycA.png", width=710, height=710)
trajectory_heatmap(models[["eward"]], c("CyclinA_mRNA"), draw_whole_population = T, draw_population_average = F, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_cycB.png", width=710, height=710)
trajectory_heatmap(models[["eward"]], c("CyclinB_mRNA"), draw_whole_population = T, draw_population_average = F, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_cycE_popavg.png", width=710, height=60)
trajectory_heatmap(models[["eward"]], c("CyclinE_mRNA"), draw_whole_population = F, draw_population_average = T, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_cycA_popavg.png", width=710, height=60)
trajectory_heatmap(models[["eward"]], c("CyclinA_mRNA"), draw_whole_population = F, draw_population_average = T, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_cycB_popavg.png", width=710, height=60)
trajectory_heatmap(models[["eward"]], c("CyclinB_mRNA"), draw_whole_population = F, draw_population_average = T, show_heatmap_legend = F)
dev.off()
models[["eward"]] <- bcm3.release.cpp(models[["eward"]])

models[["eward_noCycArepression"]] <- bcm3.init.cpp(models[["eward_noCycArepression"]], "", numthreads)
pdf("figures/trajectory_heatmaps_eward_noCycArepression.pdf", width=5, height=5)
trajectory_heatmap(models[["eward_noCycArepression"]], c("CyclinE_mRNA"))
trajectory_heatmap(models[["eward_noCycArepression"]], c("CyclinA_mRNA"))
trajectory_heatmap(models[["eward_noCycArepression"]], c("CyclinB_mRNA"))
dev.off()

png("figures/trajectory_heatmaps_eward_noCycArepression_cycE.png", width=710, height=710)
trajectory_heatmap(models[["eward_noCycArepression"]], c("CyclinE_mRNA"), draw_whole_population = T, draw_population_average = F, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_noCycArepression_cycA.png", width=710, height=710)
trajectory_heatmap(models[["eward_noCycArepression"]], c("CyclinA_mRNA"), draw_whole_population = T, draw_population_average = F, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_noCycArepression_cycB.png", width=710, height=710)
trajectory_heatmap(models[["eward_noCycArepression"]], c("CyclinB_mRNA"), draw_whole_population = T, draw_population_average = F, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_noCycArepression_cycE_popavg.png", width=710, height=60)
trajectory_heatmap(models[["eward_noCycArepression"]], c("CyclinE_mRNA"), draw_whole_population = F, draw_population_average = T, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_noCycArepression_cycA_popavg.png", width=710, height=60)
trajectory_heatmap(models[["eward_noCycArepression"]], c("CyclinA_mRNA"), draw_whole_population = F, draw_population_average = T, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_noCycArepression_cycB_popavg.png", width=710, height=60)
trajectory_heatmap(models[["eward_noCycArepression"]], c("CyclinB_mRNA"), draw_whole_population = F, draw_population_average = T, show_heatmap_legend = F)
dev.off()
models[["eward_noCycArepression"]] <- bcm3.release.cpp(models[["eward_noCycArepression"]])

models[["eward_noE2F7"]] <- bcm3.init.cpp(models[["eward_noE2F7"]], "", numthreads)
pdf("figures/trajectory_heatmaps_eward_noE2F7d.pdf", width=5, height=5)
trajectory_heatmap(models[["eward_noE2F7"]], c("CyclinE_mRNA"))
trajectory_heatmap(models[["eward_noE2F7"]], c("CyclinA_mRNA"))
trajectory_heatmap(models[["eward_noE2F7"]], c("CyclinB_mRNA"))
dev.off()

png("figures/trajectory_heatmaps_eward_noE2F7_cycE.png", width=710, height=710)
trajectory_heatmap(models[["eward_noE2F7"]], c("CyclinE_mRNA"), draw_whole_population = T, draw_population_average = F, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_noE2F7_cycA.png", width=710, height=710)
trajectory_heatmap(models[["eward_noE2F7"]], c("CyclinA_mRNA"), draw_whole_population = T, draw_population_average = F, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_noE2F7_cycB.png", width=710, height=710)
trajectory_heatmap(models[["eward_noE2F7"]], c("CyclinB_mRNA"), draw_whole_population = T, draw_population_average = F, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_noE2F7_cycE_popavg.png", width=710, height=60)
trajectory_heatmap(models[["eward_noE2F7"]], c("CyclinE_mRNA"), draw_whole_population = F, draw_population_average = T, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_noE2F7_cycA_popavg.png", width=710, height=60)
trajectory_heatmap(models[["eward_noE2F7"]], c("CyclinA_mRNA"), draw_whole_population = F, draw_population_average = T, show_heatmap_legend = F)
dev.off()
png("figures/trajectory_heatmaps_eward_noE2F7_cycB_popavg.png", width=710, height=60)
trajectory_heatmap(models[["eward_noE2F7"]], c("CyclinB_mRNA"), draw_whole_population = F, draw_population_average = T, show_heatmap_legend = F)
dev.off()
models[["eward_noE2F7"]] <- bcm3.release.cpp(models[["eward_noE2F7"]])

for (j in 1:length(models)) {
  pdf(paste("figures/posterior_densities_", names(models)[j], "_all.pdf", sep=""), width=200/25.4, height=160/25.4)
  par(mfrow=c(4,5))
  par(cex=0.5, mgp=c(2.0,0.8,0), mar=c(3,3,3,1))
  for (i in 1:models[[j]]$nvar) {
    plot_variable_distribution(models[[j]], var_ix=i)
  }
  res <- dev.off()
}