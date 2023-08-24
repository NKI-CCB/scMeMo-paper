source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))
library(pals)
library(RcppHungarian)
library(beeswarm)

model <- bcm3.load.results(".", "../2-mitosis-model/output_gavetpines", prior_file="../2-mitosis-model/prior_gavetpines.xml", likelihood_file="likelihood_gavetpines_dense.xml", load_sampler_adaptation = F)
temperature_ix <- dim(model$posterior$samples)[2]
map_sample_ix <- which.max(model$posterior$llikelihood[temperature_ix,] + model$posterior$lprior[temperature_ix,])
sample <- model$posterior$samples[,temperature_ix,map_sample_ix]
names(sample) <- model$variables

dir.create("figures", showWarnings = F)

correct_fit <- bcm3.load.results(".", "output_simdata_correctspecified_t24_n100", prior_file="prior_simulated_data.xml", likelihood_file="likelihood_simulated_data.xml")
misspec_fit <- bcm3.load.results(".", "output_simdata_misspecified_t24_n100", prior_file="prior_simulated_data_misspecified.xml", likelihood_file="likelihood_simulated_data_misspecified.xml")
misspec_fit_with_timings <- bcm3.load.results(".", "output_simdata_misspecified_with_mitosis_timings_t24_n100", prior_file="prior_simulated_data_misspecified.xml", likelihood_file="likelihood_simulated_data_misspecified_with_mitosis_timings.xml")

logloss <- readRDS("../3-mitosis-model-test-prediction/logloss.rds")
summed <- t(logloss[[1]]$logloss) + t(logloss[[2]]$logloss) + t(logloss[[3]]$logloss) + t(logloss[[4]]$logloss)
Rlogloss_log <- apply(log10(t(summed)), 1, mean)
calibrated_temperature_ix <- which.min(Rlogloss_log)

# Parameter densities

pdf("figures/all_densities_correctspecified.pdf", width=160/(25.4*0.75), height=80/(25.4*0.75))
par(mfrow=c(3,5))
par(mar=c(3,2,2,1))
par(cex=0.5)
for (i in 1:15) {
  plot_variable_distribution(correct_fit, var_ix=i)
  ix <- match(correct_fit$variables[i], names(sample))
  if (!is.na(ix)) {
    abline(v = sample[ix], col='black', lwd=2)
  }
}
par(mfrow=c(1,1))
res <- dev.off()

pdf("figures/all_densities_misspecified.pdf", width=160/(25.4*0.75), height=80/(25.4*0.75))
par(mfrow=c(3,5))
par(mar=c(3,2,2,1))
par(cex=0.5)
for (i in 1:12) {
  plot_variable_distribution(misspec_fit, var_ix=i)
  ix <- match(misspec_fit$variables[i], names(sample))
  if (!is.na(ix)) {
    abline(v = sample[ix], col='black', lwd=2)
  }
}
par(mfrow=c(1,1))
res <- dev.off()

pdf("figures/all_densities_misspecified_calibrated.pdf", width=160/(25.4*0.75), height=80/(25.4*0.75))
par(mfrow=c(3,5))
par(mar=c(3,2,2,1))
par(cex=0.5)
for (i in 1:12) {
  plot_variable_distribution(misspec_fit, var_ix=i, temperature_ix=calibrated_temperature_ix)
  ix <- match(misspec_fit$variables[i], names(sample))
  if (!is.na(ix)) {
    abline(v = sample[ix], col='black', lwd=2)
  }
}
par(mfrow=c(1,1))
res <- dev.off()

pdf("figures/all_densities_misspecified_fithwithtimings.pdf", width=160/(25.4*0.75), height=80/(25.4*0.75))
par(mfrow=c(3,5))
par(mar=c(3,2,2,1))
par(cex=0.5)
for (i in 1:12) {
  plot_variable_distribution(misspec_fit_with_timings, var_ix=i)
  ix <- match(misspec_fit_with_timings$variables[i], names(sample))
  if (!is.na(ix)) {
    abline(v = sample[ix], col='black', lwd=2)
  }
}
par(mfrow=c(1,1))
res <- dev.off()

pdf("figures/densities_comparison.pdf", width=100/(25.4*0.75), height=80/(25.4*0.75))
par(mfcol=c(4,4))
par(mar=c(3,2,2,1))
par(cex=0.5)
for (i in c(4,6,5,3)) {
  plot_variable_distribution(correct_fit, var_ix=i)
  ix <- match(correct_fit$variables[i], names(sample))
  if (!is.na(ix)) {
    abline(v = sample[ix], col='black', lwd=2)
  }
  
  plot_variable_distribution(misspec_fit, var_ix=i)
  ix <- match(misspec_fit$variables[i], names(sample))
  if (!is.na(ix)) {
    abline(v = sample[ix], col='black', lwd=2)
  }
  
  plot_variable_distribution(misspec_fit, var_ix=i, temperature_ix=calibrated_temperature_ix)
  ix <- match(misspec_fit$variables[i], names(sample))
  if (!is.na(ix)) {
    abline(v = sample[ix], col='black', lwd=2)
  }
  
  plot_variable_distribution(misspec_fit_with_timings, var_ix=i)
  ix <- match(misspec_fit_with_timings$variables[i], names(sample))
  if (!is.na(ix)) {
    abline(v = sample[ix], col='black', lwd=2)
  }
}
res <- dev.off()

# Posterior data predictive

nthreads <- 4
numppdsamples <- 200
ppd_sample_ix <- sample((dim(correct_fit$posterior$samples)[3]/2+1):dim(correct_fit$posterior$samples)[3], numppdsamples)

correct_fit <- bcm3.init.cpp(correct_fit, "", nthreads)
misspec_fit <- bcm3.init.cpp(misspec_fit, "", nthreads)

pdf("figures/ppd_comparison.pdf", width=120/(25.4*0.75), height=120/(25.4*0.75))
par(mfrow=c(3,3))

experiment_ix <- 1
data_ix <- 1
ppd_color <- "#00b25d"
ppd_symbol <- 4
ylim <- c(-0.4,1.4)

obsdata <- bcm3.cellpop.get.observed.data(correct_fit, correct_fit$likelihood$experiments[[experiment_ix]]$name)
merged <- list()
for (j in 1:length(obsdata)) {
  merged[[j]] <- list()
  for (l in 1:dim(obsdata[[j]]$data)[2]) {
    merged[[j]][[l]] <- matrix(NA, numppdsamples, length(obsdata[[j]]$time))
  }
}
for (i in 1:numppdsamples) {
  sample <- correct_fit$posterior$samples[,temperature_ix,ppd_sample_ix[i]]
  simdata <- bcm3.cellpop.get.simulated.data(correct_fit, correct_fit$likelihood$experiments[[experiment_ix]]$name, sample)
  for (j in 1:length(obsdata)) {
    for (l in 1:dim(obsdata[[j]]$data)[2]) {
      merged[[j]][[l]][i,] <- simdata[[j]]$data[,l,1]
    }
  }
}
stdev_ix <- match(correct_fit$likelihood$experiments[[experiment_ix]]$data[[data_ix]]$stdev, correct_fit$variables)
sds <- 10^as.numeric(correct_fit$posterior$samples[stdev_ix, temperature_ix, ppd_sample_ix])
t <- obsdata[[data_ix]]$time / 60
for (cell_ix in c(2,8,14)) {
  lower <- rep(NA, length(obsdata[[data_ix]]$time))
  upper <- rep(NA, length(obsdata[[data_ix]]$time))
  for (i in 1:length(obsdata[[data_ix]]$time)) {
    x <- rnorm(1000, merged[[data_ix]][[cell_ix]][,i], sds * correct_fit$likelihood$experiments[[experiment_ix]]$data[[data_ix]]$stdev_multiplication_factor)
    lower[i] <- quantile(x, 0.05, na.rm=T)
    upper[i] <- quantile(x, 0.95, na.rm=T)
  }
  plot(t, obsdata[[data_ix]]$data[,cell_ix,1], pch=ppd_symbol, ylab="CDK1 sensor level", xlab="Time (minutes)", ylim=ylim, type='n')
  have_plot <- which(!is.na(lower) & !is.na(upper))
  polygon(c(t[have_plot], rev(t[have_plot])), c(lower[have_plot], rev(upper[have_plot])), col=adjustcolor(ppd_color, alpha.f = 0.4), border=F)
  lines(t[have_plot], lower[have_plot], col=ppd_color, lwd=2)
  lines(t[have_plot], upper[have_plot], col=ppd_color, lwd=2)
  points(t, obsdata[[data_ix]]$data[,cell_ix,1], pch=ppd_symbol)
}

obsdata <- bcm3.cellpop.get.observed.data(misspec_fit, misspec_fit$likelihood$experiments[[experiment_ix]]$name)
merged <- list()
for (j in 1:length(obsdata)) {
  merged[[j]] <- list()
  for (l in 1:dim(obsdata[[j]]$data)[2]) {
    merged[[j]][[l]] <- matrix(NA, numppdsamples, length(obsdata[[j]]$time))
  }
}
for (i in 1:numppdsamples) {
  sample <- misspec_fit$posterior$samples[,temperature_ix,ppd_sample_ix[i]]
  simdata <- bcm3.cellpop.get.simulated.data(misspec_fit, misspec_fit$likelihood$experiments[[experiment_ix]]$name, sample)
  for (j in 1:length(obsdata)) {
    for (l in 1:dim(obsdata[[j]]$data)[2]) {
      merged[[j]][[l]][i,] <- simdata[[j]]$data[,l,1]
    }
  }
}
stdev_ix <- match(misspec_fit$likelihood$experiments[[experiment_ix]]$data[[data_ix]]$stdev, misspec_fit$variables)
sds <- 10^as.numeric(misspec_fit$posterior$samples[stdev_ix, temperature_ix, ppd_sample_ix])
t <- obsdata[[data_ix]]$time / 60
for (cell_ix in c(2,8,14)) {
  lower <- rep(NA, length(obsdata[[data_ix]]$time))
  upper <- rep(NA, length(obsdata[[data_ix]]$time))
  for (i in 1:length(obsdata[[data_ix]]$time)) {
    x <- rnorm(1000, merged[[data_ix]][[cell_ix]][,i], sds * misspec_fit$likelihood$experiments[[experiment_ix]]$data[[data_ix]]$stdev_multiplication_factor)
    lower[i] <- quantile(x, 0.05, na.rm=T)
    upper[i] <- quantile(x, 0.95, na.rm=T)
  }
  plot(t, obsdata[[data_ix]]$data[,cell_ix,1], pch=ppd_symbol, ylab="CDK1 sensor level", xlab="Time (minutes)", ylim=ylim, type='n')
  have_plot <- which(!is.na(lower) & !is.na(upper))
  polygon(c(t[have_plot], rev(t[have_plot])), c(lower[have_plot], rev(upper[have_plot])), col=adjustcolor(ppd_color, alpha.f = 0.4), border=F)
  lines(t[have_plot], lower[have_plot], col=ppd_color, lwd=2)
  lines(t[have_plot], upper[have_plot], col=ppd_color, lwd=2)
  points(t, obsdata[[data_ix]]$data[,cell_ix,1], pch=ppd_symbol)
}

merged <- list()
for (j in 1:length(obsdata)) {
  merged[[j]] <- list()
  for (l in 1:dim(obsdata[[j]]$data)[2]) {
    merged[[j]][[l]] <- matrix(NA, numppdsamples, length(obsdata[[j]]$time))
  }
}
for (i in 1:numppdsamples) {
  sample <- misspec_fit$posterior$samples[,calibrated_temperature_ix,ppd_sample_ix[i]]
  simdata <- bcm3.cellpop.get.simulated.data(misspec_fit, misspec_fit$likelihood$experiments[[experiment_ix]]$name, sample)
  for (j in 1:length(obsdata)) {
    for (l in 1:dim(obsdata[[j]]$data)[2]) {
      merged[[j]][[l]][i,] <- simdata[[j]]$data[,l,1]
    }
  }
}
stdev_ix <- match(misspec_fit$likelihood$experiments[[experiment_ix]]$data[[data_ix]]$stdev, misspec_fit$variables)
sds <- 10^as.numeric(misspec_fit$posterior$samples[stdev_ix, calibrated_temperature_ix, ppd_sample_ix])
t <- obsdata[[data_ix]]$time / 60
for (cell_ix in c(2,8,14)) {
  lower <- rep(NA, length(obsdata[[data_ix]]$time))
  upper <- rep(NA, length(obsdata[[data_ix]]$time))
  for (i in 1:length(obsdata[[data_ix]]$time)) {
    x <- rnorm(1000, merged[[data_ix]][[cell_ix]][,i], sds * misspec_fit$likelihood$experiments[[experiment_ix]]$data[[data_ix]]$stdev_multiplication_factor)
    lower[i] <- quantile(x, 0.05, na.rm=T)
    upper[i] <- quantile(x, 0.95, na.rm=T)
  }
  plot(t, obsdata[[data_ix]]$data[,cell_ix,1], pch=ppd_symbol, ylab="CDK1 sensor level", xlab="Time (minutes)", ylim=ylim, type='n')
  have_plot <- which(!is.na(lower) & !is.na(upper))
  polygon(c(t[have_plot], rev(t[have_plot])), c(lower[have_plot], rev(upper[have_plot])), col=adjustcolor(ppd_color, alpha.f = 0.4), border=F)
  lines(t[have_plot], lower[have_plot], col=ppd_color, lwd=2)
  lines(t[have_plot], upper[have_plot], col=ppd_color, lwd=2)
  points(t, obsdata[[data_ix]]$data[,cell_ix,1], pch=ppd_symbol)
}
res <- dev.off()
