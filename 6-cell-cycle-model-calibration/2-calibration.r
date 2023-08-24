library(beeswarm)
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))
source("../scripts/cellpop_functions.r")

dir.create("figures", showWarnings = F)

models <- list()
models[["akopyan"]]    <- bcm3.load.results("../5-cell-cycle-model", "output_akopyan_t24n50e2_j12k16", prior_file="prior_akopyan.xml", likelihood_file="../6-cell-cycle-model-calibration/likelihood_akopyan_sphasepred.xml", load_sampler_adaptation = F)
models[["cappell"]]    <- bcm3.load.results("../5-cell-cycle-model", "output_cappell_t24n20e2_j12k16", prior_file="prior_cappell.xml", likelihood_file="../6-cell-cycle-model-calibration/likelihood_cappell_sphasepred.xml", load_sampler_adaptation = F)
models[["westendorp"]] <- bcm3.load.results("../5-cell-cycle-model", "output_westendorp_t24n20e2_j12k16", prior_file="prior_westendorp.xml", likelihood_file="../6-cell-cycle-model-calibration/likelihood_westendorp_sphasepred.xml", load_sampler_adaptation = F)

timings <- list()
burgess <- H5File$new("data/Burgess2012.nc", 'r')
timings[["burgess"]] <- burgess[["Burgess_Fig2B/Sphase_duration"]][]
burgess$close_all()
grant <- H5File$new("data/Grant2018_15cells.nc", 'r')
timings[["grant"]] <- grant[["Grant2018_Fig4C/Sphase_duration_U2OS"]][]
grant$close_all()

numppdsamples <- 1000
numthreads <- 8

loglosses <- list()
loglosses[["akopyan"]] <- readRDS("log_loss_akopyan.rds")
loglosses[["cappell"]] <- readRDS("log_loss_cappell.rds")
loglosses[["westendorp"]] <- readRDS("log_loss_westendorp.rds")

pdf("figures/log_loss.pdf", width=200/25.4, height=70/25.4)
par(cex=0.5, mgp=c(2.0,0.8,0), mar=c(4,3,3,1))
par(mfcol=c(1,3))
for (i in 1:3) {
  Rlogloss <- apply(log10(loglosses[[i]]$logloss), 1, mean)
  calibrated_temperature_ix <- which.min(Rlogloss)
  
  subsample <- sample(ncol(loglosses[[i]]$logloss), 200)
  
  beeswarm(data.frame(log10(t(loglosses[[i]]$logloss[,subsample]))), pch=19, cex=0.5, spacing=0.1, xlab="Learning rate", ylab="Log loss", las=2, col=rgb(0,0,0,0.05), xaxt='n', yaxt='n', ylim=c(1.5,6.5), main=paste(names(models)[i]))
  axis(2, at=seq(2,6), labels=format(10^(2:6)))
  interpolating_spline <- spline(models[[i]]$posterior$temperatures, y=1:length(models[[i]]$posterior$temperatures), xout=c(0.01, 0.1, 0.5))
  axis(1, at=c(1,interpolating_spline$y, 24), labels=c(0, interpolating_spline$x, 1))
  lines(1:24, apply(log10(loglosses[[i]]$logloss), 1, mean), lwd=2, col="#4053d3")
  points(calibrated_temperature_ix, Rlogloss[calibrated_temperature_ix], pch=8, col="#b51d14", lwd=2, cex=1.5)
  
  cat(names(loglosses)[i], "min ix =", calibrated_temperature_ix, "; temperature =", models[[i]]$posterior$temperatures[calibrated_temperature_ix], "\n")
}
par(mfrow=c(1,1))
res <- dev.off()


which_timings <- list()
which_timings[["akopyan"]] <- "grant"
which_timings[["cappell"]] <- "burgess"
which_timings[["westendorp"]] <- "burgess"

for (i in 1:3) {
  cat(names(models)[i], "\n")
  Rlogloss <- apply(log10(loglosses[[i]]$logloss), 1, mean)
  calibrated_temperature_ix <- which.max(Rlogloss)
  full_temperature_ix <- 24

  models[[i]] <- bcm3.init.cpp(models[[i]], "", numthreads)
  
  pps <- list()
  pps[[calibrated_temperature_ix]] <- get_posterior_predictive(models[[i]], 1, calibrated_temperature_ix, numppdsamples = numppdsamples)
  pps[[full_temperature_ix]] <- get_posterior_predictive(models[[i]], 1, full_temperature_ix, numppdsamples = numppdsamples)
  
  pdf(paste("figures/Sphase_timing_predictions_", names(models)[i], ".pdf", sep=""), width=100/25.4, height=50/25.4)
  par(mfrow=c(1,2))
  par(cex=0.5, mgp=c(2.0,0.8,0))
  
  for (temperature_ix in c(full_temperature_ix, calibrated_temperature_ix)) {
    lower <- sapply(pps[[temperature_ix]]$merged[[1]], quantile, probs=0.05, na.rm=T)
    upper <- sapply(pps[[temperature_ix]]$merged[[1]], quantile, probs=0.95, na.rm=T)
    
    nsamples <- 15
    
    plot(0, type='n', xlim=c(0,33), ylim=c(0, 20), xlab="", ylab="Time (hours)", xaxt='n')
    
    ordering <- order(lower)
    for (j in 1:length(lower)) {
      lines(c(j,j), c(lower[ordering[j]]/3600, upper[ordering[j]]/3600))
    }
    timname <- which_timings[[names(models)[i]]]
    points(nsamples+2 + 1:length(timings[[timname]]), sort(timings[[timname]])/3600, pch=4)
    
    if (temperature_ix == full_temperature_ix) {
      text(x=8, y=-2, labels="Prediction", xpd=T)
    } else {
      text(x=8, y=-3, labels="Calibrated\nprediction", xpd=T)
    }
    text(x=25, y=-2, labels=timname, xpd=T)
  }
  
  dev.off()
}
