library(beeswarm)
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))

model <- bcm3.load.results("../2-mitosis-model",
                           output_folder="output_gavetpines",
                           prior_file="prior_gavetpines.xml",
                           likelihood_file="../3-mitosis-model-test-prediction/likelihood_gavetpines_timings.xml")

timings <- list()
meraldi <- H5File$new("../2-mitosis-model/data/Meraldi2004.nc", 'r')
timings[["meraldi"]] <- meraldi[["meraldi_fig1b/NEBD_to_AO_duration"]][]
meraldi$close_all()
zhou <- H5File$new("../2-mitosis-model/data/Zhou2017_20cells.nc", 'r')
timings[["zhou"]] <- zhou[["Zhou_FigS1E/NEBD_to_AO_duration"]][]
zhou$close_all()
lu <- H5File$new("../2-mitosis-model/data/Lu2013_20cells.nc", 'r')
timings[["lu"]] <- lu[["Lu_Fig5E/NEBD_to_AO_duration"]][]
lu$close_all()
liu <- H5File$new("../2-mitosis-model/data/Liu.nc", 'r')
timings[["liu"]] <- liu[["Liu_fig4e/NEBD_to_AO_duration"]][]
liu$close_all()

logloss <- readRDS("logloss.rds")
pps <- readRDS("pps.rds")

summed <- t(logloss[[1]]$logloss) + t(logloss[[2]]$logloss) + t(logloss[[3]]$logloss) + t(logloss[[4]]$logloss)
Rlogloss_log <- apply(log10(t(summed)), 1, mean)
calibrated_temperature_ix <- which.min(Rlogloss_log)

calibrated_temperature_ix
model$posterior$temperatures[calibrated_temperature_ix]

dir.create("figures", showWarnings = F)

pdf("figures/log_loss_beeswarm.pdf", width=60/25.4, height=70/25.4)
par(cex=0.5, mgp=c(2.0,0.8,0))
plot_sample_ix <- sample(nrow(summed), 200)
beeswarm(data.frame(log10(summed[plot_sample_ix,])), pch=19, cex=0.5, spacing=0.1, xlab="Learning rate", ylab="Log loss", las=2, col=rgb(0,0,0,0.05), xaxt='n', yaxt='n', ylim=c(2.5,6.0))
axis(2, at=seq(3,6), labels=format(10^(3:6)))
interpolating_spline <- spline(model$posterior$temperatures, y=1:length(model$posterior$temperatures), xout=c(0.01, 0.1, 0.5))
axis(1, at=c(1,interpolating_spline$y, length(model$posterior$temperatures)), labels=c(0, interpolating_spline$x, 1))
lines(1:length(model$posterior$temperatures), Rlogloss_log, lwd=2, col="#4053d3")
points(calibrated_temperature_ix, Rlogloss_log[calibrated_temperature_ix], pch=8, col="#b51d14", lwd=2, cex=1.5)
res <- dev.off()

pdf("figures/predictions_from_gavetpines.pdf", width=120/25.4, height=250/25.4)
par(mfrow=c(4,1))
par(cex=0.5, mgp=c(2.0,0.8,0))
for (temperature_ix in c(calibrated_temperature_ix, dim(model$posterior$samples)[2])) {
  tmp <- matrix(unlist(pps[[temperature_ix]]$merged[[1]]), nrow=length(pps[[temperature_ix]]$ppd_sample_ix))
  
  lower <- apply(tmp, 2, quantile, probs=0.05, na.rm=T)
  mean <- apply(tmp, 2, mean, na.rm=T)
  upper <- apply(tmp, 2, quantile, probs=0.95, na.rm=T)
  
  for (maxy in c(150, 360)) {
    plot(0, type='n', xlim=c(0,200), ylim=c(0, maxy), xlab="", ylab="Time (minutes)", xaxt='n')
    ordering <- order(mean)
    for (i in 1:ncol(tmp)) {
      lines(c(i,i), c(lower[ordering[i]], upper[ordering[i]])/60)
    }
    abline(h=min(lower)/60, lty=2, col='grey')
    for (j in 1:length(timings)) {
      xinc <- 36 / length(timings[[j]])
      points(40 + (j-1)*40 + 1:length(timings[[j]]) * xinc, sort(timings[[j]])/60, pch=4)
    }
  }
}
res <- dev.off()


pdf("figures/densities.pdf", width=180/25.4, height=100/25.4)
par(mfrow=c(2,4))
par(cex=0.5, mgp=c(2.0,0.8,0))

plot_variable_distribution(model, var_ix=3, temperature_ix=dim(model$posterior$samples)[2])
plot_variable_distribution(model, var_ix=13, temperature_ix=dim(model$posterior$samples)[2])
plot_variable_distribution(model, var_ix=7, temperature_ix=dim(model$posterior$samples)[2])
plot_variable_distribution(model, var_ix=4, temperature_ix=dim(model$posterior$samples)[2])

plot_variable_distribution(model, var_ix=3, temperature_ix=calibrated_temperature_ix)
plot_variable_distribution(model, var_ix=13, temperature_ix=calibrated_temperature_ix)
plot_variable_distribution(model, var_ix=7, temperature_ix=calibrated_temperature_ix)
plot_variable_distribution(model, var_ix=4, temperature_ix=calibrated_temperature_ix)
res <- dev.off()
