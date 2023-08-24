library(pals)
library(ellipse)
source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/stats.r", sep=""))

model1 <- bcm3.load.results(".", "output_globalcov")
model2 <- bcm3.load.results(".", "output_autoblock")
model3 <- bcm3.load.results(".", "output_gmm")

cols <- brewer.set1(5)

sample_ix <- 6001:12000
temperature_ix <- 8
var_ix <- 2

lag_max <- 70
acfs <- list()
acfs[[1]] <- acf(t(model1$posterior$samples[var_ix,,sample_ix]), lag.max = lag_max, plot=F)
acfs[[2]] <- acf(t(model2$posterior$samples[var_ix,,sample_ix]), lag.max = lag_max, plot=F)
acfs[[3]] <- acf(t(model3$posterior$samples[var_ix,,sample_ix]), lag.max = lag_max, plot=F)

boundary <- qnorm((1 + 0.95)/2)/sqrt(length(sample_ix))

dir.create("figures", showWarnings = F)

pdf("figures/acf_plots.pdf", width=5, height=5)
plot(0, type='n', xlim=c(0,lag_max*1.1), ylim=c(-0.1,1), xlab="Lag", ylab="Correlation")
lines(acfs[[1]]$lag[,temperature_ix,temperature_ix], acfs[[1]]$acf[,temperature_ix,temperature_ix], col=cols[1], lwd=2)
lines(acfs[[2]]$lag[,temperature_ix,temperature_ix], acfs[[2]]$acf[,temperature_ix,temperature_ix], col=cols[2], lwd=2)
lines(acfs[[3]]$lag[,temperature_ix,temperature_ix], acfs[[3]]$acf[,temperature_ix,temperature_ix], col=cols[3], lwd=2)
legend("topright", legend=c("Global covariance", "Autoblock", "Gaussian mixture"), col=cols, lwd=2)
abline(h=0, lty=1)
abline(h=c(boundary,-boundary), lty=2)
res <- dev.off()

pdf("figures/sample_traces.pdf", width=10, height=8, useDingbats=F)
par(mfrow=c(3,1))
plot(model1$posterior$samples[var_ix,temperature_ix,], pch=20, cex=0.1, xlim=c(0, 12000), ylim=c(-1.3,-0.7), main="Global covariance", xlab="Sample iteration", ylab="Parameter value")
abline(v=c(1999.5,3999.5), lty=2, col='grey')
abline(v=5999.5, lty=2, col='grey', lwd=2)
plot(model2$posterior$samples[var_ix,temperature_ix,], pch=20, cex=0.1, xlim=c(0, 12000), ylim=c(-1.3,-0.7), main="Autoblock", xlab="Sample iteration", ylab="Parameter value")
abline(v=c(1999.5,3999.5), lty=2, col='grey')
abline(v=5999.5, lty=2, col='grey', lwd=2)
plot(model3$posterior$samples[var_ix,temperature_ix,], pch=20, cex=0.1, xlim=c(0, 12000), ylim=c(-1.3,-0.7), main="Gaussian mixture", xlab="Sample iteration", ylab="Parameter value")
abline(v=c(1999.5,3999.5), lty=2, col='grey')
abline(v=5999.5, lty=2, col='grey', lwd=2)
res <- dev.off()

clustcols <- brewer.set1(13)
iter <- length(model3$sampler_adaptation)
adapt_sample_ix <- 4001:6000
name <- names(model3$sampler_adaptation)[iter]
group <- model3$sampler_adaptation[[name]][[1]]
nclusters <- (length(group)-2)/3

pdf("figures/sample_2d_traces.pdf", width=8, height=4, useDingbats=F)
par(mfrow=c(1,2))
i <- 1
for (j in c(2,3)) {
  x <- model3$posterior$samples[i,temperature_ix,adapt_sample_ix]
  y <- model3$posterior$samples[j,temperature_ix,adapt_sample_ix]

  plot(x, y, pch='.', main=paste(group$variable_indices[i]+1, group$variable_indices[j]+1, sep="-"),
       xlab=model3$variables[i], ylab=model3$variables[j])
  for (clusti in 1:nclusters) {
    mean <- group[[paste("cluster", clusti-1, "_mean", sep="")]][c(i,j)]
    cov <- group[[paste("cluster", clusti-1, "_covariance", sep="")]][c(i,j),c(i,j)]
    ell <- ellipse(cov, centre = mean, level=0.6, draw=F)
    lines(ell, col=clustcols[clusti], lwd=2)
  }
}
res <- dev.off()

variable_summary(model1)
variable_summary(model2)
variable_summary(model3)
