source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))
source("../scripts/cellpop_functions.r")

dir.create("figures", showWarnings = F)

model <- bcm3.load.results(".", "output_gookin_t12n25e2_j12k16", prior_file="prior_gookin.xml", likelihood_file="likelihood_gookin.xml", load_sampler_adaptation = F)

numppdsamples <- 1000
numthreads <- 8
model <- bcm3.init.cpp(model, "", numthreads)

pp <- get_posterior_predictive(model, 1, 12, numppdsamples)

pdf("figures/posterior_predictive_gookin.pdf", width=160/25.4, height=80/25.4)
par(mfrow=c(2,4))
par(cex=0.5, mgp=c(2.0,0.8,0))
posterior_predictive_plot(model, pp, 1, data_ix=1, cell_ixs=28, ylim=c(-0.6,0.6), ppd_point_symbol=4, ppd_color = "#00b25d")

posterior_predictive_plot(model, pp, 1, data_ix=1, cell_ixs=48, ylim=c(-0.6,0.6), ppd_point_symbol=4, ppd_color = "#00b25d")
posterior_predictive_plot(model, pp, 1, data_ix=1, cell_ixs=49, ylim=c(-0.6,0.6), ppd_point_symbol=4, ppd_color = "#00b25d")
posterior_predictive_plot(model, pp, 1, data_ix=1, cell_ixs=50, ylim=c(-0.6,0.6), ppd_point_symbol=4, ppd_color = "#00b25d")

posterior_predictive_plot(model, pp, 1, data_ix=1, cell_ixs=33, ylim=c(-0.6,0.6), ppd_point_symbol=4, ppd_color = "#00b25d")

posterior_predictive_plot(model, pp, 1, data_ix=1, cell_ixs=76, ylim=c(-0.6,0.6), ppd_point_symbol=4, ppd_color = "#00b25d")
posterior_predictive_plot(model, pp, 1, data_ix=1, cell_ixs=77, ylim=c(-0.6,0.6), ppd_point_symbol=4, ppd_color = "#00b25d")
posterior_predictive_plot(model, pp, 1, data_ix=1, cell_ixs=78, ylim=c(-0.6,0.6), ppd_point_symbol=4, ppd_color = "#00b25d")
res <- dev.off()

load("prediction.rda")
epitopes <- names(observed)

pdf("figures/predictions_model_quantitative.pdf", width=150/25.4, height=100/25.4)
par(mfrow=c(2,3))
par(cex=0.5, mgp=c(2.0,0.8,0))
for (i in 1:6) {
  use_ix <- which(!is.na(observed[[i]]) & observed[[i]] != 0)
  x <- predicted[[i]][use_ix,1]
  y <- observed[[i]][use_ix]
  res <- cor.test(x, y, method="pearson")
  print(res)
  plot(x, y, main=paste(names(observed)[[i]], "; r=", format(res$estimate, digits=2), ", ", ifelse(res$p.value < 1e-16, "p < 1e-16", paste("p=", format(res$p.value, digits=2), sep="")), sep=""),
       xlab="Predicted", ylab="Observed (a.u.)", pch=20, col=rgb(0,0,0,0.2))
}
par(mfrow=c(1,1))
res <- dev.off()