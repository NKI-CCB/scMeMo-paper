library(RcppHungarian)
source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))
source("../scripts/cellpop_functions.r")

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
numthreads <- 16

for (i in 1:length(models)) {
  models[[i]] <- bcm3.init.cpp(models[[i]], "", numthreads)
}

ll <- calculate_timing_log_loss(models[["akopyan"]], timings[["grant"]], stdevfactor=7802.144, numppdsamples=numppdsamples)
saveRDS(ll, "log_loss_akopyan.rds")

ll <- calculate_timing_log_loss(models[["cappell"]], timings[["burgess"]], stdevfactor=3773.131, numppdsamples=numppdsamples)
saveRDS(ll, "log_loss_cappell.rds")

ll <- calculate_timing_log_loss(models[["westendorp"]], timings[["burgess"]], stdevfactor=3773.131, numppdsamples=numppdsamples)
saveRDS(ll, "log_loss_westendorp.rds")

for (i in 1:length(models)) {
  models[[i]] <- bcm3.release.cpp(models[[i]])
}
