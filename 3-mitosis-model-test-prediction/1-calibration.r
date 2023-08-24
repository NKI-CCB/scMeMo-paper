source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))
source("../scripts/cellpop_functions.r")

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

sample_ix <- (dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3])
numppdsamples <- dim(model$posterior$samples)[3]/2
experiment_ix <- 1

model <- bcm3.init.cpp(model, threads=16)

stdevfactors <- rep(1305,481.7302,1207.623,990.0,244.5993)
logloss <- list()
for (dsi in 1:length(timings)) {
  logloss[[dsi]] <- calculate_timing_log_loss(model, timings[[dsi]], stdevfactor=stdevfactors[dsi], numppdsamples=numppdsamples)
}
saveRDS(logloss, "logloss.rds")

summed <- t(logloss[[1]]$logloss) + t(logloss[[2]]$logloss) + t(logloss[[3]]$logloss) + t(logloss[[4]]$logloss)
Rlogloss_log <- apply(log10(t(summed)), 1, mean)
calibrated_temperature_ix <- which.min(Rlogloss_log)

pps <- list()
for (temperature_ix in c(calibrated_temperature_ix, dim(model$posterior$samples)[2])) {
  cat("Temperature", temperature_ix, "\n")
  pps[[temperature_ix]] <- get_posterior_predictive(model, experiment_ix, temperature_ix, numppdsamples = numppdsamples)
}
saveRDS(pps, "pps.rds")
