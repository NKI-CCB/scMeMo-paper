source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))
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

numthreads <- 16
numppdsamples <- 1500

pps <- list()
for (i in 1:length(models)) {
  cat(names(models)[i], "\n")
  models[[i]] <- bcm3.init.cpp(models[[i]], "", numthreads)
  pps[[i]] <- get_posterior_predictive(models[[i]], 1, 24, numppdsamples)
  models[[i]] <- bcm3.release.cpp(models[[i]])
  cat("\n")
}
names(pps) <- names(models)
saveRDS(pps, file="pps.rds")
