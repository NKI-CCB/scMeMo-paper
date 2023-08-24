library(XML)
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))

input_model_filename <- "../5-cell-cycle-model/model.xml"
output_model_filename <- "model_setparams.xml"

modelxml <- xmlTreeParse(input_model_filename)
root <- xmlRoot(modelxml)
npar <- length(root[["model"]][["listOfParameters"]])

bcmres <- list()
bcmres[["akopyan"]]  <- bcm3.load.results("../5-cell-cycle-model", "output_akopyan_t24n50e2_j12k16", prior_file="prior_akopyan.xml", likelihood_file="likelihood_akopyan.xml", load_sampler_adaptation = F)
bcmres[["barr2016"]] <- bcm3.load.results("../5-cell-cycle-model", "output_barr2016_t24n50e2_j12k16", prior_file="prior_barr2016.xml", likelihood_file="likelihood_barr2016.xml", load_sampler_adaptation = F)
bcmres[["barr2017"]] <- bcm3.load.results("../5-cell-cycle-model", "output_barr2017_t24n20e2_j12k16", prior_file="prior_barr2017.xml", likelihood_file="likelihood_barr2017.xml", load_sampler_adaptation = F)
bcmres[["cappell"]]  <- bcm3.load.results("../5-cell-cycle-model", "output_cappell_t24n20e2_j12k16", prior_file="prior_cappell.xml", likelihood_file="likelihood_cappell.xml", load_sampler_adaptation = F)

set_param_values <- c()

temperature_ix <- 14

# Values from Barr2016
max_posterior_sample <- which.max(bcmres[["barr2016"]]$posterior$lfracposterior[temperature_ix,])
set_param_values["p27_transcription"]          <- 10^bcmres[["barr2016"]]$posterior$samples[which(bcmres[["barr2016"]]$variables == "p27_transcription"), temperature_ix, max_posterior_sample]
set_param_values["CyclinE_transcription"]      <- 10^bcmres[["barr2016"]]$posterior$samples[which(bcmres[["barr2016"]]$variables == "CyclinE_transcription"), temperature_ix, max_posterior_sample]
set_param_values["CyclinE_degradation"]        <- 10^bcmres[["barr2016"]]$posterior$samples[which(bcmres[["barr2016"]]$variables == "CyclinE_degradation"), temperature_ix, max_posterior_sample]
set_param_values["active_CyclinE_degradation"] <- 10^bcmres[["barr2016"]]$posterior$samples[which(bcmres[["barr2016"]]$variables == "active_CyclinE_degradation"), temperature_ix, max_posterior_sample]
set_param_values["CyclinA_transcription"]      <- 10^bcmres[["barr2016"]]$posterior$samples[which(bcmres[["barr2016"]]$variables == "CyclinA_transcription"), temperature_ix, max_posterior_sample]

# Values from Barr2017
max_posterior_sample <- which.max(bcmres[["barr2017"]]$posterior$lfracposterior[temperature_ix,])
set_param_values["p21_transcription"]          <- 10^bcmres[["barr2017"]]$posterior$samples[which(bcmres[["barr2017"]]$variables == "p21_transcription"), temperature_ix, max_posterior_sample]
set_param_values["p21_degradation"]            <- 10^bcmres[["barr2017"]]$posterior$samples[which(bcmres[["barr2017"]]$variables == "p21_degradation"), temperature_ix, max_posterior_sample]

# Values from Akopyan
max_posterior_sample <- which.max(bcmres[["akopyan"]]$posterior$lfracposterior[temperature_ix,])
set_param_values["CyclinB_transcription"]      <- 10^bcmres[["akopyan"]]$posterior$samples[which(bcmres[["akopyan"]]$variables == "CyclinB_transcription"), temperature_ix, max_posterior_sample]

# Values from Cappell
max_posterior_sample <- which.max(bcmres[["cappell"]]$posterior$lfracposterior[temperature_ix,])
set_param_values["Geminin_transcription"]      <- 10^bcmres[["cappell"]]$posterior$samples[which(bcmres[["cappell"]]$variables == "Geminin_transcription"), temperature_ix, max_posterior_sample]
set_param_values["Emi1_transcription"]         <- 10^bcmres[["cappell"]]$posterior$samples[which(bcmres[["cappell"]]$variables == "Emi1_transcription"), temperature_ix, max_posterior_sample]

for (i in 1:npar) {
  varname <- xmlAttrs(root[["model"]][["listOfParameters"]][[i]])[["id"]]
  
  variable_ix <- match(varname, names(set_param_values))
  if (!is.na(variable_ix)) {
    cat("Set parameter", varname, "to", set_param_values[variable_ix], "\n", sep=" ")
    xmlAttrs(root[["model"]][["listOfParameters"]][[i]])[["value"]] <- set_param_values[variable_ix]
  }
}

saveXML(root, output_model_filename, prefix = '<?xml version="1.0" encoding="UTF-8"?>\n')
