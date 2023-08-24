source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))

model <- bcm3.load.results(".", "output_gookin_t12n25e2_j12k16", prior_file="prior_gookin.xml", likelihood_file="likelihood_gookin_predict.xml")

epitopes <- c("CycB1", "Geminin", "p21", "CycE", "Cdt1", "CycA2")

predictions <- list()

data <- H5File$new("data/gookin2017_all_cells.nc", 'r')
sensor_data <- data[["Gookin_all_cells"]][["CDK2_sensor_log"]][1,,]
if_data <- list()
parent_data <- data[["Gookin_all_cells"]][["parent"]][]
for (e in epitopes) {
  if_data[[e]] <- data[["Gookin_all_cells"]][[paste("IF_", e, sep="")]][1,,1]
}
data$close_all()

temperature_ix <- 12
experiment_ix <- 1
nthreads <- 16

model <- bcm3.init.cpp(model, "--cellpop.use_only_cell_ix=0", nthreads)

nspecies <- bcm3.cellpop.get.num.species(model, model$likelihood$experiments[[1]]$name)
species_names <- rep(NA, nspecies)
for (i in 1:nspecies) {
  species_names[i] <- bcm3.cellpop.get.species.name(model, model$likelihood$experiments[[1]]$name, i)
}

model_species <- list()
model_species[["CycB1"]]   <- c("CyclinB", "active_CyclinB_CDK1")
model_species[["Geminin"]] <- c("Geminin")
model_species[["p21"]]     <- c("p21")
model_species[["CycE"]]    <- c("CyclinE", "active_CyclinE_CDK2")
model_species[["Cdt1"]]    <- c("Cdt1")
model_species[["CycA2"]]   <- c("CyclinA", "active_CyclinA_CDK1", "active_CyclinA_CDK2")

observed <- list()
predicted <- list()

for (e in epitopes) {
  observed[[e]] <- rep(NA, length(parent_data))
  predicted[[e]] <- matrix(NA, length(parent_data), 1)
  
  species_ix <- which(species_names %in% model_species[[e]])

  sample_ix <- which.max(model$posterior$lprior[temperature_ix,] + model$posterior$llikelihood[temperature_ix,])
  sample <- model$posterior$samples[,temperature_ix,dim(model$posterior$samples)[3]]
  names(sample) <- model$variables
  
  which_data_cell_ix <- which(!is.na(if_data[[e]]))
  which_data_cell_ix <- sort(unique(c(which_data_cell_ix, parent_data[which_data_cell_ix])))
  
  cat("Predicting", e, "\n")
  pb <- txtProgressBar(max=length(which_data_cell_ix), style=3)
  j <- which_data_cell_ix[1]
  count <- 0
  while (T) {
    if (j > tail(which_data_cell_ix, n=1)) {
      break
    }
    
    cell_ix <- j
    j <- j + 1
    
    cells <- cell_ix
    for (k in 1:2) {
      if (j < length(if_data[[e]])) {
        if (!is.na(parent_data[j]) && parent_data[j] == cell_ix) {
          cells <- c(cells, cell_ix+1)
          j <- j + 1
          count <- count + 1
        }
      }
    }
    
    model <- bcm3.reinit.cpp(model, paste("--cellpop.use_only_cell_ix=", paste(cells-1, collapse=","), sep=""), nthreads)

    obsdata <- bcm3.cellpop.get.observed.data(model, model$likelihood$experiments[[experiment_ix]]$name)
    traj <- bcm3.cellpop.get.matched.simulation(model, model$likelihood$experiments[[experiment_ix]]$name, sample)
    
    prediction <- rep(NA, dim(traj[[1]]$cells)[3])
    for (k in 1:length(prediction)) {
      prediction[k] <- sum(traj[[1]]$cells[species_ix,dim(traj[[1]]$cells)[2],k])
    }
  
    observed[[e]][cells] <- if_data[[e]][cells]
    predicted[[e]][cells,1] <- prediction
    
    setTxtProgressBar(pb, count)
  }
  cat("\n")
  
  epitope_result <- list()
  epitope_result$observed <- observed[[e]]
  epitope_result$predicted <- predicted[[e]]
}

save(observed, predicted, file="prediction.rda")

model <- bcm3.release.cpp(model)
