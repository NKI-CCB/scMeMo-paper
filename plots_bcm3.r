library(pals)
library(ellipse)
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))

#model_name <- "2-mitosis-model"
model_name <- "4-mitosis-model-posterior-calibration"
#model_name <- "5-cell-cycle-model"
#model_name <- "7-test-gookin-prediction"

# output_folder <- "output_akopyan"
# prior_file <- "prior_akopyan.xml"
# likelihood_file <- "likelihood_akopyan.xml"

# output_folder <- "output_gavetpines"
# prior_file <- "prior_gavetpines.xml"
# likelihood_file <- "likelihood_gavetpines.xml"

# output_folder <- "output_simdata_correctspecified_t24_n100"
# prior_file <- "prior_simulated_data.xml"
# likelihood_file <- "likelihood_simulated_data.xml"

# output_folder <- "output_simdata_misspecified_t24_n100"
# prior_file <- "prior_simulated_data_misspecified.xml"
# likelihood_file <- "likelihood_simulated_data_misspecified.xml"

output_folder <- "output_simdata_misspecified_with_mitosis_timings_t24_n100"
prior_file <- "prior_simulated_data_misspecified.xml"
likelihood_file <- "likelihood_simulated_data_misspecified_with_mitosis_timings.xml"

# output_folder <- "output_akopyan_t24n50e2_j12k16"
# prior_file <- "prior_akopyan.xml"
# likelihood_file <- "likelihood_akopyan.xml"

# output_folder <- "output_barr2016_t24n50e2_j12k16"
# prior_file <- "prior_barr2016.xml"
# likelihood_file <- "likelihood_barr2016.xml"

# output_folder <- "output_barr2017_t24n20e2_j12k16"
# prior_file <- "prior_barr2017.xml"
# likelihood_file <- "likelihood_barr2017.xml"

# output_folder <- "output_cappell_t24n20e2_j12k16"
# prior_file <- "prior_cappell.xml"
# likelihood_file <- "likelihood_cappell.xml"

# output_folder <- "output_eward_t24n20e2_j12k16"
# prior_file <- "prior_eward.xml"
# likelihood_file <- "likelihood_eward.xml"

# output_folder <- "output_eward_noCycArepression_t24n20e2_j12k16"
# prior_file <- "prior_eward.xml"
# likelihood_file <- "likelihood_eward_noCycArepression.xml"

# output_folder <- "output_eward_noE2F7_t24n20e2_j12k16"
# prior_file <- "prior_eward.xml"
# likelihood_file <- "likelihood_eward_noE2F7.xml"

# output_folder <- "output_westendorp_t24n20e2_j12k16"
# prior_file <- "prior_westendorp.xml"
# likelihood_file <- "likelihood_westendorp.xml"

# output_folder <- "output_gookin_t12n25e2_j12k16"
# prior_file <- "prior_gookin.xml"
# likelihood_file <- "likelihood_gookin.xml"

output_filename <- "output.nc"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) { model_name <- args[1] }
if (length(args) >= 2) { output_folder <- args[2] }
if (length(args) >= 3) { prior_file <- args[3] }
if (length(args) >= 4) { likelihood_file <- args[4] }
if (length(args) >= 5) { output_filename <- args[5] }

model <- bcm3.load.results(model_name, output_folder, prior_file=prior_file, likelihood_file=likelihood_file, output_filename=output_filename)

sample_ix <- (dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3])
temperature_ix <- dim(model$posterior$samples)[2]
numppdsamples <- 100
ppd_sample_ix <- sample(sample_ix, numppdsamples)
nthreads <- 4

print("Plotting densities")
png_tile(paste(model$output_folder, "/densities.png", sep=""), width=1600, height=1000, model$nvar)
par(mar=c(3,2,2,1))
for (i in 1:model$nvar) {
  plot_variable_distribution(model, var_ix=i, sample_ix=sample_ix)
}
par(mfrow=c(1,1))
res <- dev.off()

print("Plotting traces...")
png_tile(paste(model$output_folder, "/traces.png", sep=""), width=2000, height=1400, model$nvar)
par(mar=c(2,2,2,2))
for (i in 1:model$nvar) {
  plot_trace(model, var_ix=i)
}
par(mfrow=c(1,1))
res <- dev.off();

opt_ix <- NA

if (!is.na(opt_ix)) {
  print("Plotting densities")
  png_tile(paste(model$output_folder, "/densities_opt.png", sep=""), width=1600, height=1000, model$nvar)
  par(mar=c(3,2,2,1))
  for (i in 1:model$nvar) {
    plot_variable_distribution(model, var_ix=i, temperature_ix=opt_ix, sample_ix=sample_ix)
  }
  par(mfrow=c(1,1))
  res <- dev.off()
  
  print("Plotting traces...")
  png_tile(paste(model$output_folder, "/traces_opt.png", sep=""), width=2000, height=1400, model$nvar)
  par(mar=c(2,2,2,2))
  for (i in 1:model$nvar) {
    plot_trace(model, var_ix=i, temperature_ix=opt_ix)
  }
  par(mfrow=c(1,1))
  res <- dev.off();
}

if (model$nvar <= 10) {
  print("Plotting 2D scatter plots...")
  #for (k in 1:dim(model$posterior$samples)[2]) {
  k <- dim(model$posterior$samples)[2]
  png(paste(model$output_folder, "/2d_traces_", k, ".png", sep=""), width=3000, height=3000);
  par(mfrow=c(model$nvar, model$nvar))
  par(mar=c(1,1,1,1))
  for (i in 1:model$nvar) {
    for (j in 1:model$nvar) {
      if (i == j) {
        plot_variable_distribution(model, var_ix=i, temperature_ix=k)
      } else {
        xb <- prior_bounds(model, j)
        yb <- prior_bounds(model, i)
        plot(t(model$posterior$samples[c(j,i), k,sample_ix]), xlim=xb, ylim=yb, col=rgb(0,0,0,0.2), pch=20)
      }
    }
  }
  par(mfrow=c(1,1))
  res <- dev.off();
  #}
}

model <- bcm3.init.cpp(model, "", nthreads)

for (experiment_ix in 1:length(model$likelihood$experiments)) {
  
  #for (temperature_ix in seq(2, 32, by=2)) {
    
  sample <- model$posterior$samples[,temperature_ix,dim(model$posterior$samples)[3]]
  #sample <- model$posterior$samples[,temperature_ix,2569]
  names(sample) <- model$variables
  #logl <- bcm3.cellpop.get.likelihood(model, "", sample)
  #logl
  
  obsdata <- bcm3.cellpop.get.observed.data(model, model$likelihood$experiments[[experiment_ix]]$name)
  simdata <- bcm3.cellpop.get.simulated.data(model, model$likelihood$experiments[[experiment_ix]]$name, sample)
  
  merged <- list()
  for (j in 1:length(obsdata)) {
    merged[[j]] <- list()
    for (l in 1:dim(obsdata[[j]]$data)[2]) {
      merged[[j]][[l]] <- matrix(NA, numppdsamples, length(obsdata[[j]]$time))
    }
  }
  for (i in 1:numppdsamples) {
    sample <- model$posterior$samples[,temperature_ix,ppd_sample_ix[i]]
    names(sample) <- model$variables
    #sample["DNA_licensing"] <- sample["DNA_licensing"] + 0.1
    #logl1 <- bcm3.cellpop.get.likelihood(model, "", sample)
    simdata <- bcm3.cellpop.get.simulated.data(model, model$likelihood$experiments[[experiment_ix]]$name, sample)

    for (j in 1:length(obsdata)) {
      for (l in 1:dim(obsdata[[j]]$data)[2]) {
        vals <- simdata[[j]]$data[,l,1]
        vals[vals == -Inf] <- -2
        merged[[j]][[l]][i,] <- vals
      }
    }
  }
  
  plot_count <- 0
  for (k in 1:length(obsdata)) {
    if (model$likelihood$experiments[[experiment_ix]]$data[[k]]$type == "duration") {
      plot_count <- plot_count + 1
    } else {
      plot_count <- plot_count + dim(obsdata[[k]]$data)[2]
    }
  }
  
  png_tile(paste(model$output_folder, "/postpred_t", temperature_ix, "_", gsub("/", "_", model$likelihood$experiments[[experiment_ix]]$name), ".png", sep=""), 2000, 1200, plot_count)

  for (k in 1:length(obsdata)) {
    stdev_ix <- match(model$likelihood$experiments[[experiment_ix]]$data[[k]]$stdev, model$variables)
    if (is.na(stdev_ix)) {
      vars <- as.numeric(model$likelihood$experiments[[experiment_ix]]$data[[k]]$stdev)
    } else {
      vars <- 10^as.numeric(model$posterior$samples[stdev_ix, temperature_ix, ppd_sample_ix])
    }
    if (model$likelihood$experiments[[experiment_ix]]$data[[k]]$type == "duration") {
      tmp <- matrix(unlist(merged[[k]]), nrow=numppdsamples)
    
      ppd_barplot(model, tmp, obsdata[[k]]$data[1,,1], labels=1:ncol(tmp),
                  error_model = "normal", sd_samples = vars * model$likelihood$experiments[[experiment_ix]]$data[[k]]$stdev_multiplication_factor,
                  xlab="Cells", ylab="NEBD to anaphase time (seconds)", ylim=c(0,6000), pointsize=1.0)
      
    } else {
      t <- obsdata[[k]]$time
      for (cell_ix in 1:dim(obsdata[[k]]$data)[2]) {
        lower <- rep(NA, length(obsdata[[k]]$time))
        upper <- rep(NA, length(obsdata[[k]]$time))
        for (i in 1:length(obsdata[[k]]$time)) {
          x <- rnorm(1000, merged[[k]][[cell_ix]][,i], vars * model$likelihood$experiments[[experiment_ix]]$data[[k]]$stdev_multiplication_factor)
          lower[i] <- quantile(x, 0.05, na.rm=T)
          upper[i] <- quantile(x, 0.95, na.rm=T)
        }
        plot(obsdata[[k]]$time, obsdata[[k]]$data[,cell_ix,1], pch=19, ylab=model$likelihood$experiments[[experiment_ix]]$data[[k]]$data_name, xlab="Time (seconds)", 
             #ylim=c(-0.5, 1.2))
             ylim=c(min(obsdata[[k]]$data, lower, na.rm=T),
                    max(obsdata[[k]]$data, upper, na.rm=T)+0.2))
        ppd_color <- rgb(43,131,186, max=255)
        have_plot <- which(!is.na(lower) & !is.na(upper))
        polygon(c(t[have_plot], rev(t[have_plot])), c(lower[have_plot], rev(upper[have_plot])), col=adjustcolor(ppd_color, alpha.f = 0.4), border=F)
        lines(t[have_plot], lower[have_plot], col=ppd_color, lwd=2)
        lines(t[have_plot], upper[have_plot], col=ppd_color, lwd=2)
        for (i in 1:numppdsamples) {
          lines(t, merged[[k]][[cell_ix]][i,], col=rgb(0,0,0,0.1))
        }
        points(obsdata[[k]]$time, obsdata[[k]]$data[,cell_ix,1], pch=19)
      }
    }
  }
  
  res <- dev.off()
  #}
}

#model <- bcm3.release.cpp(model)
#model <- bcm3.init.cpp(model, nthreads)

for (experiment_ix in 1:length(model$likelihood$experiments)) {
  plot_sample_ix <- which.max(model$posterior$lposterior[temperature_ix,])
  sample <- model$posterior$samples[,temperature_ix,plot_sample_ix]
  names(sample) <- model$variables
  #sample["CyclinA_derepression_by_CDK2"] <- -4.5
  
  #bcm3.cellpop.get.likelihood(model, "", sample)
  traj <- bcm3.cellpop.get.simulated.trajectories(model, model$likelihood$experiments[[experiment_ix]]$name, sample)
  
  nspecies <- dim(traj$cells)[1]
  
  png_tile(paste(model$output_folder, "/trajectories_", gsub("/", "_", model$likelihood$experiments[[experiment_ix]]$name), ".png", sep=""), 2000, 1200, nspecies)
  par(mgp=c(1.6,0.8,0))
  par(mar=c(3,2,2,1))
  for (i in 1:nspecies) {
    plot(traj$time, traj$cells[i,,1], type='n', main=bcm3.cellpop.get.species.name(model, model$likelihood$experiments[[experiment_ix]]$name, i), ylim=c(0, 2), xlab='Time (seconds)', ylab='', xlim=range(traj$time))
    for (j in 1:dim(traj$cells)[3]) {
      lines(traj$time, traj$cells[i,,j])
    }
  }
  res <- dev.off()
  
  # i <- 9
  # plot(traj$time, traj$cells[i,,1], type='n', main=bcm3.cellpop.get.species.name(model, model$likelihood$experiments[[experiment_ix]]$name, i), ylim=c(0, 2), xlab='Time (seconds)', ylab='', xlim=range(traj$time))
  # for (j in 1:dim(traj$cells)[3]) {
  #   lines(traj$time, traj$cells[9,,j]+traj$cells[54,,j]+traj$cells[57,,j])
  # }
  
  # sample["mitogenic_signal_rate]"] <- -3
  # sample["mitogenic_signal_rate"] <- sample["mitogenic_signal_rate"] + 0.001
  # bcm3.cellpop.get.likelihood(model, "", sample)
  # traj <- bcm3.cellpop.get.simulated.trajectories(model, model$likelihood$experiments[[experiment_ix]]$name, sample)
  # #
  # plot(traj$time, traj$cells[23,,1]+traj$cells[24,,1], type='l', main=bcm3.cellpop.get.species.name(model, model$likelihood$experiments[[experiment_ix]]$name, i), ylim=c(min(0, sapply(c(traj$cells), max, 0, na.rm=T), na.rm=T),1.2), xlab='Time (seconds)', ylab='Concentration')
  # for (j in 1:dim(traj$cells)[3]) {
  #   lines(traj$time, traj$cells[23,,j]+traj$cells[24,,j])
  # }
  # plot(traj$time, traj$cells[25,,1]+traj$cells[26,,1], type='n', main=bcm3.cellpop.get.species.name(model, model$likelihood$experiments[[experiment_ix]]$name, i), ylim=c(min(0, sapply(c(traj$cells), max, 0, na.rm=T), na.rm=T),1.2), xlab='Time (seconds)', ylab='Concentration')
  # for (j in 1:64) {
  #   lines(traj$time, traj$cells[25,,j]+traj$cells[26,,j])S
  # }
  # 
  # plot(traj$time, traj$cells[23,,1]+traj$cells[24,,1], xlim=c(0, 50000), type='n', main=bcm3.cellpop.get.species.name(model, model$likelihood$experiments[[experiment_ix]]$name, i), ylim=c(min(0, sapply(c(traj$cells), max, 0, na.rm=T), na.rm=T),1.2), xlab='Time (seconds)', ylab='Concentration')
  # for (j in 1:dim(traj$cells)[3]) {
  #   ana_onset <- head(which(traj$cells[11,,j] > 1e-3), n=1)
  #   if (length(ana_onset) > 0) {
  #     offset <- traj$time[ana_onset]
  #     lines(traj$time-offset, traj$cells[23,,j]+traj$cells[24,,j])
  #   }
  # }
  # Heatmap(t(traj$cells[29,,]), cluster_rows = F, cluster_columns = F, col=cividis())
  # Heatmap(t(traj$cells[30,,]+traj$cells[37,,]), cluster_rows = F, cluster_columns = F, col=cividis())
  
  #Heatmap(t(log10(traj$cells[29,,]/traj$cells[28,,])), cluster_rows=F, cluster_columns=F, col=cividis())
}
model <- bcm3.release.cpp(model)


if (!is.null(model$sampler_adaptation)) {
  clustcols <- brewer.set1(13)
  
  for (iter in 1:length(model$sampler_adaptation)) {
    #adapt_sample_ix <- seq((iter-1)*2000+1, iter*2000, by=3)
    adapt_sample_ix <- seq((iter-1)*2000+1, iter*2000)
    name <- names(model$sampler_adaptation)[iter]
    for (block in 1:length(model$sampler_adaptation[[iter]])) {
      blockname <- names(model$sampler_adaptation[[iter]])[block]
      group <- model$sampler_adaptation[[name]][[blockname]]
      
      png(paste(model$output_folder, "/sampler_adaptation_", name, "_", blockname, ".png", sep=""), width=model$nvar*250, height=model$nvar*250)
      par(mfcol=c(model$nvar, model$nvar))
      par(mar=c(4,4,4,1))
      for (i in 1:model$nvar) {
        for (j in 1:model$nvar) {
          x <- model$posterior$samples[i,temperature_ix,adapt_sample_ix]
          y <- model$posterior$samples[j,temperature_ix,adapt_sample_ix]
          nclusters <- (length(group)-2)/3
          
          plot(x, y, pch='.', main=paste(group$variable_indices[i]+1, group$variable_indices[j]+1, sep="-"))
          for (clusti in 1:nclusters) {
            mean <- group[[paste("cluster", clusti-1, "_mean", sep="")]][c(i,j)]
            cov <- group[[paste("cluster", clusti-1, "_covariance", sep="")]][c(i,j),c(i,j)]
            ell <- ellipse(cov, centre = mean, level=0.6, draw=F)
            lines(ell, col=clustcols[clusti], lwd=2)
          }
        }
      }
      res <- dev.off()
    }
  }
}

marginal_likelihood(model)


