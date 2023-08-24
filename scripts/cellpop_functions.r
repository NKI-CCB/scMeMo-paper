library(RcppHungarian)

get_posterior_predictive <- function(bcm3, experiment_ix, temperature_ix, numppdsamples=dim(model$posterior$samples)[3]/2)
{
  results <- list()
  
  model <- bcm3
  results$obsdata <- bcm3.cellpop.get.observed.data(model, model$likelihood$experiments[[experiment_ix]]$name)
  
  sample_ix <- (dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3])
  if (numppdsamples == dim(model$posterior$samples)[3]/2) {
    results$ppd_sample_ix <- sample_ix
  } else {
    results$ppd_sample_ix <- sample(sample_ix, numppdsamples)
  }
  
  results$merged <- list()
  for (j in 1:length(results$obsdata)) {
    results$merged[[j]] <- list()
    for (l in 1:dim(results$obsdata[[j]]$data)[2]) {
      results$merged[[j]][[l]] <- matrix(NA, numppdsamples, length(results$obsdata[[j]]$time))
    }
  }
  
  pb <- txtProgressBar(max=length(results$ppd_sample_ix), style=3)
  for (i in 1:numppdsamples) {
    sample <- model$posterior$samples[,temperature_ix,results$ppd_sample_ix[i]]
    names(sample) <- model$variables
    simdata <- bcm3.cellpop.get.simulated.data(model, model$likelihood$experiments[[experiment_ix]]$name, sample)
    
    for (j in 1:length(results$obsdata)) {
      for (l in 1:dim(results$obsdata[[j]]$data)[2]) {
        vals <- simdata[[j]]$data[,l,1]
        vals[vals == -Inf] <- -2
        results$merged[[j]][[l]][i,] <- vals
      }
    }
    setTxtProgressBar(pb, i)
  }
  cat("\n")
  
  results$sds <- list()
  for (k in 1:length(results$obsdata)) {
    stdev_ix <- match(model$likelihood$experiments[[experiment_ix]]$data[[k]]$stdev, model$variables)
    if (is.na(stdev_ix)) {
      sds <- as.numeric(model$likelihood$experiments[[experiment_ix]]$data[[k]]$stdev)
    } else {
      sds <- 10^as.numeric(model$posterior$samples[stdev_ix, temperature_ix, results$ppd_sample_ix])
    }
    results$sds[[k]] <- sds
  }
  return(results)
}

calculate_Rsq <- function(bcm3, results, experiment_ix)
{
  model <- bcm3
  results$Rsq <- list()
  for (k in 1:length(results$obsdata)) {
    if (model$likelihood$experiments[[experiment_ix]]$data[[k]]$type == "duration") {
    } else {
      numppd <- nrow(results$merged[[k]][[cell_ix]])
      results$Rsq[[k]] <- matrix(NA, dim(results$obsdata[[k]]$data)[2], numppd)
      for (cell_ix in 1:dim(results$obsdata[[k]]$data)[2]) {
        #all_y <- c(results$obsdata[[k]]$data[,cell_ix,1], results$obsdata[[k+1]]$data[,cell_ix,1])
        all_y <- results$obsdata[[k]]$data[,cell_ix,1]
        avg_y <- mean(all_y, na.rm=T)
        SStot <- sum((all_y - avg_y)^2, na.rm=T)
        SSres <- rep(NA, numppd)
        for (i in 1:numppd) {
          res <- results$merged[[k]][[cell_ix]][i,] - results$obsdata[[k]]$data[,cell_ix,1]
          SSres[i] <- sum(res^2, na.rm=T)
        }
        Rsq <- 1 - (SSres / SStot)
        results$Rsq[[k]][cell_ix,] <- Rsq
      }
    }
  }
  return(results)
}

# A bit hacky but want to include the two replicate data points
calculate_Rsq_replicated_data <- function(bcm3, results, experiment_ix)
{
  model <- bcm3
  results$Rsq <- list()
  for (k in seq(1,length(results$obsdata), by=2)) {
    if (model$likelihood$experiments[[experiment_ix]]$data[[k]]$type == "duration") {
    } else {
      numppd <- nrow(results$merged[[k]][[1]])
      results$Rsq[[as.integer(k/2+0.5)]] <- matrix(NA, dim(results$obsdata[[k]]$data)[2], numppd)
      for (cell_ix in 1:dim(results$obsdata[[k]]$data)[2]) {
        all_y <- c(results$obsdata[[k]]$data[,cell_ix,1], results$obsdata[[k+1]]$data[,cell_ix,1])
        avg_y <- mean(all_y, na.rm=T)
        SStot <- sum((all_y - avg_y)^2, na.rm=T)
        SSres <- rep(NA, numppd)
        for (i in 1:numppd) {
          res <- c(results$merged[[k  ]][[cell_ix]][i,] - results$obsdata[[k  ]]$data[,cell_ix,1],
                   results$merged[[k+1]][[cell_ix]][i,] - results$obsdata[[k+1]]$data[,cell_ix,1])
          SSres[i] <- sum(res^2, na.rm=T)
        }
        Rsq <- 1 - (SSres / SStot)
        results$Rsq[[as.integer(k/2+0.5)]][cell_ix,] <- Rsq
      }
    }
  }
  return(results)
}

posterior_predictive_plot <- function(bcm3, results, experiment_ix, data_ix=1, cell_ixs=NULL, xlim=NULL, ylim=c(-0.2, 1.2), use_time="hours",
                                      ppd_point_symbol=4, ppd_color = "#00b25d",
                                      ...)
{
  model <- bcm3
  
  if (use_time == "hours") {
    t <- results$obsdata[[data_ix]]$time / 3600
  } else if (use_time == "minutes") {
    t <- results$obsdata[[data_ix]]$time / 60
  }
  if (is.null(cell_ixs)) {
    cell_ixs <- 1:dim(results$obsdata[[data_ix]]$data)[2]
  }
  for (cell_ix in cell_ixs) {
    lower <- rep(NA, length(results$obsdata[[data_ix]]$time))
    upper <- rep(NA, length(results$obsdata[[data_ix]]$time))
    for (i in 1:length(results$obsdata[[data_ix]]$time)) {
      x <- rnorm(5000, results$merged[[data_ix]][[cell_ix]][,i], results$sds[[data_ix]] * model$likelihood$experiments[[experiment_ix]]$data[[data_ix]]$stdev_multiplication_factor)
      lower[i] <- quantile(x, 0.05, na.rm=T)
      upper[i] <- quantile(x, 0.95, na.rm=T)
    }
    if (is.null(xlim)) {
      xlim <- range(t)
    }
    plot(t, results$obsdata[[data_ix]]$data[,cell_ix,1], pch=ppd_point_symbol, ylab=model$likelihood$experiments[[experiment_ix]]$data[[data_ix]]$data_name, xlab=paste("Time (", use_time, ")", sep=""), xlim=xlim, ylim=ylim, type='n', yaxt='n', ...)
    axis(2, at=seq(round(min(-0.7)),round(max(ylim)),by=0.5))
    rug(seq(ylim[1], ylim[2], by=0.1), side=2, ticksize=-0.03)
    if (use_time == "hours") {
      rug(seq(as.integer(min(t)), as.integer(max(t)), by=1), side=1, ticksize=-0.03)
    }
    have_plot <- which(!is.na(lower) & !is.na(upper))
    polygon(c(t[have_plot], rev(t[have_plot])), c(lower[have_plot], rev(upper[have_plot])), col=adjustcolor(ppd_color, alpha.f = 0.4), border=F)
    lines(t[have_plot], lower[have_plot], col=ppd_color, lwd=2)
    lines(t[have_plot], upper[have_plot], col=ppd_color, lwd=2)
    points(t, results$obsdata[[data_ix]]$data[,cell_ix,1], pch=ppd_point_symbol)
  }
}

# A bit hacky but want to get the two replicate data in the same plot
posterior_predictive_plot_replicate_data <- function(bcm3, results, experiment_ix, data_ix, cell_ix, ylim=c(-0.2, 1.2),
                                                     ppd_point_symbol=4, ppd_color = "#00b25d",
                                                     ...)
{
  model <- bcm3

  t <- results$obsdata[[data_ix]]$time / 3600
  
  lower <- rep(NA, length(results$obsdata[[data_ix]]$time))
  upper <- rep(NA, length(results$obsdata[[data_ix]]$time))
  for (i in 1:length(results$obsdata[[data_ix]]$time)) {
    x <- rnorm(5000, results$merged[[data_ix]][[cell_ix]][,i], results$sds[[data_ix]] * model$likelihood$experiments[[experiment_ix]]$data[[data_ix]]$stdev_multiplication_factor)
    lower[i] <- quantile(x, 0.05, na.rm=T)
    upper[i] <- quantile(x, 0.95, na.rm=T)
  }
  plot(t, results$obsdata[[data_ix]]$data[,cell_ix,1], pch=ppd_point_symbol, ylab=model$likelihood$experiments[[experiment_ix]]$data[[data_ix]]$data_name, xlab="Time (hours)", ylim=ylim, type='n', yaxt='n', ...)
  axis(2, at=seq(0,as.integer(max(ylim)),by=0.5))
  rug(seq(ylim[1], ylim[2], by=0.1), side=2, ticksize=-0.03)
  rug(seq(as.integer(min(t)), as.integer(max(t)), by=1), side=1, ticksize=-0.03)
  have_plot <- which(!is.na(lower) & !is.na(upper))
  polygon(c(t[have_plot], rev(t[have_plot])), c(lower[have_plot], rev(upper[have_plot])), col=adjustcolor(ppd_color, alpha.f = 0.4), border=F)
  lines(t[have_plot], lower[have_plot], col=ppd_color, lwd=2)
  lines(t[have_plot], upper[have_plot], col=ppd_color, lwd=2)
  points(t, results$obsdata[[data_ix]]$data[,cell_ix,1], pch=ppd_point_symbol)
  points(t, results$obsdata[[data_ix+1]]$data[,cell_ix,1], pch=ppd_point_symbol)
}

trajectory_heatmap <- function(bcm3, species_names, experiment_ix=1, temperature_ix=dim(bcm3$posterior$samples)[2], sample_ix=which.max(bcm3$posterior$lposterior[temperature_ix,]), range=c(0,1), draw_whole_population = T, draw_population_average = F, ...)
{
  sample <- bcm3$posterior$samples[,temperature_ix,sample_ix]
  
  traj <- bcm3.cellpop.get.simulated.trajectories(bcm3, bcm3$likelihood$experiments[[experiment_ix]]$name, sample)
  
  species_ix <- which(dimnames(traj$cells)[[1]] %in% species_names)
  mitogen_ix <- which(dimnames(traj$cells)[[1]] == "mitogenic_signal")
    
  plotdm <- t(apply(traj$cells[species_ix,,1:32,drop=F], c(2,3), sum))
  mitosis_values <- matrix(NA, nrow(plotdm), ncol(plotdm))
  parent_used <- rep(F, length(traj$parents))
  cell_remap <- rep(NA, length(traj$parents))
  cell_remap[1:32] <- 1:32
  new_cell_start_time <- rep(NA, 32)
  if (length(traj$parents) > 32) {
    for (i in 33:length(traj$parents)) {
      parent_ix <- traj$parents[i]+1
      if (parent_used[parent_ix]) {
        plotdm <- rbind(plotdm, t(apply(traj$cells[species_ix,,i,drop=F], c(2,3), sum)))
        mitosis_values <- rbind(mitosis_values, NA)
        cell_remap[i] <- nrow(plotdm)
      } else {
        parent_used[parent_ix] <- T
        time_ix <- which(!is.na(traj$cells[species_ix[1],,i]))
        plotdm[cell_remap[parent_ix],time_ix] <- apply(traj$cells[species_ix,time_ix,i,drop=F], c(2,3), sum)
        cell_remap[i] <- parent_ix
        if (is.na(new_cell_start_time[cell_remap[parent_ix]])) {
          new_cell_start_time[cell_remap[parent_ix]] <- time_ix[1]
          new_cell_start_time <- c(new_cell_start_time, time_ix[1])
        }
        
        mitosis_values[cell_remap[parent_ix],time_ix[1]] <- 100
        mitosis_values[cell_remap[parent_ix],time_ix[1]-1] <- 100
      }
    }
    ordering <- c(order(traj$cells[mitogen_ix,1,], decreasing = T)[1:32], 32+order(new_cell_start_time[33:nrow(plotdm)]))
    
    initial_cells <- new_cell_start_time[1:32]
    ordering <- c(intersect(order(initial_cells),which(!is.na(initial_cells))),
                  setdiff(order(traj$cells[mitogen_ix,1,1:32], decreasing = T), which(!is.na(initial_cells))),
                  32+order(new_cell_start_time[33:nrow(plotdm)]))
  } else {
    ordering <- order(traj$cells[mitogen_ix,1,], decreasing = T)
  }
  
  make_color_function <- function(breaks, colors, transparency = 0, space = "LAB", mitosis_value = 100)
  {
    colorRamp2_col_fun = colorRamp2(breaks, colors, transparency, space)
    new_col_fun <- function(x = NULL, return_rgb = FALSE, max_value = 1) {
      l_na = is.na(x)
      if(all(l_na)) {
        return(rep(NA, length(l_na)))
      }
      x2 = x[!l_na]
      colors <- character(length(x2))
      colors[x2 == mitosis_value] <- rep(rgb(1,0,0,1), sum(x2 == mitosis_value))
      colors[x2 != mitosis_value] <- colorRamp2_col_fun(x2[x2 != mitosis_value], return_rgb, max_value)
      
      res_col2 = character(length(x))
      res_col2[l_na] = NA
      res_col2[!l_na] = colors
      return(res_col2)
    }
    attr <- list(colorRamp2_col_fun = colorRamp2_col_fun, mitotis_value = mitosis_value, breaks=breaks)
    attributes(new_col_fun) = attr
    return(new_col_fun)
  }
  col_fun <- make_color_function(seq(range[1],range[2], len=25), viridis(25))
  
  plotdm_with_mitosis <- plotdm
  plotdm_with_mitosis[which(!is.na(mitosis_values))] <- mitosis_values[!is.na(mitosis_values)]
  
  if (draw_whole_population) {
    draw(Heatmap(plotdm_with_mitosis[ordering,], col=col_fun, cluster_rows = F, cluster_columns = F, ...))
  }
  
  if (draw_population_average) {
    population_average <- apply(plotdm, 2, mean, na.rm=T)
    draw(Heatmap(t(population_average), col=colorRamp2(seq(range[1],range[2], len=25), viridis(25)), cluster_rows = F, cluster_columns = F, ...))
  }
}

calculate_timing_log_loss <- function(model, timings, stdevfactor, experiment_ix=1, data_ix=1, numppdsamples = dim(model$posterior$samples)[3]/2, max_simulation_time = 86400)
{
  results <- list()
  
  sample_ix <- (dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3])
  if (numppdsamples == dim(model$posterior$samples)[3]/2) {
    results$ll_sample_ix <- sample_ix
  } else {
    results$ll_sample_ix <- sample(sample_ix, numppdsamples)
  }
  
  results$logloss <- matrix(NA, dim(model$posterior$samples)[2], numppdsamples)
  for (temp_ix in 1:nrow(results$logloss)) {
    cat("Temperature", temp_ix, "\n")
    pb <- txtProgressBar(min = 0, max = length(results$ll_sample_ix), style = 3)
    for (i in 1:length(results$ll_sample_ix)) {
      sample <- model$posterior$samples[,temp_ix,results$ll_sample_ix[i]]
      
      stdev_ix <- match(model$likelihood$experiments[[experiment_ix]]$data[[data_ix]]$stdev, model$variables)
      if (is.na(stdev_ix)) {
        sd <- as.numeric(model$likelihood$experiments[[experiment_ix]]$data[[data_ix]]$stdev) * stdevfactor
      } else {
        sd <- (10^as.numeric(sample[which(model$variables == "stdev")])) * stdevfactor
      }
      
      simdata <- bcm3.cellpop.get.simulated.data(model, model$likelihood$experiments[[experiment_ix]]$name, sample)
      times <- simdata[[data_ix]]$data[1,,1]
      times[is.na(times)] <- max_simulation_time
      
      distmat <- matrix(NA, length(times), length(timings))
      for (j in 1:length(times)) {
        for (k in 1:length(timings)) {
          distmat[j,k] <- (times[j] - timings[k])^2
        }
      }
      
      stopifnot(!is.na(distmat))
      distmat <- distmat - min(distmat)
      matching <- HungarianSolver(t(distmat))
      
      mapping <- matching$pairs[match(1:length(timings), matching$pairs[,1]),2]
      #plot(timings, times[mapping])
      logpdf <- dnorm(timings, times[mapping], sd, log=T)
      
      results$logloss[temp_ix, i] <- -sum(logpdf)
      
      setTxtProgressBar(pb, i)
    }
    cat("\n")
  }
  return(results)
}
