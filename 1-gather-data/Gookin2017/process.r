library(R.matlab)
library(pals)
library(ncdf4)
library(ComplexHeatmap)
library(circlize)

epitopes <- list.files("Fig_3_5/")
epitopes <- c(epitopes[2], epitopes[1], epitopes[3:9])

all_data <- list()
for (e in epitopes) {
  directory <- paste("Fig_3_5/", e, "/", sep="")
  files <- list.files(directory)
  all_wells <- NULL
  for (f in files) {
    res <- readMat(paste(directory, f, sep=""))
    
    if (f == files[1]) {
      all_wells <- list()
      all_wells$cdk2traces <- rbind(res$cdk2emergtraces, res$cdk2lowtraces, res$cdk2inctraces)
      if ("cdk2lowIF.cyto" %in% names(res)) {
        all_wells$if_measurements <- cbind(res$cdk2emergIF.cyto, res$cdk2lowIF.cyto, res$cdk2incIF.cyto)
      } else {
        all_wells$if_measurements <- cbind(res$cdk2emergIF, res$cdk2lowIF, res$cdk2incIF)
      }
      all_wells$frame_of_mitosis <- cbind(res$cdk2emergframeofmitosis, res$cdk2lowframeofmitosis, res$cdk2incframeofmitosis)
    } else {
      all_wells$cdk2traces <- rbind(all_wells$cdk2traces, res$cdk2emergtraces, res$cdk2lowtraces, res$cdk2inctraces)
      if ("cdk2lowIF.cyto" %in% names(res)) {
        all_wells$if_measurements <- cbind(all_wells$if_measurements, res$cdk2emergIF.cyto, res$cdk2lowIF.cyto, res$cdk2incIF.cyto)
      } else {
        all_wells$if_measurements <- cbind(all_wells$if_measurements, res$cdk2emergIF, res$cdk2lowIF, res$cdk2incIF)
      }
      all_wells$frame_of_mitosis <- cbind(all_wells$frame_of_mitosis, res$cdk2emergframeofmitosis, res$cdk2lowframeofmitosis, res$cdk2incframeofmitosis)
    }
  }
  
  # Some datasets have more than 120 timepoints, drop timepoints at the beginning to make all data structures the same
  ntimepoints <- ncol(all_wells$cdk2traces)
  if (ntimepoints > 120) {
    all_wells$cdk2traces <- all_wells$cdk2traces[,(ntimepoints-119):ntimepoints]
  }
  
  all_data[[e]] <- all_wells
}

# Make one huge list of all cells, to select a subset at random for model parameter inference
all_data_concatenated <- all_data[[1]]
all_data_concatenated$if_epitope <- rep(names(all_data)[1], ncol(all_data[[1]]$if_measurements))
for (i in 2:length(all_data)) {
  all_data_concatenated$cdk2traces <- rbind(all_data_concatenated$cdk2traces, all_data[[i]]$cdk2traces)
  all_data_concatenated$if_measurements <- cbind(all_data_concatenated$if_measurements, all_data[[i]]$if_measurements)
  all_data_concatenated$frame_of_mitosis <- cbind(all_data_concatenated$frame_of_mitosis, all_data[[i]]$frame_of_mitosis)
  all_data_concatenated$if_epitope <- c(all_data_concatenated$if_epitope, rep(names(all_data)[i], ncol(all_data[[i]]$if_measurements)))
}

split_cells_at_mitosis <- function(traces, frame_of_mitosis, if_measurements, if_epitope)
{
  split_into_parents <- matrix(NA, nrow(traces) * 3, ncol(traces))
  split_into_parents_IF <- rep(NA, nrow(traces) * 3)
  split_into_parents_IF_epitope <- rep(NA, nrow(traces) * 3)
  cell_ix <- 1
  out_cell_ix <- 1
  parents <- c()
  parent_cells <- c()
  while(cell_ix <= nrow(traces)) {
    have_two_children <- F
    
    if (cell_ix != nrow(traces)) {
      # Does the next cell have part of the trajectory exactly duplicated?
      diff <- traces[cell_ix+1,] - traces[cell_ix,]
      
      if (sum(diff == 0, na.rm=T) > 1) {
        shared_parent_trajectory <- traces[cell_ix,1:(frame_of_mitosis[cell_ix]-1)]
        if (frame_of_mitosis[cell_ix] <= 61) {
          # We're only using the second half of the timepoints, so just keep the two daughters as if they have no parents
          split_into_parents[out_cell_ix+0,(length(shared_parent_trajectory)+1):120] <- traces[cell_ix  ,frame_of_mitosis[cell_ix]:120]
          split_into_parents[out_cell_ix+1,(length(shared_parent_trajectory)+1):120] <- traces[cell_ix+1,frame_of_mitosis[cell_ix]:120]
          split_into_parents_IF[out_cell_ix+0] <- if_measurements[cell_ix  ]
          split_into_parents_IF[out_cell_ix+1] <- if_measurements[cell_ix+1]
          split_into_parents_IF_epitope[out_cell_ix+0] <- if_epitope[cell_ix  ]
          split_into_parents_IF_epitope[out_cell_ix+1] <- if_epitope[cell_ix+1]
          
          parents[out_cell_ix+0] <- NA
          parents[out_cell_ix+1] <- NA
          parent_cells <- c(parent_cells, out_cell_ix, out_cell_ix+1)
          
          cell_ix <- cell_ix + 2
          out_cell_ix <- out_cell_ix + 2
          have_two_children <- T
        } else {
          split_into_parents[out_cell_ix  ,1:length(shared_parent_trajectory)] <- shared_parent_trajectory
          split_into_parents[out_cell_ix+1,(length(shared_parent_trajectory)+1):120] <- traces[cell_ix  ,frame_of_mitosis[cell_ix]:120]
          split_into_parents[out_cell_ix+2,(length(shared_parent_trajectory)+1):120] <- traces[cell_ix+1,frame_of_mitosis[cell_ix]:120]
          split_into_parents_IF[out_cell_ix+1] <- if_measurements[cell_ix  ]
          split_into_parents_IF[out_cell_ix+2] <- if_measurements[cell_ix+1]
          split_into_parents_IF_epitope[out_cell_ix+1] <- if_epitope[cell_ix  ]
          split_into_parents_IF_epitope[out_cell_ix+2] <- if_epitope[cell_ix+1]
          
          parents[out_cell_ix] <- NA
          parents[out_cell_ix+1] <- out_cell_ix
          parents[out_cell_ix+2] <- out_cell_ix
          parent_cells <- c(parent_cells, out_cell_ix)
          
          cell_ix <- cell_ix + 2
          out_cell_ix <- out_cell_ix + 3
          have_two_children <- T
        }
      }
    }
    
    if (!have_two_children) {
      # Don't have a duplicate parent trajectory, just split this cell at the mitosis point if needed
      if (frame_of_mitosis[cell_ix] <= 61) {
        # We're only using the second half of the timepoints, so don't need to split this cell
        split_into_parents[out_cell_ix,] <- traces[cell_ix, ]
        split_into_parents_IF[out_cell_ix] <- if_measurements[cell_ix]
        split_into_parents_IF_epitope[out_cell_ix] <- if_epitope[cell_ix]
        
        parents[out_cell_ix] <- NA
        parent_cells <- c(parent_cells, out_cell_ix)
        
        cell_ix <- cell_ix + 1
        out_cell_ix <- out_cell_ix + 1
      } else {
        parent_trajectory <- traces[cell_ix,1:(frame_of_mitosis[cell_ix]-1)]
        split_into_parents[out_cell_ix  ,1:length(parent_trajectory)] <- parent_trajectory
        split_into_parents[out_cell_ix+1,(length(parent_trajectory)+1):120] <- traces[cell_ix, frame_of_mitosis[cell_ix]:120]
        split_into_parents_IF[out_cell_ix+1] <- if_measurements[cell_ix]
        split_into_parents_IF_epitope[out_cell_ix+1] <- if_epitope[cell_ix]
        
        parents[out_cell_ix] <- NA
        parents[out_cell_ix+1] <- out_cell_ix
        parent_cells <- c(parent_cells, out_cell_ix)
        
        cell_ix <- cell_ix + 1
        out_cell_ix <- out_cell_ix + 2
      }
    }
  }
  split_into_parents <- split_into_parents[1:(out_cell_ix-1),]
  split_into_parents_IF <- split_into_parents_IF[1:(out_cell_ix-1)]
  
  res <- list()
  res$split_into_parents <- split_into_parents
  res$split_into_parents_IF <- split_into_parents_IF
  res$split_into_parents_IF_epitope <- split_into_parents_IF_epitope
  res$parents <- parents
  res$which_cells_are_parents <- parent_cells
  return(res)
}

split_res <- split_cells_at_mitosis(all_data_concatenated$cdk2traces,
                                    all_data_concatenated$frame_of_mitosis,
                                    all_data_concatenated$if_measurements,
                                    all_data_concatenated$if_epitope)

# Heatmap(split_res$split_into_parents[,], col=cividis(30), cluster_columns = F, cluster_rows = F)
Heatmap(split_res$split_into_parents[sample(nrow(split_res$split_into_parents), 100),], col=cividis(30), cluster_columns = F, cluster_rows = F)

# Some cells had multiple mitosis, but these have not been marked
# Cut off the first 12 hours to get rid of this
second_half <- split_res$split_into_parents[,61:120]
Heatmap(second_half[sample(nrow(split_res$split_into_parents), 100),], col=cividis(30), cluster_columns = F, cluster_rows = F, use_raster = F)
Heatmap(second_half[1:10,], col=cividis(30), cluster_columns = F, cluster_rows = F, use_raster = F)

ncells <- 40
selected_ix <- sample(split_res$which_cells_are_parents, ncells)

# Bit convoluted but want to have parents & daughters next to each other
# union(selected_ix, which(split_res$parents %in% selected_ix))
selected_cells_with_daughters <- c()
for (i in 1:ncells) {
  selected_cells_with_daughters <- c(selected_cells_with_daughters, selected_ix[i], which(split_res$parents == selected_ix[i]))
}

pdf("selected_cells_heatmap.pdf", width=8, height=8)
Heatmap(log10(second_half[selected_cells_with_daughters,]), col=colorRamp2(seq(-0.5, 0.5, len=30), cividis(30)), cluster_columns = F, cluster_rows = F)
res <- dev.off()


dim1 <- ncdim_def("Gookin_random_subset/time", "seconds", (0:59) * 12 * 60)
dim2 <- ncdim_def("Gookin_random_subset/cell_id", "id", selected_cells_with_daughters)
dim3 <- ncdim_def("Gookin_random_subset/marker", "marker", 1)
parvar <- ncvar_def("Gookin_random_subset/parent", "id", dim2, prec="integer")
signalvar <- ncvar_def("Gookin_random_subset/CDK2_sensor", "signal", list(dim3, dim2, dim1), NA, prec="float")
signalvarlog <- ncvar_def("Gookin_random_subset/CDK2_sensor_log", "signal", list(dim3, dim2, dim1), NA, prec="float")

ncnew <- nc_create("gookin2017_random_subset.nc", list(parvar, signalvar, signalvarlog), force_v4=T)
ncvar_put(ncnew, parvar, split_res$parents[selected_cells_with_daughters])
for (i in 1:length(selected_cells_with_daughters)) {
  ncvar_put(ncnew, signalvar   ,       split_res$split_into_parents[selected_cells_with_daughters[i],61:120] , start=c(1,i,1), count=c(1, 1, 60))
  ncvar_put(ncnew, signalvarlog, log10(split_res$split_into_parents[selected_cells_with_daughters[i],61:120]), start=c(1,i,1), count=c(1, 1, 60))
}
nc_close(ncnew)



output_cell_ix <- 1:length(split_res$parents)

dim1 <- ncdim_def("Gookin_all_cells/time", "seconds", (0:59) * 12 * 60)
dim2 <- ncdim_def("Gookin_all_cells/cell_id", "id", output_cell_ix)
dim3 <- ncdim_def("Gookin_all_cells/marker", "marker", 1)
dim4 <- ncdim_def("Gookin_all_cells/timeIF", "seconds", 59 * 12 * 60)
parvar <- ncvar_def("Gookin_all_cells/parent", "id", dim2, prec="integer")
signalvar <- ncvar_def("Gookin_all_cells/CDK2_sensor", "signal", list(dim3, dim2, dim1), NA, prec="float")
signalvarlog <- ncvar_def("Gookin_all_cells/CDK2_sensor_log", "signal", list(dim3, dim2, dim1), NA, prec="float")
ifvar <- list()
for (e in epitopes) {
  ifvar[[e]] <- ncvar_def(paste("Gookin_all_cells/IF_", e, sep=""), "signal", list(dim3, dim2, dim4), NA, prec="float")
}

ncnew <- nc_create("gookin2017_all_cells.nc", c(list(parvar, signalvar, signalvarlog), ifvar), force_v4=T)
ncvar_put(ncnew, parvar, split_res$parents[output_cell_ix])
for (i in 1:length(output_cell_ix)) {
  ncvar_put(ncnew, signalvar   ,       split_res$split_into_parents[output_cell_ix[i],61:120] , start=c(1,i,1), count=c(1, 1, 60))
  ncvar_put(ncnew, signalvarlog, log10(split_res$split_into_parents[output_cell_ix[i],61:120]), start=c(1,i,1), count=c(1, 1, 60))
  epi <- split_res$split_into_parents_IF_epitope[output_cell_ix[i]]
  if (!is.na(epi)) {
    ncvar_put(ncnew, ifvar[[epi]], split_res$split_into_parents_IF[output_cell_ix[i]], start=c(1,i,1), count=c(1,1,1))
  }
}
nc_close(ncnew)










res <- all_wells

plot(1, type='n', xlim=c(0,30), ylim=c(0.2, 5), log='y', xlab='Time (hours)', ylab='CDK2 sensor')
for (i in 1:nrow(res$cdk2lowtraces)) {
  color <- brewer.set1(3)[1]
  lines((1:120)*(12/60), res$cdk2lowtraces[i,], col=color)
  mito_frame <- res$cdk2lowframeofmitosis[1,i]
  points(mito_frame*(12/60), res$cdk2lowtraces[i,mito_frame], col=color, pch=15)

  points(30, res$cdk2lowIF[1,i]/12, pch=1, col=color)
}
for (i in 1:nrow(res$cdk2emergtraces)) {
  color <- brewer.set1(3)[2]
  lines((1:120)*(12/60), res$cdk2emergtraces[i,], col=color)
  mito_frame <- res$cdk2emergframeofmitosis[1,i]
  points(mito_frame*(12/60), res$cdk2emergtraces[i,mito_frame], col=color, pch=15)

  points(30, res$cdk2emergIF[1,i]/12, pch=1, col=color)
}
for (i in 1:nrow(res$cdk2inctraces)) {
  color <- brewer.set1(3)[3]
  lines((1:120)*(12/60), res$cdk2inctraces[i,], col=color)
  mito_frame <- res$cdk2incframeofmitosis[1,i]
  points(mito_frame*(12/60), res$cdk2inctraces[i,mito_frame], col=color, pch=15)

  points(30, res$cdk2incIF[1,i]/12, pch=1, col=color)
}


plot(1, type='n', xlim=c(-20,15), ylim=c(0.2, 2.5), log='y', xlab='Time (hours)', ylab='CDK2 sensor')
for (i in 1:15) {
  color <- brewer.set1(3)[1]
  mito_frame <- res$cdk2lowframeofmitosis[1,i]
  lines((1:120-mito_frame)*(12/60), res$cdk2lowtraces[i,], col=color)
}
for (i in 1:10) {
  color <- brewer.set1(3)[2]
  mito_frame <- res$cdk2emergframeofmitosis[1,i]
  lines((1:120-mito_frame)*(12/60), res$cdk2emergtraces[i,], col=color)
}

for (i in 1:29) {
  color <- brewer.set1(3)[3]
  mito_frame <- res$cdk2incframeofmitosis[1,i]
  lines((1:120-mito_frame)*(12/60), res$cdk2inctraces[i,], col=color)
}


merged_cdk2traces <- rbind(res$cdk2lowtraces,
                           res$cdk2emergtraces,
                           res$cdk2inctraces)
merged_measurement <- c(res$cdk2lowIF.cyto,
                        res$cdk2emergIF.cyto,
                        res$cdk2incIF.cyto) / 12
merged_frame_of_mitosis <- c(res$cdk2lowframeofmitosis, res$cdk2emergframeofmitosis, res$cdk2incframeofmitosis)

Heatmap(merged_cdk2traces, col=colorRamp2(seq(0, 2, len=30), cividis(30)), cluster_columns = F, cluster_rows = F)
Heatmap(merged_cdk2traces[,61:120], col=colorRamp2(seq(0, 2, len=30), cividis(30)), cluster_columns = F, cluster_rows = F)







ix <- which(apply(is.na(merged_cdk2traces), 1, sum) == 0)
Heatmap(merged_cdk2traces[ix,], col=cividis(30), cluster_columns = F, cluster_rows = F)

selected_ix <- ix[c(2,4,6,7,8,9,10,11,12,13,15,17,18,20,21)]
Heatmap(merged_cdk2traces[selected_ix,], col=cividis(30), cluster_columns = F, cluster_rows = F)


split_into_parents <- matrix(NA, nrow(merged_cdk2traces) * 3, ncol(merged_cdk2traces))
split_into_parents_IF <- rep(NA, nrow(merged_cdk2traces) * 3)
cell_ix <- 1
out_cell_ix <- 1
parents <- c()
parent_cells <- c()
while(cell_ix <= nrow(merged_cdk2traces)) {
  have_two_children <- F

  if (cell_ix != nrow(merged_cdk2traces)) {
    # Does the next cell have part of the trajectory exactly duplicated?
    diff <- merged_cdk2traces[cell_ix+1,] - merged_cdk2traces[cell_ix,]

    if (sum(diff == 0, na.rm=T) > 1) {
      shared_parent_trajectory <- merged_cdk2traces[cell_ix,1:(merged_frame_of_mitosis[cell_ix]-1)]
      split_into_parents[out_cell_ix  ,1:length(shared_parent_trajectory)] <- shared_parent_trajectory
      split_into_parents[out_cell_ix+1,(length(shared_parent_trajectory)+1):120] <- merged_cdk2traces[cell_ix  ,merged_frame_of_mitosis[cell_ix]:120]
      split_into_parents[out_cell_ix+2,(length(shared_parent_trajectory)+1):120] <- merged_cdk2traces[cell_ix+1,merged_frame_of_mitosis[cell_ix]:120]
      split_into_parents_IF[out_cell_ix+1] <- merged_measurement[cell_ix  ]
      split_into_parents_IF[out_cell_ix+2] <- merged_measurement[cell_ix+1]

      parents[out_cell_ix] <- NA
      parents[out_cell_ix+1] <- out_cell_ix
      parents[out_cell_ix+2] <- out_cell_ix
      parent_cells <- c(parent_cells, out_cell_ix)

      cell_ix <- cell_ix + 2
      out_cell_ix <- out_cell_ix + 3
      have_two_children <- T
    }
  }

  if (!have_two_children) {
    # Don't have a duplicate parent trajectory, just split this cell at the mitosis point if needed
    if (merged_frame_of_mitosis[cell_ix] <= 60) {
      # We're only using the second half of the timepoints, so don't need to split this cell
      split_into_parents[out_cell_ix,] <- merged_cdk2traces[cell_ix, ]
      split_into_parents_IF[out_cell_ix] <- merged_measurement[cell_ix]

      parents[out_cell_ix] <- NA
      parent_cells <- c(parent_cells, out_cell_ix)

      cell_ix <- cell_ix + 1
      out_cell_ix <- out_cell_ix + 1
    } else {
      parent_trajectory <- merged_cdk2traces[cell_ix,1:(merged_frame_of_mitosis[cell_ix]-1)]
      split_into_parents[out_cell_ix  ,1:length(parent_trajectory)] <- parent_trajectory
      split_into_parents[out_cell_ix+1,(length(parent_trajectory)+1):120] <- merged_cdk2traces[cell_ix, merged_frame_of_mitosis[cell_ix]:120]
      split_into_parents_IF[out_cell_ix+1] <- merged_measurement[cell_ix+1]

      parents[out_cell_ix] <- NA
      parents[out_cell_ix+1] <- out_cell_ix
      parent_cells <- c(parent_cells, out_cell_ix)

      cell_ix <- cell_ix + 1
      out_cell_ix <- out_cell_ix + 2
    }
 }
}
split_into_parents <- split_into_parents[1:(out_cell_ix-1),]
split_into_parents_IF <- split_into_parents_IF[1:(out_cell_ix-1)]
Heatmap(split_into_parents, col=cividis(30), cluster_columns = F, cluster_rows = F)

# Some cells had multiple mitosis, but these have not been marked
# Cut off the first 12 hours to get rid of this
second_half <- split_into_parents[,61:120]
Heatmap(second_half, col=cividis(30), cluster_columns = F, cluster_rows = F)


ncells <- nrow(split_into_parents)

dim1 <- ncdim_def("Gookin_cycB/time", "seconds", (0:59) * 12 * 60)
dim2 <- ncdim_def("Gookin_cycB/cell_id", "id", 1:ncells)
dim3 <- ncdim_def("Gookin_cycB/marker", "marker", 1)
dim4 <- ncdim_def("Gookin_cycB/timeIF", "seconds", 59 * 12 * 60)
parvar <- ncvar_def("Gookin_cycB/parent", "id", dim2, prec="integer")
signalvar <- ncvar_def("Gookin_cycB/CDK2_sensor", "signal", list(dim3, dim2, dim1), NA, prec="float")
signalvarlog <- ncvar_def("Gookin_cycB/CDK2_sensor_log", "signal", list(dim3, dim2, dim1), NA, prec="float")
IFvar <- ncvar_def("Gookin_cycB/CycB_IF_cyto", "signal", list(dim3, dim2, dim4), NA, prec="float")

ncnew <- nc_create("gookin2017_allcells_2ndhalf.nc", list(parvar, signalvar, signalvarlog, IFvar), force_v4=T)
ncvar_put(ncnew, parvar, parents)
for (i in 1:ncells) {
  ncvar_put(ncnew, signalvar   ,       second_half[i,] , start=c(1,i,1), count=c(1, 1, 60))
  ncvar_put(ncnew, signalvarlog, log10(second_half[i,]), start=c(1,i,1), count=c(1, 1, 60))
  ncvar_put(ncnew, IFvar, split_into_parents_IF[i], start=c(1,i,1), count=c(1, 1, 1))
}
nc_close(ncnew)


# dim1 <- ncdim_def("Gookin_cycB/time", "seconds", (0:119) * 12 * 60)
# dim2 <- ncdim_def("Gookin_cycB/cell_id", "id", 1:2)
# dim3 <- ncdim_def("Gookin_cycB/marker", "marker", 1)
# parvar <- ncvar_def("Gookin_cycB/parent", "id", dim2, prec="integer")
# signalvar <- ncvar_def("Gookin_cycB/CDK2_sensor", "signal", list(dim3, dim2, dim1), NA, prec="float")
# signalvarlog <- ncvar_def("Gookin_cycB/CDK2_sensor_log", "signal", list(dim3, dim2, dim1), NA, prec="float")
#
# ncnew <- nc_create("gookin2017_1cell.nc", list(parvar, signalvar, signalvarlog), force_v4=T)
# ncvar_put(ncnew, parvar, c(NA, 1))
# ncvar_put(ncnew, signalvar, res$cdk2emergtraces[i,1:res$cdk2emergframeofmitosis[1,i]], start=c(1,1,1), count=c(1, 1, res$cdk2emergframeofmitosis[1,i]))
# ncvar_put(ncnew, signalvar, res$cdk2emergtraces[i,(res$cdk2emergframeofmitosis[1,i]+1):120], start=c(1,2,res$cdk2emergframeofmitosis[1,i]+1), count=c(1, 1, 120-res$cdk2emergframeofmitosis[1,i]))
# ncvar_put(ncnew, signalvarlog, log10(res$cdk2emergtraces[i,1:res$cdk2emergframeofmitosis[1,i]]), start=c(1,1,1), count=c(1, 1, res$cdk2emergframeofmitosis[1,i]))
# ncvar_put(ncnew, signalvarlog, log10(res$cdk2emergtraces[i,(res$cdk2emergframeofmitosis[1,i]+1):120]), start=c(1,2,res$cdk2emergframeofmitosis[1,i]+1), count=c(1, 1, 120-res$cdk2emergframeofmitosis[1,i]))
# nc_close(ncnew)
#
#
# ncells <- 4
#
# dim1 <- ncdim_def("Gookin_cycB/time", "seconds", (0:119) * 12 * 60)
# dim2 <- ncdim_def("Gookin_cycB/cell_id", "id", 1:(ncells*2))
# dim3 <- ncdim_def("Gookin_cycB/marker", "marker", 1)
# parvar <- ncvar_def("Gookin_cycB/parent", "id", dim2, prec="integer")
# signalvar <- ncvar_def("Gookin_cycB/CDK2_sensor", "signal", list(dim3, dim2, dim1), NA, prec="float")
# signalvarlog <- ncvar_def("Gookin_cycB/CDK2_sensor_log", "signal", list(dim3, dim2, dim1), NA, prec="float")
#
# ncnew <- nc_create("gookin2017_4cells.nc", list(parvar, signalvar, signalvarlog), force_v4=T)
# ncvar_put(ncnew, parvar, c(rbind(rep(NA, ncells),seq(1, ncells*2, by=2))))
# for (i in 1:ncells) {
#   ncvar_put(ncnew, signalvar, res$cdk2emergtraces[i,1:res$cdk2emergframeofmitosis[1,i]], start=c(1,(i*2)-1,1), count=c(1, 1, res$cdk2emergframeofmitosis[1,i]))
#   ncvar_put(ncnew, signalvar, res$cdk2emergtraces[i,(res$cdk2emergframeofmitosis[1,i]+1):120], start=c(1,i*2,res$cdk2emergframeofmitosis[1,i]+1), count=c(1, 1, 120-res$cdk2emergframeofmitosis[1,i]))
#   ncvar_put(ncnew, signalvarlog, log10(res$cdk2emergtraces[i,1:res$cdk2emergframeofmitosis[1,i]]), start=c(1,(i*2)-1,1), count=c(1, 1, res$cdk2emergframeofmitosis[1,i]))
#   ncvar_put(ncnew, signalvarlog, log10(res$cdk2emergtraces[i,(res$cdk2emergframeofmitosis[1,i]+1):120]), start=c(1,i*2,res$cdk2emergframeofmitosis[1,i]+1), count=c(1, 1, 120-res$cdk2emergframeofmitosis[1,i]))
# }
# nc_close(ncnew)
#
#
# i <- 1
# res$cdk2emergframeofmitosis[1,i]
#
# ncells <- length(selected_ix)
#
# dim1 <- ncdim_def("Gookin_cycB/time", "seconds", (0:59) * 12 * 60)
# dim2 <- ncdim_def("Gookin_cycB/cell_id", "id", 1:(ncells*2))
# dim3 <- ncdim_def("Gookin_cycB/marker", "marker", 1)
# dim4 <- ncdim_def("Gookin_cycB/timeIF", "seconds", 59 * 12 * 60)
# parvar <- ncvar_def("Gookin_cycB/parent", "id", dim2, prec="integer")
# signalvar <- ncvar_def("Gookin_cycB/CDK2_sensor", "signal", list(dim3, dim2, dim1), NA, prec="float")
# signalvarlog <- ncvar_def("Gookin_cycB/CDK2_sensor_log", "signal", list(dim3, dim2, dim1), NA, prec="float")
# IFvar <- ncvar_def("Gookin_cycB/CycB_IF_cyto", "signal", list(dim3, dim2, dim4), NA, prec="float")
#
# ncnew <- nc_create("gookin2017_15cells_2ndhalf.nc", list(parvar, signalvar, signalvarlog, IFvar), force_v4=T)
# ncvar_put(ncnew, parvar, c(rbind(rep(NA, ncells),seq(1, ncells*2, by=2))))
# for (i in 1:ncells) {
#   x1 <- merged_cdk2traces[selected_ix[i],1:merged_frame_of_mitosis[selected_ix[i]]]
#   x2 <- merged_cdk2traces[selected_ix[i],(merged_frame_of_mitosis[selected_ix[i]]+1):120]
#
#   ncvar_put(ncnew, signalvar, x1[-c(1:60)], start=c(1,(i*2)-1,1), count=c(1, 1, length(x1)-60))
#   ncvar_put(ncnew, signalvar, x2, start=c(1,i*2,length(x1)-60+1), count=c(1, 1, length(x2)))
#
#   ncvar_put(ncnew, signalvarlog, log10(x1[-c(1:60)]), start=c(1,(i*2)-1, 1              ), count=c(1, 1, length(x1)-60))
#   ncvar_put(ncnew, signalvarlog, log10(x2          ), start=c(1, i*2,    length(x1)-60+1), count=c(1, 1, length(x2)  ))
#
#   ncvar_put(ncnew, IFvar, merged_measurement[selected_ix[i]], start=c(1,i*2,1), count=c(1, 1, 1))
# }
# nc_close(ncnew)


