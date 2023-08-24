source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))

model <- bcm3.load.results("../2-mitosis-model",
                           output_folder="output_gavetpines",
                           prior_file="prior_gavetpines.xml",
                           likelihood_file="likelihood_gavetpines.xml")
model <- bcm3.init.cpp(model)

temperature_ix <- dim(model$posterior$samples)[2]
MAP_sample_ix <- which.max((model$posterior$lprior + model$posterior$llikelihood)[temperature_ix,])
num_uncertainty_samples <- 200
ppd_sample_ix <- sample((dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3]), num_uncertainty_samples)

active_cyclinB_ix <- 11
assembled_spindle_ix <- 15
num_cells <- 16

colors <- c("#4053d3","#ddb310","#b51d14","#00beff","#fb49b0","#00b25d","#cacaca")
linewidth <- 1.0

dir.create("figures", showWarnings = F)

pdf("figures/simulations.pdf", width=160/(25.4*0.75), height=50/(25.4*0.75))

par(mfcol=c(1,3))
par(mar=c(5,6,4,1))

sample <- model$posterior$samples[,temperature_ix,MAP_sample_ix]
names(sample) <- model$variables
sample["spindle_assembly_variability"] <- 1.0
sample["FOXM1_synthesis_variability"] <- 0.01
sample["mitotic_entry_variability"] <- 0.01
traj <- bcm3.cellpop.get.simulated.trajectories(model, model$likelihood$experiments[[1]]$name, sample)

use_time_ix <- which(!is.na(traj$cells[active_cyclinB_ix,,1]))
t <- (traj$time[use_time_ix] - traj$time[use_time_ix[1]])/60
plot(t, traj$cells[active_cyclinB_ix,use_time_ix,1], type='l', col=colors[1], lwd=linewidth, ylim=c(-0.1,1.1), xlab='Time (minutes)', ylab='Concentration\n(arbitrary units)', xaxt='n', yaxt='n')
axis(1, at=c(0,100,200,300), tck=-0.02)
rug(c(50,150,250), ticksize=-0.02, side=1)
axis(2, at=c(0,1), labels=format(c(0.0, 1.0), nsmall=1), tck=-0.02)
rug(c(0.2,0.4,0.6,0.8), ticksize=-0.02, side=2)
lines(t, traj$cells[assembled_spindle_ix,use_time_ix,1], type='l', col=colors[2], lwd=linewidth)

plot(t, traj$cells[active_cyclinB_ix,use_time_ix,1], type='l', col=colors[1], lwd=linewidth, ylim=c(-0.1,1.1), xlab='Time (minutes)', ylab='Concentration\n(arbitrary units)', xaxt='n', yaxt='n')
axis(1, at=c(0,100,200,300), tck=-0.02)
rug(c(50,150,250), ticksize=-0.02, side=1)
axis(2, at=c(0,1), labels=format(c(0.0, 1.0), nsmall=1), tck=-0.02)
rug(c(0.2,0.4,0.6,0.8), ticksize=-0.02, side=2)
for (i in 1:num_cells) {
  lines(t, traj$cells[active_cyclinB_ix,use_time_ix,i], type='l', col=colors[1], lwd=linewidth)
  lines(t, traj$cells[assembled_spindle_ix,use_time_ix,i], type='l', col=colors[2], lwd=linewidth)
}

sample <- model$posterior$samples[,temperature_ix,MAP_sample_ix]
names(sample) <- model$variables
sample["spindle_assembly_variability"] <- 0.01
sample["FOXM1_synthesis_variability"] <- 0.01
sample["mitotic_entry_variability"] <- 1.0
traj <- bcm3.cellpop.get.simulated.trajectories(model, model$likelihood$experiments[[1]]$name, sample)

plot(t, traj$cells[active_cyclinB_ix,use_time_ix,1], type='l', col=colors[1], lwd=linewidth, ylim=c(-0.1,1.1), xlab='Time (minutes)', ylab='Concentration\n(arbitrary units)', xaxt='n', yaxt='n')
axis(1, at=c(0,100,200,300), tck=-0.02)
rug(c(50,150,250), ticksize=-0.02, side=1)
axis(2, at=c(0,1), labels=format(c(0.0, 1.0), nsmall=1), tck=-0.02)
rug(c(0.2,0.4,0.6,0.8), ticksize=-0.02, side=2)
for (i in 1:num_cells) {
  lines(t, traj$cells[active_cyclinB_ix,use_time_ix,i], type='l', col=colors[1], lwd=linewidth)
  lines(t, traj$cells[assembled_spindle_ix,use_time_ix,i], type='l', col=colors[2], lwd=linewidth)
}

res <- dev.off()



pdf("figures/simulations_uncertainty.pdf", width=160/(25.4*0.75), height=50/(25.4*0.75))

par(mfcol=c(1,3))
par(mar=c(5,6,4,1))

species <- list()
species_all <- list()

pb <- txtProgressBar(max=length(ppd_sample_ix), style=3)
for (j in 1:length(ppd_sample_ix)) {
  setTxtProgressBar(pb, j)
  sample <- model$posterior$samples[,temperature_ix,ppd_sample_ix[j]]
  names(sample) <- model$variables
  sample["spindle_assembly_variability"] <- 1.0
  sample["FOXM1_synthesis_variability"] <- 0.01
  sample["mitotic_entry_variability"] <- 0.01
  traj <- bcm3.cellpop.get.simulated.trajectories(model, model$likelihood$experiments[[1]]$name, sample)
  
  if (j == 1) {
    species[[1]] <- traj$cells[active_cyclinB_ix,,1]
    species[[2]] <- traj$cells[assembled_spindle_ix,,1]
    
    for (k in 1:num_cells) {
      species_all[[1]] <- list()
      species_all[[1]][[k]] <- traj$cells[active_cyclinB_ix,,k]
      species_all[[2]] <- list()
      species_all[[2]][[k]] <- traj$cells[assembled_spindle_ix,,k]
    }
  } else {
    species[[1]] <- rbind(species[[1]], traj$cells[active_cyclinB_ix,,1])
    species[[2]] <- rbind(species[[2]], traj$cells[assembled_spindle_ix,,1])
    
    for (k in 1:num_cells) {
      species_all[[1]][[k]] <- rbind(species_all[[1]][[k]], traj$cells[active_cyclinB_ix,,k])
      species_all[[2]][[k]] <- rbind(species_all[[2]][[k]], traj$cells[assembled_spindle_ix,,k])
    }
  }
}

use_time_ix <- which(!is.na(traj$cells[active_cyclinB_ix,,1]))
t <- (traj$time[use_time_ix] - traj$time[use_time_ix[1]]) / 60
plot(t, traj$cells[active_cyclinB_ix,use_time_ix,1], type='n', ylim=c(-0.1,1.1), xlab='Time (minutes)', ylab='Concentration\n(arbitrary units)', xaxt='n', yaxt='n')
axis(1, at=c(0,100,200,300), tck=-0.02)
rug(c(50,150,250), ticksize=-0.02, side=1)
axis(2, at=c(0,1), labels=format(c(0.0, 1.0), nsmall=1), tck=-0.02)
rug(c(0.2,0.4,0.6,0.8), ticksize=-0.02, side=2)

lower <- apply(species[[1]], 2, quantile, probs=0.05, na.rm=T)
upper <- apply(species[[1]], 2, quantile, probs=0.95, na.rm=T)
lines(t, lower[use_time_ix], col=colors[1], lwd=0.5)
lines(t, upper[use_time_ix], col=colors[1], lwd=0.5)
coloralpha <- rgb(col2rgb(colors[1])[1,1],col2rgb(colors[1])[2,1],col2rgb(colors[1])[3,1], 0.1*255, max=255)
polygon(c(t, rev(t)), c(lower[use_time_ix], rev(upper[use_time_ix])), border=F, col=coloralpha)

lower <- apply(species[[2]], 2, quantile, probs=0.05, na.rm=T)
upper <- apply(species[[2]], 2, quantile, probs=0.95, na.rm=T)
lines(t, lower[use_time_ix], col=colors[2], lwd=0.5)
lines(t, upper[use_time_ix], col=colors[2], lwd=0.5)
coloralpha <- rgb(col2rgb(colors[2])[1,1],col2rgb(colors[2])[2,1],col2rgb(colors[2])[3,1], 0.1*255, max=255)
polygon(c(t, rev(t)), c(lower[use_time_ix], rev(upper[use_time_ix])), border=F, col=coloralpha)


plot(t, traj$cells[1,use_time_ix,1], type='n', col=colors[1], lwd=linewidth, ylim=c(-0.1,1.1), xlab='Time (minutes)', ylab='Concentration\n(arbitrary units)', xaxt='n', yaxt='n')
axis(1, at=c(0,100,200,300), tck=-0.02)
rug(c(50,150,250), ticksize=-0.02, side=1)
axis(2, at=c(0,1), labels=format(c(0.0, 1.0), nsmall=1), tck=-0.02)
rug(c(0.2,0.4,0.6,0.8), ticksize=-0.02, side=2)
for (k in 1:num_cells) {
  lower <- apply(species_all[[1]][[k]], 2, quantile, probs=0.05, na.rm=T)
  upper <- apply(species_all[[1]][[k]], 2, quantile, probs=0.95, na.rm=T)
  lines(t, lower[use_time_ix], col=colors[1], lwd=0.5)
  lines(t, upper[use_time_ix], col=colors[1], lwd=0.5)
  coloralpha <- rgb(col2rgb(colors[1])[1,1],col2rgb(colors[1])[2,1],col2rgb(colors[1])[3,1], 0.1*255, max=255)
  polygon(c(t, rev(t)), c(lower[use_time_ix], rev(upper[use_time_ix])), border=F, col=coloralpha)
}
for (k in 1:num_cells) {
  lower <- apply(species_all[[2]][[k]], 2, quantile, probs=0.05, na.rm=T)
  upper <- apply(species_all[[2]][[k]], 2, quantile, probs=0.95, na.rm=T)
  lines(t, lower[use_time_ix], col=colors[2], lwd=0.5)
  lines(t, upper[use_time_ix], col=colors[2], lwd=0.5)
  coloralpha <- rgb(col2rgb(colors[2])[1,1],col2rgb(colors[2])[2,1],col2rgb(colors[2])[3,1], 0.1*255, max=255)
  polygon(c(t, rev(t)), c(lower[use_time_ix], rev(upper[use_time_ix])), border=F, col=coloralpha)
}

species <- list()
species_all <- list()

pb <- txtProgressBar(max=length(ppd_sample_ix), style=3)
for (j in 1:length(ppd_sample_ix)) {
  setTxtProgressBar(pb, j)
  sample <- model$posterior$samples[,temperature_ix,ppd_sample_ix[j]]
  names(sample) <- model$variables
  sample["spindle_assembly_variability"] <- 0.01
  sample["FOXM1_synthesis_variability"] <- 0.01
  sample["mitotic_entry_variability"] <- 1.0
  traj <- bcm3.cellpop.get.simulated.trajectories(model, model$likelihood$experiments[[1]]$name, sample)
  
  if (j == 1) {
    species[[1]] <- traj$cells[active_cyclinB_ix,,1]
    species[[2]] <- traj$cells[assembled_spindle_ix,,1]
    
    for (k in 1:num_cells) {
      species_all[[1]] <- list()
      species_all[[1]][[k]] <- traj$cells[active_cyclinB_ix,,k]
      species_all[[2]] <- list()
      species_all[[2]][[k]] <- traj$cells[assembled_spindle_ix,,k]
    }
  } else {
    species[[1]] <- rbind(species[[1]], traj$cells[active_cyclinB_ix,,1])
    species[[2]] <- rbind(species[[2]], traj$cells[assembled_spindle_ix,,1])
    
    for (k in 1:num_cells) {
      species_all[[1]][[k]] <- rbind(species_all[[1]][[k]], traj$cells[active_cyclinB_ix,,k])
      species_all[[2]][[k]] <- rbind(species_all[[2]][[k]], traj$cells[assembled_spindle_ix,,k])
    }
  }
}

plot(t, traj$cells[1,use_time_ix,1], type='n', col=colors[1], lwd=linewidth, ylim=c(-0.1,1.1), xlab='Time (minutes)', ylab='Concentration\n(arbitrary units)', xaxt='n', yaxt='n')
axis(1, at=c(0,100,200,300), tck=-0.02)
rug(c(50,150,250), ticksize=-0.02, side=1)
axis(2, at=c(0,1), labels=format(c(0.0, 1.0), nsmall=1), tck=-0.02)
rug(c(0.2,0.4,0.6,0.8), ticksize=-0.02, side=2)
for (k in 1:num_cells) {
  lower <- apply(species_all[[1]][[k]], 2, quantile, probs=0.05, na.rm=T)
  upper <- apply(species_all[[1]][[k]], 2, quantile, probs=0.95, na.rm=T)
  lines(t, lower[use_time_ix], col=colors[1], lwd=0.5)
  lines(t, upper[use_time_ix], col=colors[1], lwd=0.5)
  coloralpha <- rgb(col2rgb(colors[1])[1,1],col2rgb(colors[1])[2,1],col2rgb(colors[1])[3,1], 0.1*255, max=255)
  polygon(c(t, rev(t)), c(lower[use_time_ix], rev(upper[use_time_ix])), border=F, col=coloralpha)
}
  
for (k in 1:num_cells) {
  lower <- apply(species_all[[2]][[k]], 2, quantile, probs=0.05, na.rm=T)
  upper <- apply(species_all[[2]][[k]], 2, quantile, probs=0.95, na.rm=T)
  lines(t, lower[use_time_ix], col=colors[2], lwd=0.5)
  lines(t, upper[use_time_ix], col=colors[2], lwd=0.5)
  coloralpha <- rgb(col2rgb(colors[2])[1,1],col2rgb(colors[2])[2,1],col2rgb(colors[2])[3,1], 0.1*255, max=255)
  polygon(c(t, rev(t)), c(lower[use_time_ix], rev(upper[use_time_ix])), border=F, col=coloralpha)
}

res <- dev.off()


model <- bcm3.release.cpp(model)
