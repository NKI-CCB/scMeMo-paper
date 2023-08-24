library(ncdf4)
source("../extract_data_functions.r")

ncells <- 17
fig1a <- list()
for (i in 1:ncells) {
  cat(i, "\n")
  fig1a[[i]] <- extract_data_from_SVG("Fig1a.svg", seq(-5, 10, by=1), seq(1,0, by=-0.1), point_layer_name = paste("cell", i, sep=""))
}
plot(fig1a[[1]]$x, fig1a[[1]]$y, ylim=c(0,1), type='l')
iqrs <- rep(NA, length(fig1a))
for (i in 1:ncells) {
  lines(fig1a[[i]]$x, fig1a[[i]]$y)
  iqrs[i] <- diff(as.numeric(quantile(fig1a[[i]]$y, c(0.25, 0.75), na.rm=T)))
}
mean(iqrs)

ncells <- 15 # 16th cell doesn't parse well.. dashed line?
fig1b <- list()
for (i in 1:ncells) {
  cat(i, "\n")
  fig1b[[i]] <- extract_data_from_SVG("Fig1b.svg", seq(-5, 10, by=1), seq(1,0, by=-0.1), point_layer_name = paste("cell", i, sep=""))
}
plot(fig1b[[1]]$x, fig1b[[1]]$y, ylim=c(0,1), type='l')
iqrs <- rep(NA, length(fig1b))
for (i in 1:ncells) {
  lines(fig1b[[i]]$x, fig1b[[i]]$y)
  iqrs[i] <- diff(as.numeric(quantile(fig1b[[i]]$y, c(0.25, 0.75), na.rm=T)))
}
mean(iqrs)

ncells <- 11
fig1c <- list()
for (i in 1:ncells) {
  cat(i, "\n")
  fig1c[[i]] <- extract_data_from_SVG("Fig1c.svg", seq(-5, 10, by=1), seq(1,0, by=-0.1), point_layer_name = paste("cell", i, sep=""))
}
# The outlier cell messes up the relative values, scale everything by 0.6; and a factor 2 because cycA is modelled from [0-2]
for (i in 1:ncells) {
  fig1c[[i]]$y <- 2.0 * fig1c[[i]]$y / 0.6
}
plot(fig1c[[1]]$x, fig1c[[1]]$y, ylim=c(0,1), type='l')
iqrs <- rep(NA, length(fig1c))
for (i in 1:ncells) {
  lines(fig1c[[i]]$x, fig1c[[i]]$y)
  iqrs[i] <- diff(as.numeric(quantile(fig1c[[i]]$y, c(0.25, 0.75), na.rm=T)))
}
mean(iqrs)

ncells <- 9
fig1d <- list()
for (i in 1:ncells) {
  cat(i, "\n")
  fig1d[[i]] <- extract_data_from_SVG("Fig1d.svg", seq(-5, 10, by=1), seq(3.4, 0.2, by=-0.2), point_layer_name = paste("cell", i, sep=""))
}
plot(fig1d[[1]]$x, fig1d[[1]]$y, ylim=c(0,4), type='l')
for (i in 1:ncells) {
  lines(fig1d[[i]]$x, fig1d[[i]]$y)
}
plot(fig1d[[1]]$x, log10(fig1d[[1]]$y), ylim=c(-0.5,0.5), type='l')
iqrs <- rep(NA, length(fig1d))
for (i in 1:ncells) {
  first_two_frames <- head(which(!is.na(fig1d[[i]]$y)), n=2)
  fig1d[[i]]$y[first_two_frames] <- NA
  lines(fig1d[[i]]$x, log10(fig1d[[i]]$y))
  iqrs[i] <- diff(as.numeric(quantile(log10(fig1d[[i]]$y), c(0.25, 0.75), na.rm=T)))
}
mean(iqrs)

ncells <- 10
#selected_cells_p27 <- sample(length(fig1b), ncells)
#selected_cells_cycE <- sample(length(fig1b), ncells)
#selected_cells_cycA <- sample(length(fig1c), ncells)
# Forgot to store the seed, but these are the indices that were sampled
selected_cells_p27 <- c(5,14,11,15,6,13,8,12,7,9)
selected_cells_cycE <- c(11,8,9,13,7,15,12,14,4,10)
selected_cells_cycA <- c(3,8,7,6,4,5,2,9,10,11)
selected_cells_cdk2sensor <- 1:length(fig1d) # Only 9 traces available, use them all

dim1 <- ncdim_def("Barr_fig1/time", "seconds", fig1c[[1]]$x * 60 * 60)
dim2 <- ncdim_def("Barr_fig1/cells", "id", 1:ncells)
var1 <- ncvar_def("Barr_fig1/p27_gfp", "signal", list(dim2, dim1), NA)
var2 <- ncvar_def("Barr_fig1/CyclinE_gfp", "signal", list(dim2, dim1), NA)
var3 <- ncvar_def("Barr_fig1/CyclinA_gfp", "signal", list(dim2, dim1), NA)
var4 <- ncvar_def("Barr_fig1/CDK2_sensor", "signal", list(dim2, dim1), NA)
var5 <- ncvar_def("Barr_fig1/CDK2_sensor_log", "signal", list(dim2, dim1), NA)
ncnew <- nc_create("Barr2016.nc", list(var1, var2, var3, var4, var5), force_v4=T)
for (i in 1:length(selected_cells_p27)) {
  ncvar_put(ncnew, var1, fig1a[[selected_cells_p27[i]]]$y, start=c(i,1), count=c(1, length(fig1a[[selected_cells_p27[i]]]$y)))
}
for (i in 1:length(selected_cells_cycE)) {
  ncvar_put(ncnew, var2, fig1b[[selected_cells_cycE[i]]]$y, start=c(i,1), count=c(1, length(fig1b[[selected_cells_cycE[i]]]$y)))
}
for (i in 1:length(selected_cells_cycA)) {
  ncvar_put(ncnew, var3, fig1c[[selected_cells_cycA[i]]]$y, start=c(i,1), count=c(1, length(fig1c[[selected_cells_cycA[i]]]$y)))
}
for (i in 1:length(selected_cells_cdk2sensor)) {
  ncvar_put(ncnew, var4, fig1d[[selected_cells_cdk2sensor[i]]]$y, start=c(i,1), count=c(1, length(fig1d[[selected_cells_cdk2sensor[i]]]$y)))
  ncvar_put(ncnew, var5, log10(fig1d[[selected_cells_cdk2sensor[i]]]$y), start=c(i,1), count=c(1, length(fig1d[[selected_cells_cdk2sensor[i]]]$y)))
}
nc_close(ncnew)


pdf("barr_2016_data_plot.pdf", width=120/25.4, height=40/25.4)
par(mfrow=c(1,3))
par(cex=0.5, mgp=c(2.0,0.8,0))
plot(0, xlim=c(-5, 10), ylim=c(0, 1), type='n', xlab="Time (hours)", ylab="Nuclear fluorescence (a.u.)", main="p27")
for (i in 1:length(selected_cells_p27)) {
  lines(fig1a[[1]]$x, fig1a[[selected_cells_p27[i]]]$y)
}
plot(0, xlim=c(-5, 10), ylim=c(0, 1), type='n', xlab="Time (hours)", ylab="Nuclear fluorescence (a.u.)", main="cyclin E")
for (i in 1:length(selected_cells_cycE)) {
  lines(fig1b[[1]]$x, fig1b[[selected_cells_cycE[i]]]$y)
}
plot(0, xlim=c(-5, 10), ylim=c(0, 2), type='n', xlab="Time (hours)", ylab="Nuclear fluorescence (a.u.)", main="cyclin A")
for (i in 1:length(selected_cells_cycA)) {
  lines(fig1c[[1]]$x, fig1c[[selected_cells_cycA[i]]]$y)
}
res <- dev.off()
