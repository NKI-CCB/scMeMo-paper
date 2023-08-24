source("../extract_data_functions.r")

ncells <- 53
fig1f <- list()
for (i in 1:ncells) {
  cat(i, "\n")
  fig1f[[i]] <- extract_data_from_SVG("Figure1F.svg", seq(0, 1200, by=200), seq(120,0, by=-20), point_layer_name = paste("cell", i, sep=""))
}

# Cut off everything before mitosis
for (i in 1:ncells) {
  ix <- fig1f[[i]]$x >= 0
  fig1f[[i]] <- fig1f[[i]][ix,]
}

# Scale between 0 and 1
for (i in 1:ncells) {
  fig1f[[i]]$y <- fig1f[[i]]$y / 80
}

# Merge time points
timepoints <- seq(0, 1200, by=12)
max_diff <- 0
for (i in 1:ncells) {
  for (j in 1:length(fig1f[[i]]$x)) {
    ix <- which.min(abs(timepoints - fig1f[[i]]$x[j]))
    max_diff <- max(abs(timepoints[ix] - fig1f[[i]]$x[j]), max_diff)
    fig1f[[i]]$x[j] <- timepoints[ix]
  }
}

plot(fig1f[[1]]$x, fig1f[[1]]$y, ylim=c(0,1.2), type='l')
iqrs <- rep(NA, length(fig1f))
for (i in 1:ncells) {
  lines(fig1f[[i]]$x, fig1f[[i]]$y)
  iqrs[i] <- diff(as.numeric(quantile(fig1f[[i]]$y, c(0.25, 0.75), na.rm=T)))
}
mean(iqrs)
mean(iqrs[iqrs > 0.2])

#selected_cell_lines <- sample(1:length(fig1f), 10)
# Forgot to store the seed, but these are the indices that were sampled
#selected_cell_lines <- c(26,40,20,48,1,10,31,2,13,42)

set.seed(1)
selected_cell_lines <- sample(1:length(fig1f), 16)

plot(fig1f[[1]]$x, fig1f[[1]]$y, ylim=c(0,1.2), type='n')
for (i in selected_cell_lines) {
  lines(fig1f[[i]]$x, fig1f[[i]]$y)
}

library(ncdf4)

dim1 <- ncdim_def("Barr_fig1f/time", "seconds", timepoints * 60)
dim2 <- ncdim_def("Barr_fig1f/cells", "id", 1:length(selected_cell_lines))
var1 <- ncvar_def("Barr_fig1f/p21_gfp", "signal", list(dim2, dim1), NA)

ncnew <- nc_create("Barr2017.nc", list(var1), force_v4=T)
for (i in 1:length(selected_cell_lines)) {
  cell_ix <- selected_cell_lines[i]
  time_ix <- match(fig1f[[cell_ix]]$x, timepoints)
  values <- rep(NA, length(timepoints))
  values[time_ix] <- fig1f[[cell_ix]]$y
  ncvar_put(ncnew, var1, values, start=c(i,1), count=c(1, length(values)))
}
nc_close(ncnew)

