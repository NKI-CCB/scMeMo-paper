df <- read.delim("fig1a_quantification.tsv", stringsAsFactors = F, check.names = F)

time <- as.numeric(df[,1]) * 60 * 60
e2f1 <- as.numeric(df[,2]) / max(as.numeric(df[,2]))
cycf <- as.numeric(df[,3]) / max(as.numeric(df[,3]))

plot(time, e2f1, pch=19, ylim=c(0,1))
points(time, cycf, col='red', pch=19)

diff(as.numeric(quantile(e2f1, c(0.25, 0.75), na.rm=T)))
diff(as.numeric(quantile(cycf, c(0.25, 0.75), na.rm=T)))

library(ncdf4)

dim <- ncdim_def("clijsters_fig1a/time", "seconds", time)
vars <- list(ncvar_def("clijsters_fig1a/E2F1", "signal", dim, NA),
             ncvar_def("clijsters_fig1a/CycF", "signal", dim, NA))

ncnew <- nc_create("Clijsters2019.nc", vars, force_v4=T)
ncvar_put(ncnew, vars[[1]], e2f1)
ncvar_put(ncnew, vars[[2]], cycf)
nc_close(ncnew)
