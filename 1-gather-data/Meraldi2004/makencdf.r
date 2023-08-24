# Manually read off from Figure 1B
durations <- c(18,21,21,21,21,21,24,24,24,27,27,27,30,30,42,45,45,60,75,81) * 60

library(ncdf4)
dim1 <- ncdim_def("meraldi_fig1b/cells", "id", 1:length(durations))
var1 <- ncvar_def("meraldi_fig1b/NEBD_to_AO_duration", "seconds", list(dim1), NA)
ncnew <- nc_create("Meraldi2004.nc", list(var1), force_v4=T)
ncvar_put(ncnew, var1, durations, start=1, count=length(durations))
nc_close(ncnew)

diff(as.numeric(quantile(durations, c(0.25, 0.75))))
