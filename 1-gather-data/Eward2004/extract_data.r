source("../extract_data_functions.r")

cycEexp1 <- extract_data_from_SVG("Figure2_exp1cycE.svg", seq(0, 42, len=15), seq(8, 0, by=-2), point_layer_name = "CyclinE")
cycAexp1 <- extract_data_from_SVG("Figure2_exp1cycA.svg", seq(0, 45, len=16), seq(8, 0, by=-2), point_layer_name = "CyclinA")
cycBexp1 <- extract_data_from_SVG("Figure2_exp1cycB.svg", seq(0, 45, len=16), seq(8, 0, by=-2), point_layer_name = "CyclinB")
cycEexp2 <- extract_data_from_SVG("Figure2_exp2cycE.svg", seq(0, 45, len=16), seq(8, 0, by=-2), point_layer_name = "CyclinE")
cycAexp2 <- extract_data_from_SVG("Figure2_exp2cycA.svg", seq(0, 45, len=16), seq(8, 0, by=-2), point_layer_name = "CyclinA")
cycBexp2 <- extract_data_from_SVG("Figure2_exp2cycB.svg", seq(0, 45, len=16), seq(8, 0, by=-2), point_layer_name = "CyclinB")

min(cycAexp1$y)
min(cycEexp1$y)
min(cycBexp1$y)
min(cycAexp2$y)
min(cycEexp2$y)
min(cycBexp2$y)

df <- data.frame(time = seq(0, 45, by=3) * 3600)
df$cycA1 <- (cycAexp1$y - 1) / 9
df$cycE1 <- NA
df$cycE1[1:15] <- (cycEexp1$y - 1) / 3
df$cycB1 <- (cycBexp1$y - 1) / 9
df$cycA2 <- (cycAexp2$y - 1) / 9
df$cycE2 <- (cycEexp2$y - 1) / 3
df$cycB2 <- (cycBexp2$y - 1) / 9

mean(diff(as.numeric(quantile(df$cycA1, c(0.25, 0.75), na.rm=T))),
     diff(as.numeric(quantile(df$cycA2, c(0.25, 0.75), na.rm=T))))
mean(diff(as.numeric(quantile(df$cycE1, c(0.25, 0.75), na.rm=T))),
     diff(as.numeric(quantile(df$cycE2, c(0.25, 0.75), na.rm=T))))
mean(diff(as.numeric(quantile(df$cycB1, c(0.25, 0.75), na.rm=T))),
     diff(as.numeric(quantile(df$cycB2, c(0.25, 0.75), na.rm=T))))

# Use only the first two cell cycles; cells might have diverged after that,
# and simulating large populations take too long for not too much additional information
# Also exclude the first timepoint - cells are still exiting the cell cycle, which would be
# fine if we simulate the preceding cell cycle, but then we would again have another cell cycle
# which is too costly
use_ix <- 2:12
plot(df$time[use_ix]/3600, df$cycA1[use_ix])
plot(df$time[use_ix]/3600, df$cycE1[use_ix])
plot(df$time[use_ix]/3600, df$cycB1[use_ix])

library(ncdf4)
dim <- ncdim_def("Eward_Fig2/time", "seconds", df$time[use_ix])
vars <- list(ncvar_def("Eward_Fig2/cycAexp1", "signal", dim, NA),
             ncvar_def("Eward_Fig2/cycAexp2", "signal", dim, NA),
             ncvar_def("Eward_Fig2/cycEexp1", "signal", dim, NA),
             ncvar_def("Eward_Fig2/cycEexp2", "signal", dim, NA),
             ncvar_def("Eward_Fig2/cycBexp1", "signal", dim, NA),
             ncvar_def("Eward_Fig2/cycBexp2", "signal", dim, NA))

ncnew <- nc_create("Eward2004.nc", vars, force_v4=T)
ncvar_put(ncnew, vars[[1]], df$cycA1[use_ix])
ncvar_put(ncnew, vars[[2]], df$cycA2[use_ix])
ncvar_put(ncnew, vars[[3]], df$cycE1[use_ix])
ncvar_put(ncnew, vars[[4]], df$cycE2[use_ix])
ncvar_put(ncnew, vars[[5]], df$cycB1[use_ix])
ncvar_put(ncnew, vars[[6]], df$cycB2[use_ix])
nc_close(ncnew)


