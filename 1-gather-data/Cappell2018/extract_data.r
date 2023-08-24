library(ncdf4)
source("../extract_data_functions.r")

fig1c <- extract_data_from_SVG("Figure1c.svg", seq(-3, 10, by=1), seq(1, 0, by=-0.2), point_layer_name="median_x5F_apc_x5F_activity")

plot(fig1c$x, fig1c$y)


# fig2d_emi1 <- extract_data_from_SVG("Figure2d.svg", seq(0, 12, by=1), seq(1.5,0, by=-0.25), point_layer_name="emi1")
# fig2d_geminin <- extract_data_from_SVG("Figure2d.svg", seq(0, 12, by=1), seq(1.5,0, by=-0.25), point_layer_name="geminin")
# fig2d_cycline <- extract_data_from_SVG("Figure2d.svg", seq(0, 12, by=1), seq(1.5,0, by=-0.25), point_layer_name="cycline")
# fig2d_edu <- extract_data_from_SVG("Figure2d.svg", seq(0, 12, by=1), seq(1.5,0, by=-0.25), point_layer_name="edu")
# 
# plot(fig2d_emi1$x, fig2d_emi1$y, pch=19, col='green')
# points(fig2d_geminin$x, fig2d_geminin$y, pch=19, col='red')
# points(fig2d_cycline$x, fig2d_cycline$y, pch=19, col='blue')
# points(fig2d_edu$x, fig2d_edu$y, pch=19, col='black')
# 
# library(ncdf4)
# dims <- list(ncdim_def("Cappell_fig1c/time", "seconds", fig1c$x * 3600),
#              ncdim_def("Cappell_fig2d/time", "seconds", fig2d_emi1$x * 3600))
# vars <- list(ncvar_def("Cappell_fig1c/APC_activity", "activity", dims[[1]], NA),
#              ncvar_def("Cappell_fig2d/Emi1", "normalized_protein_abundance", dims[[2]], NA),
#              ncvar_def("Cappell_fig2d/Geminin", "normalized_protein_abundance", dims[[2]], NA),
#              ncvar_def("Cappell_fig2d/CyclinE", "normalized_protein_abundance", dims[[2]], NA),
#              ncvar_def("Cappell_fig2d/edu", "signal", dims[[2]], NA),
#              ncvar_def("Cappell_fig2d/edu_plus_one", "signal", dims[[2]], NA))
# ncnew <- nc_create("Cappell2018.nc", vars, force_v4=T)
# ncvar_put(ncnew, vars[[1]], fig1c$y)
# ncvar_put(ncnew, vars[[2]], fig2d_emi1$y)
# ncvar_put(ncnew, vars[[3]], fig2d_geminin$y)
# ncvar_put(ncnew, vars[[4]], fig2d_cycline$y)
# ncvar_put(ncnew, vars[[5]], fig2d_edu$y)
# ncvar_put(ncnew, vars[[6]], fig2d_edu$y+1)
# nc_close(ncnew)

df2c <- read.delim("fig2c_data.txt", stringsAsFactors = F)
df2d <- read.delim("fig2d_data.txt", stringsAsFactors = F)

plot(df2c$Time.since.mitosis..hrs., df2c$E2F1, type='l')
lines(df2c$Time.since.mitosis..hrs., df2c$E2F1.1, type='l')

diff(as.numeric(quantile(c(df2d$CycE, df2c$CycE.1), c(0.25, 0.75), na.rm=T)))
diff(as.numeric(quantile(c(df2d$Emi1, df2d$Emi1.1), c(0.25, 0.75), na.rm=T)))
diff(as.numeric(quantile(c(df2d$GMNN, df2d$GMNN.1), c(0.25, 0.75), na.rm=T)))
diff(as.numeric(quantile(c(df2d$EdU, df2d$EdU.1), c(0.25, 0.75), na.rm=T)))



# dims <- list(ncdim_def("Cappell_fig1c/time", "seconds", fig1c$x * 3600),
#              ncdim_def("Cappell_fig2d/time", "seconds", df$TIme.since.mitosis..hrs. * 3600),
#              ncdim_def("Cappell_fig2d/cells", "id", 1:2))
# vars <- list(ncvar_def("Cappell_fig1c/APC_activity", "activity", dims[1], NA),
#              ncvar_def("Cappell_fig2d/Emi1_rep1", "normalized_protein_abundance", dims[3:2], NA),
#              ncvar_def("Cappell_fig2d/Geminin", "normalized_protein_abundance", dims[3:2], NA),
#              ncvar_def("Cappell_fig2d/CyclinE", "normalized_protein_abundance", dims[3:2], NA),
#              ncvar_def("Cappell_fig2d/edu", "signal", dims[3:2], NA),
#              ncvar_def("Cappell_fig2d/edu_plus_one", "signal", dims[3:2], NA))
# ncnew <- nc_create("Cappell2018.nc", vars, force_v4=T)
# ncvar_put(ncnew, vars[[1]], fig1c$y)
# ncvar_put(ncnew, vars[[2]], df$Emi1,   start=c(1,1), count=c(1, nrow(df)))
# ncvar_put(ncnew, vars[[2]], df$Emi1.1, start=c(2,1), count=c(1, nrow(df)))
# ncvar_put(ncnew, vars[[3]], df$GMNN,   start=c(1,1), count=c(1, nrow(df)))
# ncvar_put(ncnew, vars[[3]], df$GMNN.1, start=c(2,1), count=c(1, nrow(df)))
# ncvar_put(ncnew, vars[[4]], df$CycE,   start=c(1,1), count=c(1, nrow(df)))
# ncvar_put(ncnew, vars[[4]], df$CycE.1, start=c(2,1), count=c(1, nrow(df)))
# ncvar_put(ncnew, vars[[5]], df$EdU,    start=c(1,1), count=c(1, nrow(df)))
# ncvar_put(ncnew, vars[[5]], df$EdU.1,  start=c(2,1), count=c(1, nrow(df)))
# ncvar_put(ncnew, vars[[6]], df$EdU+1,   start=c(1,1), count=c(1, nrow(df)))
# ncvar_put(ncnew, vars[[6]], df$EdU.1+1, start=c(2,1), count=c(1, nrow(df)))
# nc_close(ncnew)

dims <- list(ncdim_def("Cappell_fig1c/time", "seconds", fig1c$x * 3600),
             ncdim_def("Cappell_fig2/time", "seconds", df$TIme.since.mitosis..hrs. * 3600))
vars <- list(ncvar_def("Cappell_fig1c/APC_activity", "activity", dims[1], NA),
             ncvar_def("Cappell_fig2/Emi1_mRNA_rep1", "normalized_transcript_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/Emi1_mRNA_rep2", "normalized_transcript_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/Geminin_mRNA_rep1", "normalized_transcript_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/Geminin_mRNA_rep2", "normalized_transcript_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/CyclinE_mRNA_rep1", "normalized_transcript_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/CyclinE_mRNA_rep2", "normalized_transcript_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/E2F1_mRNA_rep1", "normalized_transcript_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/E2F1_mRNA_rep2", "normalized_transcript_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/Emi1_protein_rep1", "normalized_protein_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/Emi1_protein_rep2", "normalized_protein_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/Geminin_protein_rep1", "normalized_protein_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/Geminin_protein_rep2", "normalized_protein_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/CyclinE_protein_rep1", "normalized_protein_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/CyclinE_protein_rep2", "normalized_protein_abundance", dims[2], NA),
             ncvar_def("Cappell_fig2/edu_rep1", "signal", dims[2], NA),
             ncvar_def("Cappell_fig2/edu_rep2", "signal", dims[2], NA),
             ncvar_def("Cappell_fig2/edu_plus_one_rep1", "signal", dims[2], NA),
             ncvar_def("Cappell_fig2/edu_plus_one_rep2", "signal", dims[2], NA))
ncnew <- nc_create("Cappell2018.nc", vars, force_v4=T)
ncvar_put(ncnew, vars[[ 1]], fig1c$y)
ncvar_put(ncnew, vars[[ 2]], df2c$Emi1,   start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[ 3]], df2c$Emi1.1, start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[ 4]], df2c$GMNN,   start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[ 5]], df2c$GMNN.1, start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[ 6]], df2c$CycE,   start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[ 7]], df2c$CycE.1, start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[ 8]], df2c$E2F1,   start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[ 9]], df2c$E2F1.1, start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[10]], df2d$Emi1,   start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[11]], df2d$Emi1.1, start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[12]], df2d$GMNN,   start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[13]], df2d$GMNN.1, start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[14]], df2d$CycE,   start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[15]], df2d$CycE.1, start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[16]], df2d$EdU,    start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[17]], df2d$EdU.1,  start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[18]], df2d$EdU+1,  start=c(1), count=c(nrow(df)))
ncvar_put(ncnew, vars[[19]], df2d$EdU.1+1,start=c(1), count=c(nrow(df)))
nc_close(ncnew)