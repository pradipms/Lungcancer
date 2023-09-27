# https://rdrr.io/cran/fmsb/man/radarchart.html
# Data must be given as the data frame, where the first cases show maximum.
maxmin_RMSE <- data.frame(
  DT=c(39.231, 30.341),
  GBT=c(24.830, 23.114),
  RF=c(27.238, 25.061),
  SVR=c(24.676, 24.278),
  COX=c(30.591, 27.507))

maxmin_RRMSE <- data.frame(
  DT=c(1.367, 1.057),
  GBT=c(0.865, 0.805),
  RF=c(0.949, 0.873),
  SVR=c(0.860, 0.846),
  COX=c(1.066, 0.959))

maxmin_MSE <- data.frame(
  DT=c(1539.108, 920.561),
  GBT=c(616.523, 534.250),
  RF=c(741.885, 628.058),
  SVR=c(608.923, 595.835),
  COX=c(935.804, 756.646))

maxmin_MAE <- data.frame(
  DT=c(25.645, 22.694),
  GBT=c(17.391, 17.120),
  RF=c(21.199, 18.183),
  SVR=c(16.532, 16.230),
  COX=c(25.717, 23.325))

# data for radarchart function version 1 series, minimum value must be omitted from above.
RNGkind("Mersenne-Twister")
set.seed(123)
RMSE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(31.520, 39.231, 30.341),
  GBT = c(24.146, 23.114, 24.830),
  RF = c(25.061, 27.199, 27.238),
  SVR = c(24.278, 24.676, 24.410),
  COX = c(27.507, 30.591, 28.162)
)
RMSE <- rbind(maxmin_RMSE,RMSE)


RRMSE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(1.098, 1.367, 1.057),
  GBT = c(0.841, 0.805, 0.865),
  RF = c(0.873, 0.948, 0.949),
  SVR = c(0.846, 0.860, 0.851),
  COX = c(0.959, 1.066, 0.981)
)
RRMSE <- rbind(maxmin_RRMSE,RRMSE)

MSE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(993.485, 1539.108, 920.561),
  GBT = c(583.030, 534.250, 616.523),
  RF = c(628.058, 739.771, 741.885),
  SVR = c(589.414, 608.923, 595.835),
  COX = c(756.646, 935.804, 793.122)
)
MSE <- rbind(maxmin_MSE,MSE)

MAE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(22.774, 25.645, 22.694),
  GBT = c(17.391, 17.333, 17.120),
  RF = c(18.183, 21.199, 19.397),
  SVR = c(16.244, 16.532, 16.230),
  COX = c(23.325, 25.717, 23.746)
)
MAE <- rbind(maxmin_MAE,MAE)

# Prepare color
colors_border=c( "#00AFBB", "#E7B800", "#FC4E07"  )

# Custom the radarChart For RMSE
radarchart(RMSE, axistype=2, pcol=colors_border, plty=1, axislabcol="grey", pdensity=c(5, 10, 30), 
           pangle=c(10, 45, 120), pfcol=colors_border, 
)

# Legend
legend(x=0.85, y=1, legend = c("GE", "MF", "MF+GE"), bty = "n", pch=20 ,  col=colors_border,  text.col = "black", cex=0.9, pt.cex=1.6)

# Custom the radarChart For RRMSE
radarchart(RRMSE, axistype=2, pcol=colors_border, plty=1, axislabcol="grey", pdensity=c(5, 10, 30), 
           pangle=c(10, 45, 120), pfcol=colors_border, 
)

# Legend
legend(x=0.85, y=1, legend = c("GE", "MF", "MF+GE"), bty = "n", pch=20 ,  col=colors_border,  text.col = "black", cex=0.9, pt.cex=1.6)

# Custom the radarChart For MSE
radarchart(MSE, axistype=2, pcol=colors_border, plty=1, axislabcol="grey", pdensity=c(5, 10, 30), 
           pangle=c(10, 45, 120), pfcol=colors_border, 
)

# Legend
legend(x=0.85, y=1, legend = c("GE", "MF", "MF+GE"), bty = "n", pch=20 ,  col=colors_border,  text.col = "black", cex=0.9, pt.cex=1.6)

# Custom the radarChart For MAE
radarchart(MAE, axistype=2, pcol=colors_border, plty=1, axislabcol="grey", pdensity=c(5, 10, 30), 
           pangle=c(10, 45, 120), pfcol=colors_border, 
)

# Legend
legend(x=0.85, y=1, legend = c("GE", "MF", "MF+GE"), bty = "n", pch=20 ,  col=colors_border,  text.col = "black", cex=0.9, pt.cex=1.6)

