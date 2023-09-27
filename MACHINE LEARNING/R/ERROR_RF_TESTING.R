# https://rdrr.io/cran/fmsb/man/radarchart.html
# Data must be given as the data frame, where the first cases show maximum.
maxmin_RMSE <- data.frame(
  DT=c(46.237, 45.851),
  GBT=c(36.591, 35.010),
  RF=c(38.839, 36.899),
  SVR=c(37.990, 37.758),
  COX=c(39.377, 36.287))

maxmin_RRMSE <- data.frame(
  DT=c(1.443, 1.431),
  GBT=c(1.142, 1.092),
  RF=c(1.212, 1.151),
  SVR=c(1.185, 1.178),
  COX=c(1.229, 1.132))

maxmin_MSE <- data.frame(
  DT=c(2137.889, 2102.363),
  GBT=c(1338.928, 1225.745),
  RF=c(1508.487, 1361.565),
  SVR=c(1443.253, 1425.671),
  COX=c(1550.594, 1316.796))

maxmin_MAE <- data.frame(
  DT=c(29.513, 28.484),
  GBT=c(20.965, 20.231),
  RF=c(23.814, 21.564),
  SVR=c(19.499, 19.369),
  COX=c(26.755, 22.613))

# data for radarchart function version 1 series, minimum value must be omitted from above.
RNGkind("Mersenne-Twister")
set.seed(123)
RMSE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(46.090, 45.851, 46.237),
  GBT = c(35.010, 36.591, 36.326),
  RF = c(38.512, 38.839, 36.899),
  SVR = c(37.758, 37.990, 37.878),
  COX = c(36.393, 39.377, 36.287)
)
RMSE <- rbind(maxmin_RMSE,RMSE)


RRMSE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(1.438, 1.431, 1.443),
  GBT = c(1.092, 1.142, 1.133),
  RF = c(1.202, 1.212, 1.151),
  SVR = c(1.178, 1.185, 1.182),
  COX = c(1.135, 1.229, 1.132)
)
RRMSE <- rbind(maxmin_RRMSE,RRMSE)

MSE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(2124.318, 2102.363, 2137.889),
  GBT = c(1225.745, 1338.928, 1319.631),
  RF = c(1483.187, 1508.487, 1361.565),
  SVR = c(1425.671, 1443.253, 1434.746),
  COX = c(1324.470, 1550.594, 1316.796)
)
MSE <- rbind(maxmin_MSE,MSE)

MAE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(28.484, 29.331, 29.513),
  GBT = c(20.652, 20.965, 20.231),
  RF = c(23.814, 23.269, 21.564),
  SVR = c(19.383, 19.499, 19.369),
  COX = c(22.872, 26.755, 22.613)
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

