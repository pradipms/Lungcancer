# https://rdrr.io/cran/fmsb/man/radarchart.html
# Data must be given as the data frame, where the first cases show maximum.
maxmin_RMSE <- data.frame(
  DT=c(48.404, 42.520),
  GBT=c(41.584, 40.186),
  RF=c(40.112, 37.545),
  SVR=c(38.013, 37.924),
  COX=c(39.369, 37.786))

maxmin_RRMSE <- data.frame(
  DT=c(1.510, 1.327),
  GBT=c(1.298, 1.254),
  RF=c(1.252, 1.171),
  SVR=c(1.186, 1.183),
  COX=c(1.228, 1.179))

maxmin_MSE <- data.frame(
  DT=c(2343.028, 1807.969),
  GBT=c(1729.301, 1614.962),
  RF=c(1609.025, 1409.637),
  SVR=c(1445.017, 1438.256),
  COX=c(1549.968, 1427.815))

maxmin_MAE <- data.frame(
  DT=c(28.739, 25.923),
  GBT=c(23.770, 22.224),
  RF=c(22.905, 21.720),
  SVR=c(19.533, 19.339),
  COX=c(27.120, 25.968))

# data for radarchart function version 1 series, minimum value must be omitted from above.
RNGkind("Mersenne-Twister")
set.seed(123)
RMSE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(48.404, 42.520, 47.717),
  GBT = c(41.584, 40.186, 41.365),
  RF = c(40.112, 37.545, 39.976),
  SVR = c(37.992, 38.013, 37.924),
  COX = c(38.312, 39.369, 37.786)
)
RMSE <- rbind(maxmin_RMSE,RMSE)


RRMSE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(1.510, 1.327, 1.489),
  GBT = c(1.298, 1.254, 1.291),
  RF = c(1.252, 1.171, 1.247),
  SVR = c(1.185, 1.186, 1.183),
  COX = c(1.195, 1.228, 1.179)
)
RRMSE <- rbind(maxmin_RRMSE,RRMSE)

MSE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(2343.028, 1807.969, 2276.972),
  GBT = c(1729.301, 1614.962, 1711.124),
  RF = c(1609.025, 1409.637, 1598.140),
  SVR = c(1443.439, 1445.017, 1438.256),
  COX = c(1467.861, 1549.968, 1427.815)
)
MSE <- rbind(maxmin_MSE,MSE)

MAE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(28.562, 25.923, 28.739),
  GBT = c(22.224, 23.770, 22.558),
  RF = c(22.905, 22.894, 21.720),
  SVR = c(19.351, 19.533, 19.339),
  COX = c(26.638, 27.120, 25.968)
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

