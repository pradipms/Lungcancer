# https://rdrr.io/cran/fmsb/man/radarchart.html
# Data must be given as the data frame, where the first cases show maximum.
maxmin_RMSE <- data.frame(
  DT=c(45.586, 29.739),
  GBT=c(24.022, 23.232),
  RF=c(26.631, 23.817),
  SVR=c(24.647, 24.266),
  COX=c(35.169, 29.523))

maxmin_RRMSE <- data.frame(
  DT=c(1.589, 1.036),
  GBT=c(0.837, 0.810),
  RF=c(0.928, 0.830),
  SVR=c(0.859, 0.846),
  COX=c(1.226, 1.029))

maxmin_MSE <- data.frame(
  DT=c(2078.075, 884.393),
  GBT=c(577.076, 539.735),
  RF=c(709.221, 567.249),
  SVR=c(607.470, 588.862),
  COX=c(1236.844, 871.629))

maxmin_MAE <- data.frame(
  DT=c(29.994, 21.673),
  GBT=c(17.837, 17.084),
  RF=c(20.285, 17.609),
  SVR=c(16.480, 16.292),
  COX=c(25.136, 21.475))

# data for radarchart function version 1 series, minimum value must be omitted from above.
RNGkind("Mersenne-Twister")
set.seed(123)
RMSE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(29.739, 45.586, 34.401),
  GBT = c(23.610, 24.022, 23.232),
  RF = c(23.817, 26.631, 24.900),
  SVR = c(24.266, 24.647, 24.397),
  COX = c(35.169, 29.747, 29.523)
)
RMSE <- rbind(maxmin_RMSE,RMSE)


RRMSE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(1.036, 1.589, 1.199),
  GBT = c(0.823, 0.837, 0.810),
  RF = c(0.830, 0.928, 0.868),
  SVR = c(0.846, 0.859, 0.850),
  COX = c(1.226, 1.037, 1.029)
)
RRMSE <- rbind(maxmin_RRMSE,RRMSE)

MSE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(884.393, 2078.075, 1183.451),
  GBT = c(557.411, 577.076, 539.735),
  RF = c(567.249, 709.221, 620.031),
  SVR = c(588.862, 607.470, 595.232),
  COX = c(1236.844, 884.908, 871.629)
)
MSE <- rbind(maxmin_MSE,MSE)

MAE <- data.frame(
  row.names = c("GE", "MF", "MF+GE"),
  DT = c(21.673, 29.994, 24.900),
  GBT = c(17.084, 17.837, 17.163),
  RF = c(17.609, 20.285, 17.694),
  SVR = c(16.292, 16.480, 16.359),
  COX = c(23.174, 25.136, 21.475)
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

