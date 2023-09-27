# Load Libraries
library("FactoMineR")
library("factoextra")
library("readxl")
library("corrplot")
library(ggpubr)

# Importing Dataset
p1 <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/CLUSTERING/INDIVIDUAL_GROUP/INDIVIDUAL_PATHWAYS1.csv")
p2 <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/CLUSTERING/INDIVIDUAL_GROUP/INDIVIDUAL_PATHWAYS2.csv")
p3 <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/CLUSTERING/INDIVIDUAL_GROUP/INDIVIDUAL_PATHWAYS3.csv")
p4 <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/CLUSTERING/INDIVIDUAL_GROUP/INDIVIDUAL_PATHWAYS4.csv")
p5 <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/CLUSTERING/INDIVIDUAL_GROUP/INDIVIDUAL_PATHWAYS5.csv")
p6 <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/CLUSTERING/INDIVIDUAL_GROUP/INDIVIDUAL_PATHWAYS6.csv")
p7 <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/CLUSTERING/INDIVIDUAL_GROUP/INDIVIDUAL_PATHWAYS7.csv")

p1 <- p1[2:5]
p2 <- p2[2:5]
p3 <- p3[2:5]
p4 <- p4[2:5]
p5 <- p5[2:5]
p6 <- p6[2:5]
p7 <- p7[2:5]


# Feature Scaling
p1.scaled <- scale(p1)
p2.scaled <- scale(p2)
p3.scaled <- scale(p3)
p4.scaled <- scale(p4)
p5.scaled <- scale(p5)
p6.scaled <- scale(p6)
p7.scaled <- scale(p7)


# Correlation Plot
p1_corr_plot <- corrplot(cor(p1.scaled), type = 'upper',  tl.cex = 0.9, tl.col = "black")
p2_corr_plot <- corrplot(cor(p2.scaled), type = 'upper',  tl.cex = 0.9, tl.col = "black")
p3_corr_plot <- corrplot(cor(p3.scaled), type = 'upper',  tl.cex = 0.9, tl.col = "black")
p4_corr_plot <- corrplot(cor(p4.scaled), type = 'upper',  tl.cex = 0.9, tl.col = "black")
p5_corr_plot <- corrplot(cor(p5.scaled), type = 'upper',  tl.cex = 0.9, tl.col = "black")
p6_corr_plot <- corrplot(cor(p6.scaled), type = 'upper',  tl.cex = 0.9, tl.col = "black")
p7_corr_plot <- corrplot(cor(p7.scaled), type = 'upper',  tl.cex = 0.9, tl.col = "black")


res.pca1 <- PCA(p1, graph = FALSE)
res.pca2 <- PCA(p2, graph = FALSE)
res.pca3 <- PCA(p3, graph = FALSE)
res.pca4 <- PCA(p4, graph = FALSE)
res.pca5 <- PCA(p5, graph = FALSE)
res.pca6 <- PCA(p6, graph = FALSE)
res.pca7 <- PCA(p7, graph = FALSE)


# Color by cos2 values: quality on the factor map
p1_cos2_biplot <- fviz_pca_var(res.pca1, col.var = "cos2",
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                               repel = TRUE, # Avoid text overlapping
                               ggtheme = theme_bw()+theme(
                                            plot.background = element_blank()
                                              ,panel.grid.major = element_blank()
                                                  ,panel.grid.minor = element_blank()
                                                      ) +  theme(axis.line = element_line(color = "black"))
)

p2_cos2_biplot <- fviz_pca_var(res.pca2, col.var = "cos2",
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                               repel = TRUE, # Avoid text overlapping
                               ggtheme = theme_bw()+theme(
                                              plot.background = element_blank()
                                                ,panel.grid.major = element_blank()
                                                    ,panel.grid.minor = element_blank()
                                                        ) +  theme(axis.line = element_line(color = "black"))
)

p3_cos2_biplot <- fviz_pca_var(res.pca3, col.var = "cos2",
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                               repel = TRUE, # Avoid text overlapping
                               ggtheme = theme_bw()+theme(
                                                plot.background = element_blank()
                                                  ,panel.grid.major = element_blank()
                                                      ,panel.grid.minor = element_blank()
                                                          ) +  theme(axis.line = element_line(color = "black"))
)

p4_cos2_biplot <- fviz_pca_var(res.pca4, col.var = "cos2",
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                               repel = TRUE, # Avoid text overlapping
                               ggtheme = theme_bw()+theme(
                                                  plot.background = element_blank()
                                                    ,panel.grid.major = element_blank()
                                                        ,panel.grid.minor = element_blank()
                                                            ) +  theme(axis.line = element_line(color = "black"))
)

p5_cos2_biplot <- fviz_pca_var(res.pca5, col.var = "cos2",
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                               repel = TRUE, # Avoid text overlapping
                               ggtheme = theme_bw()+theme(
                                                plot.background = element_blank()
                                                  ,panel.grid.major = element_blank()
                                                      ,panel.grid.minor = element_blank()
                                                          ) +  theme(axis.line = element_line(color = "black"))
)

p6_cos2_biplot <- fviz_pca_var(res.pca6, col.var = "cos2",
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                               repel = TRUE, # Avoid text overlapping
                               ggtheme = theme_bw()+theme(
                                            plot.background = element_blank()
                                              ,panel.grid.major = element_blank()
                                                  ,panel.grid.minor = element_blank()
                                                      ) +  theme(axis.line = element_line(color = "black"))
)

p7_cos2_biplot <- fviz_pca_var(res.pca7, col.var = "cos2",
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                               repel = TRUE, # Avoid text overlapping
                               ggtheme = theme_bw()+theme(
                                              plot.background = element_blank()
                                                ,panel.grid.major = element_blank()
                                                    ,panel.grid.minor = element_blank()
                                                        ) +  theme(axis.line = element_line(color = "black"))
)


# Export plot in Pdf
ggexport(plotlist = list(p1_cos2_biplot, p2_cos2_biplot, p3_cos2_biplot, p4_cos2_biplot, p5_cos2_biplot, p6_cos2_biplot, p7_cos2_biplot), 
         filename = "INDIVIDUAL_PATHWAYS_CORRELATION.pdf")

# Export into a TXT file
write.infile(res.pca1, "INDIVIDUAL_PATHWAYS1.txt", sep = "\t")
write.infile(res.pca2, "INDIVIDUAL_PATHWAYS2.txt", sep = "\t")
write.infile(res.pca3, "INDIVIDUAL_PATHWAYS3.txt", sep = "\t")
write.infile(res.pca4, "INDIVIDUAL_PATHWAYS4.txt", sep = "\t")
write.infile(res.pca5, "INDIVIDUAL_PATHWAYS5.txt", sep = "\t")
write.infile(res.pca6, "INDIVIDUAL_PATHWAYS6.txt", sep = "\t")
write.infile(res.pca7, "INDIVIDUAL_PATHWAYS7.txt", sep = "\t")


