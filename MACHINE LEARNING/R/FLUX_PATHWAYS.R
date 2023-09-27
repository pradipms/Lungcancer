# Load Libraries
library("FactoMineR")
library("factoextra")
library("cluster")
library("NbClust")
library("readxl")
library("corrplot")
library(ggpubr)

# Importing Dataset
pathways <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/CLUSTERING/ABSOLUTE_FLUX_PATHWAYS.csv")

X <- pathways[2:5]

X.scaled <- scale(X) # Scale the data

# Correlation Plot
flux_pathways_corr_plot <- corrplot(cor(X.scaled), type = 'upper',  tl.cex = 0.9, tl.col = "black")

flux_pathways_corr_plot


res.pca <- PCA(X, graph = FALSE)

eig.val <- get_eigenvalue(res.pca)

# Scree Plot
flux_pathways_scree_plot <- fviz_eig( res.pca, choice = c("variance", "eigenvalue"), geom = c("bar", "line"), barfill = "steelblue", barcolor = "black", linecolor = "red", ncp = 10, addlabels = TRUE, hjust = 0, main = NULL, xlab = NULL, ylab = NULL,
                                      ggtheme = theme_bw() + theme(
                                                plot.background = element_blank()
                                                  ,panel.grid.major = element_blank()
                                                      ,panel.grid.minor = element_blank()
                                                        ) +  theme(axis.line = element_line(color = "black"))
)

# Define
var <- get_pca_var(res.pca)


# Coordinates of variables
head(var$coord, 17)    

# Access cos2 values
head(var$cos2, 4)

# Visualize cos2 variables
flux_pathways_cos2_plot <- corrplot(var$cos2, is.corr=FALSE, tl.col = "black")


# Total cos2 of variables on Dim.1 and Dim.2
flux_pathways_cos2_dim_plot <- fviz_cos2(res.pca, choice = "var", axes = 1:2, ggtheme = theme_bw() + theme(
                                                          plot.background = element_blank()
                                                            ,panel.grid.major = element_blank()
                                                                ,panel.grid.minor = element_blank()
                                                                    ) +  theme(axis.line = element_line(color = "black")))

# Color by cos2 values: quality on the factor map
flux_pathways_cos2_biplot <- fviz_pca_var(res.pca, col.var = "cos2",
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                 repel = TRUE, # Avoid text overlapping
                                 ggtheme = theme_bw()+ theme(
                                          plot.background = element_blank()
                                              ,panel.grid.major = element_blank()
                                                  ,panel.grid.minor = element_blank()
                                                      ) +  theme(axis.line = element_line(color = "black"))
)


# Contributions of variables
head(var$contrib, 17)

# Correlation of contributions
flux_pathways_cont_corr_plot <- corrplot(var$contrib, is.corr=FALSE, tl.col = "black")  

# Contributions of variables to PC1
flux_pathways_cont_dim1_plot <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10, ggtheme = theme_bw() + theme(
                                                    plot.background = element_blank()
                                                        ,panel.grid.major = element_blank()
                                                            ,panel.grid.minor = element_blank()
                                                                ) +  theme(axis.line = element_line(color = "black")))
# Contributions of variables to PC2
flux_pathways_cont_dim2_plot <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10, ggtheme = theme_bw()+ theme(
                                                    plot.background = element_blank()
                                                        ,panel.grid.major = element_blank()
                                                            ,panel.grid.minor = element_blank()
                                                                ) +  theme(axis.line = element_line(color = "black")))

# Total contributions of PC1 and PC2
flux_pathways_cont_dim1_dim2_plot <- fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10, ggtheme = theme_bw()+ theme(
                                                    plot.background = element_blank()
                                                        ,panel.grid.major = element_blank()
                                                            ,panel.grid.minor = element_blank()
                                                                ) +  theme(axis.line = element_line(color = "black")))

# Plot Contributing variables
flux_pathways_cont_biplot <- fviz_pca_var(res.pca, col.var = "contrib",
                                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                          repel = TRUE, # Avoid text overlapping
                                          ggtheme = theme_bw()+ theme(
                                                    plot.background = element_blank()
                                                        ,panel.grid.major = element_blank()
                                                            ,panel.grid.minor = element_blank()
                                                                ) +  theme(axis.line = element_line(color = "black")),
)


### Method: Silhouette in Base R 

avg_sil <- function(k){
  km.res <- kmeans(X.scaled, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(res.pca))
  mean(ss[,3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15 

### Method 2: Silhouette in FactoExtra
flux_pathways_silhoutte_plot <- fviz_nbclust(X, kmeans, method = "silhouette", linecolor = "orange") + theme_bw()+ theme(
                                                                            plot.background = element_blank()
                                                                              ,panel.grid.major = element_blank()
                                                                                ,panel.grid.minor = element_blank()
                                                                                  ) +  theme(axis.line = element_line(color = "black"))


# Create a grouping variable using kmeans
# Create 2 groups of variables (centers = 2)
set.seed(123)
res.km <- kmeans(var$coord, centers = 2, nstart = 25)
grp <- as.factor(res.km$cluster)

# Bi Plot
flux_pathways_cont_survival_plot <- fviz_pca_biplot(res.pca, 
                                                    # Individuals
                                                    geom.ind = "point",
                                                    fill.ind = pathways$OverallSurvivalMonths, col.ind = "black",
                                                    pointshape = 21, pointsize = 3,
                                                    palette = "ucscgb",
                                                    addEllipses = TRUE,
                                                    # Variables
                                                    alpha.var ="contrib", col.var = "contrib",
                                                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), ggtheme = theme_bw()+ theme(
                                                                                      plot.background = element_blank()
                                                                                        ,panel.grid.major = element_blank()
                                                                                            ,panel.grid.minor = element_blank()
                                                                                                ) +  theme(axis.line = element_line(color = "black")),
                                                    repel = TRUE, # Avoid text overlapping (slow if many points)
                                                    
                                                    legend.title = list(fill = "Overall Survival Months \n(Groups)", color = "Contrib",
                                                                        alpha = "Contrib")
)

# Export plot in Pdf
ggexport(plotlist = list(flux_pathways_corr_plot, flux_pathways_scree_plot, flux_pathways_cos2_plot, flux_pathways_cos2_dim_plot, flux_pathways_cos2_biplot,
                         flux_pathways_cont_corr_plot, flux_pathways_cont_dim1_plot, flux_pathways_cont_dim2_plot, flux_pathways_cont_dim1_dim2_plot,
                         flux_pathways_cont_biplot, flux_pathways_silhoutte_plot, flux_pathways_cont_survival_plot), 
         filename = "Flux_Pathways.pdf")

# Export into a TXT file
write.infile(res.pca, "ABSOLUTE_FLUX_PATHWAYS.txt", sep = "\t")

# Export into a CSV file
write.infile(res.pca, "ABSOLUTE_FLUX_PATHWAYS", sep = ";")
