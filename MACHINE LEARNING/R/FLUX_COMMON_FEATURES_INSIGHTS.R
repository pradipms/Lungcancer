# Load Libraries
library("devtools")
library("FactoMineR")
library("factoextra")
library("cluster")
library("NbClust")
library("readxl")
library("corrplot")
library("ggpubr")  # Exporting plot in pdf


# Importing Dataset
common_features <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/ANALYTICS2/FLUX_COMMON_FEATURES.csv")

X <- common_features[2:6]

X.scaled <- scale(X) # Scale the data

# Compute PCA
res.pca <- PCA(X, scale = TRUE)

# Compute Eigen Values
eig.val <- get_eigenvalue(res.pca)

# Correlation Plot
common_features_corr_plot <- corrplot(cor(X.scaled), type = 'upper',  tl.cex = 0.9, tl.col = "black")


# Scree Plot
common_features_scree_plot <- fviz_eig( res.pca, choice = c("variance", "eigenvalue"), geom = c("bar", "line"), barfill = "steelblue", barcolor = "black", linecolor = "#FC4E07", ncp = 10, addlabels = TRUE, hjust = 0, main = NULL, xlab = NULL, ylab = NULL,
                                      ggtheme = theme_bw()+ theme(
                                              plot.background = element_blank()
                                                ,panel.grid.major = element_blank()
                                                    ,panel.grid.minor = element_blank()
                                                        ) +  theme(axis.line = element_line(color = "black"))
)

# Define Variable
var <- get_pca_var(res.pca)

# Define Individual
ind <- get_pca_ind(res.pca)

# Color by cos2 values: quality on the factor map
common_features_cos2_biplot <- fviz_pca_var(res.pca, col.var = "cos2",
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                 repel = TRUE, # Avoid text overlapping
                                 ggtheme = theme_bw()+ theme(
                                            plot.background = element_blank()
                                              ,panel.grid.major = element_blank()
                                                  ,panel.grid.minor = element_blank()
                                                      ) +  theme(axis.line = element_line(color = "black"))
)

### Method: Silhouette in Base R 

avg_sil <- function(k){
  km.res <- kmeans(res.pca, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(res.pca))
  mean(ss[,3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15 

### Method 2: Silhouette in FactoExtra
common_features_silhoutte_plot <- fviz_nbclust(X, kmeans, method = "silhouette", linecolor = "MediumSeaGreen") + theme_bw()+ theme(
                                                                                            plot.background = element_blank()
                                                                                              ,panel.grid.major = element_blank()
                                                                                                ,panel.grid.minor = element_blank()
                                                                                                    ) +  theme(axis.line = element_line(color = "black"))


# Create a grouping variable using kmeans
# Create 2 groups of variables (centers = 4)
set.seed(123)
res.km <- kmeans(var$coord, centers = 4, nstart = 25)
grp <- as.factor(res.km$cluster)

# Bi Plot
common_features_cont_plot <- fviz_pca_biplot(res.pca, 
                                                    # Individuals
                                                    geom.ind = "point",
                                                    fill.ind = common_features$OverallSurvivalMonths, col.ind = "black",
                                                    pointshape = 21, pointsize = 3,
                                                    palette = c("Aquamarine", "DeepPink", "DarkOrange", "Yellow", "DarkGreen", "Coral", "Purple"),
                                                    addEllipses = TRUE,
                                                    # Variables
                                                    alpha.var ="contrib", 
                                                    col.var = "contrib",
                                                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                                    ggtheme = theme_bw()+ theme(
                                                              plot.background = element_blank()
                                                                ,panel.grid.major = element_blank()
                                                                    ,panel.grid.minor = element_blank()
                                                                        ) +  theme(axis.line = element_line(color = "black")),
                                                    repel = TRUE, # Avoid text overlapping (slow if many points)
                                                    
                                                    legend.title = list(fill = "Overall Survival Months \n(Groups)", color = "contrib",
                                                                        alpha = "contrib")
)

#Save PCA plots to .pdf files in the workspace
pdf("common_features_scree_plot.pdf")
print(common_features_scree_plot)
dev.off()

pdf("common_features_silhoutte_plot.pdf")
print(common_features_silhoutte_plot)
dev.off()

pdf("common_features_cos2_biplot.pdf")
print(common_features_cos2_biplot)
dev.off()

pdf("common_features_cont_plot.pdf")
print(common_features_cont_plot)
dev.off()


# Export into a TXT file
write.infile(res.pca, "common_features_pcavalues.txt", sep = "\t")

# Export into a CSV file
write.infile(res.pca, "common_features_pcavalues.csv", sep = ";")

