# Load Libraries
library('devtools')
library("FactoMineR")
library("factoextra")
library("cluster")
library("NbClust")
library("readxl")
library("corrplot")
library(ggpubr)  # Exporting plot in pdf

# Importing Dataset
common_features <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/ANALYTICS2/FLUX_COMMON_FEATURES.csv")

X <- common_features[2:6]

X.scaled <- scale(X) # Scale the data

# Compute PCA
a_res.pca <- PCA(X)
res.pca <- prcomp(X)

eig.val <- get_eigenvalue(a_res.pca)
eig.val

eig.val1 <- get_eigenvalue(res.pca)

#Plot PCA graph of individuals colored according to cos2 values for each dataset
abc_PCA_plot <- fviz_pca_ind(a_res.pca, col.ind = "cos2",
                                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                     repel = TRUE # Avoid text overlapping
)

PCA_plot <- fviz_pca_ind(res.pca, col.ind = "cos2",
                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                             repel = TRUE # Avoid text overlapping
)

contributions_a <-a_res.pca$var$contrib



eig.val <- get_eigenvalue(res.pca)

# Correlation Plot
common_features_corr_plot <- corrplot(cor(a_res.pca), type = 'upper',  tl.cex = 0.9, tl.col = "black")

# Scree Plot
flux_features_scree_plot <- fviz_eig( res.pca, choice = c("variance", "eigenvalue"), geom = c("bar", "line"), barfill = "steelblue", barcolor = "black", linecolor = "#FC4E07", ncp = 10, addlabels = TRUE, hjust = 0, main = NULL, xlab = NULL, ylab = NULL,
                                      ggtheme = theme_bw()
)


flux_features_scree_plot1 <- fviz_eig( a_res.pca, choice = c("variance", "eigenvalue"), geom = c("bar", "line"), barfill = "steelblue", barcolor = "black", linecolor = "#FC4E07", ncp = 10, addlabels = TRUE, hjust = 0, main = NULL, xlab = NULL, ylab = NULL,
                                      ggtheme = theme_bw()
)

# Define
var <- get_pca_var(a_res.pca)

# Color by cos2 values: quality on the factor map
flux_cos2_biplot <- fviz_pca_var(a_res.pca, col.var = "cos2",
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                 repel = TRUE, # Avoid text overlapping
                                 ggtheme = theme_bw()
)

### Method: Silhouette in Base R 

avg_sil <- function(k){
  km.res <- kmeans(a_res.pca, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(res.pca))
  mean(ss[,3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15 

### Method 2: Silhouette in FactoExtra
flux_silhoutte_plot <- fviz_nbclust(X, kmeans, method = "silhouette") + theme_bw()


# Create a grouping variable using kmeans
# Create 2 groups of variables (centers = 2)
set.seed(123)
res.km <- kmeans(var$coord, centers = 4, nstart = 25)
grp <- as.factor(res.km$cluster)

# Color by cos2 values: quality on the factor map
fviz_pca_var(a_res.pca, col.var = "cos2",
             gradient.cols = c("#0073C2FF", "#EFC000FF", "#868686FF", "#00AFBB"), 
             repel = TRUE # Avoid text overlapping
)

# Color variables by groups
fviz_pca_var(res.pca, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF", "#00AFBB"),
             legend.title = "Cluster")


# Its possible to add the cluster assignments to the original data set 
aggregate(res.pca, by=list(cluster = res.km$cluster), mean)

df2 <- cbind(res.pca, cluster = res.km$cluster)
df2

fviz_cluster(kmeans(res.km, centers = 4), geom = "point")

res.pca
res.km
fviz_pca_var(res.pca)

# Convex hull
fviz_pca_ind(res.pca, geom.ind = "point",
             col.ind = common_features$OverallSurvivalMonths, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#00AFBB", "red", "brown", "white"),
             addEllipses = TRUE, ellipse.type = "convex",
             legend.title = "Groups"
)

# Bi Plot
flux_features_cont_survival_plot <- fviz_pca_biplot(a_res.pca, 
                                                    # Individuals
                                                    geom.ind = "point",
                                                    fill.ind = common_features$OverallSurvivalMonths, col.ind = "black",
                                                    pointshape = 21, pointsize = 2,
                                                    palette = "jco",
                                                    addEllipses = TRUE,
                                                    # Variables
                                                    alpha.var ="contrib", col.var = "contrib",
                                                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), ggtheme = theme_bw(),
                                                    repel = TRUE, # Avoid text overlapping (slow if many points)
                                                    
                                                    legend.title = list(fill = "Overall Survival Months \n(Groups)", color = "Contrib",
                                                                        alpha = "Contrib")
)

fviz_pca_biplot(res.pca, 
                col.ind = common_features$OverallSurvivalMonths, palette = c("#0073C2FF", "#EFC000FF", "#868686FF", "#00AFBB", "red", "brown", "white"), 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Species")

var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}


loadings <- a_res.pca$rotation
sdev <- a_res.pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
head(var$coord[, 1:4])
var$cos2
var$contrib

abc < -a_res.pca$var$coord
