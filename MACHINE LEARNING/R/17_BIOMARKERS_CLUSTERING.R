# Load Libraries
library("FactoMineR")
library("factoextra")
library(cluster)
library(NbClust)
library("readxl")

# Importing Dataset
biomarkers <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/CLUSTERING/ABSOLUTE_FEATURES_BIOMARKERS.csv")

X <- biomarkers[2:18]

res.pca <- PCA(X, graph = FALSE)

eig.val <- get_eigenvalue(res.pca)

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50), ggtheme = theme_light())


fviz_eig( res.pca, choice = c("variance", "eigenvalue"), geom = c("bar", "line"), barfill = "steelblue", barcolor = "black", linecolor = "red", ncp = 10, addlabels = TRUE, hjust = 0, main = NULL, xlab = NULL, ylab = NULL,
         ggtheme = theme_bw()
)


var <- get_pca_var(res.pca)


# Coordinates of variables
head(var$coord, 17)    



# Plot Variables
fviz_pca_var(res.pca, col.var = "black", ggtheme = theme_bw())

# Access cos2 values
head(var$cos2, 4)

# Visualize cos2 variables

library("corrplot")
corrplot(var$cos2, is.corr=FALSE)

# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2, ggtheme = theme_bw())

# Color by cos2 values: quality on the factor map############
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
             ggtheme = theme_bw()
             )

# Change the transparency by cos2 values
fviz_pca_var(res.pca, alpha.var = "cos2", ggtheme = theme_bw())

# Contributions of variables
head(var$contrib, 17)

# Correlation of contributions
library("corrplot")
corrplot(var$contrib, is.corr=FALSE)  

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10, ggtheme = theme_bw())
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10, ggtheme = theme_bw())

# Total contributions of PC1 and PC2
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10, ggtheme = theme_bw())

# Plot Contributing variables

fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             ggtheme = theme_bw()
)

# Change the transparency by contrib values
fviz_pca_var(res.pca, alpha.var = "contrib", ggtheme = theme_bw())


# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(res.pca, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF", "#52595D", "#837E7C", "#A0CFEC", "#16E2F5"),
             legend.title = "Cluster", ggtheme = theme_bw())

fviz_pca_ind(res.pca, 
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = biomarkers$OverallSurvivalMonths, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#52595D", "#837E7C", "#A0CFEC", "#16E2F5"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Survival Months (Groups)",
             ggtheme = theme_bw()
)

fviz_pca_biplot(res.pca, 
                col.ind = biomarkers$OverallSurvivalMonths, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Overall Survival Months", ggtheme = theme_bw()) 

fviz_pca_biplot(res.pca, 
                # Individuals
                geom.ind = "point",
                fill.ind = biomarkers$OverallSurvivalMonths, col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                addEllipses = TRUE,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdYlBu", ggtheme = theme_bw(),
                
                legend.title = list(fill = "Species", color = "Contrib",
                                    alpha = "Contrib")
)


### Method 2: Silhouette in Base R 

avg_sil <- function(k){
  km.res <- kmeans(X.scaled, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(res.pca))
  mean(ss[,3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15 

### Method 2: Silhouette in FactoExtra
fviz_nbclust(X, kmeans, method = "silhouette") + theme_bw()


# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)
res.km <- kmeans(var$coord, centers = 2, nstart = 25)
grp <- as.factor(res.km$cluster)


# Color variables by groups
fviz_pca_var(res.pca, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF", "#52595D", "#837E7C", "#A0CFEC", "#16E2F5"),
             legend.title = "Cluster", ggtheme = theme_bw())

# Color by cos2 values: quality on the factor map############
fviz_pca_var(res.pca, col.var = "cos2",
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF", "#52595D", "#837E7C", "#A0CFEC", "#16E2F5"),
             legend.title = "Cluster", 
             ggtheme = theme_bw()
)
