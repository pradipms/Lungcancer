install.packages('devtools')
devtools::install_github(repo = "ucsf-ferguson-lab/syndRomics@*release")

install.packages('remotes')
remotes::install_github(repo = "ucsf-ferguson-lab/syndRomics@*release")

install.packages("dplyr")


library("devtools")
library(syndRomics)
library(ggplot2)
library(dplyr)
library(stringr)


# Importing Dataset
sum_pca_pathways <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/PCA/SUM_PCA_LOADING/PCA_SUM_PATHWAYS.csv")
sum_rf_pathways <- read.csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/PCA/SUM_PCA_LOADING/RF_SUM_PATHWAYS.csv")

X <- sum_pca_pathways[2:17]
Y <- sum_rf_pathways[2:35]

# Compute PCA For PCA Pathways
pca_data<-as.matrix(apply(X, 2, as.numeric))
pca<-prcomp(pca_data, center = T, scale. = T)

s_per_pca<-permut_pc_test(pca, pca_data, P=1000, ndim=2, statistic = 's.loadings', perm.method = 'permV')

s_per_pca$results

pc_pca<-s_per_pca$results

# Compute PCA For RF Pathways
y_data <-as.matrix(apply(Y, 2, as.numeric))
y <-prcomp(y_data, center = T, scale. = T)

s_per_rf<-permut_pc_test(y, y_data, P=1000, ndim=2, statistic = 's.loadings', perm.method = 'permV')

s_per_rf$results

pc_rf<-s_per_rf$results


# Plot Barmap For PCA Pathways
barmap_pca_pathways<-barmap_loading(pca, pca_data, ndim=1:2, cutoff = 0.45, star_values = F,text_values = T, gradient_color = TRUE)
# Plot Barmap For RF Pathways
barmap_rf_pathways<-barmap_loading(pca=y, pca_data=y_data, ndim=1:2, cutoff = 0.45, star_values = F,text_values = T)

# Heatmap For PCA Pathways
heatmap_pca_pathways<-heatmap_loading(pca, pca_data, ndim=1:2, cutoff = 0.45, star_values = T, text_values = T,
                                      colors= c("orange2", "white","purple"))
# Heatmap For RF Pathways
heatmap_rf_pathways<-heatmap_loading(pca=y, pca_data=y_data, ndim=1:2, cutoff = 0.45, star_values = T, text_values = T,
                                     colors= c("orange2", "white","purple"))

#Save PCA and RF plots to .pdf files in the workspace
pdf("barmap_pca_pathways.pdf")
print(barmap_pca_pathways)
dev.off()

pdf("barmap_rf_pathways.pdf")
print(barmap_rf_pathways)
dev.off()

pdf("heatmap_pca_pathways.pdf")
print(heatmap_pca_pathways)
dev.off()

pdf("heatmap_rf_pathways.pdf")
print(heatmap_rf_pathways)
dev.off()




