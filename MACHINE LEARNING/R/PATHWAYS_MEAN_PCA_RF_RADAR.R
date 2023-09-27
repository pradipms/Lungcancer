library("ggradar")
library(fmsb)

#Data
pathways <- data.frame(
  row.names = c("PCA", "RF"),
  "Alanine, aspartate and glutamate metabolism" = c(495.627322, 220.9887908),
  "D-Glutamine and D-glutamate metabolism" = c(495.627322, 283.9593117),
  "Amino sugar and nucleotide sugar metabolism" = c(399.1262805, 204.6476194),
  "Butanoate metabolism" = c(196.653766, 0),
  "Cholesterol metabolism" = c(222.4233362, 134.2034188),
  "Fatty acid biosynthesis" = c(228.8449206, 149.6526607),
  "Biosynthesis of unsaturated fatty acids" = c(228.8449206, 149.6526607),
  "Glycine, serine and threonine metabolism" = c(74.71636696, 139.332702),
  "Glyoxylate and dicarboxylate metabolism" = c(218.8799096, 357.7955423),
  "Pentose and glucuronate interconversions" = c(9.499763078, 16.21349345),
  "Pentose phosphate pathway" = c(9.499763078, 259.7715163),
  "Purine metabolism" = c(637.6666941, 252.0593252),
  "Sphingolipid metabolism" = c(9.916192492, 0),
  "Tryptophan metabolism" = c(167.0732263, 0),
  "Valine, leucine and isoleucine degradation" = c(149.8265372, 321.0434734),
  "Valine, leucine and isoleucine biosynthesis" = c(149.8265372, 321.0434734),
  "Arginine and proline metabolism" = c(0, 78.34780762),
  "Cysteine and methionine metabolism" = c(0, 25.46814156),
  "Drug metabolism - cytochrome P450" = c(0, 53.01663119),
  "Drug metabolism - other enzymes" = c(0, 53.01663119),
  "Fructose and mannose metabolism" = c(0, 55.35964531),
  "Galactose metabolism" = c(0, 540.5243729),
  "Glutathione metabolism" = c(0, 240.1678502),
  "Glycerophospholipid metabolism" = c(0, 67.31333151),
  "Glycolysis / Gluconeogenesis" = c(0, 506.7897826),
  "Nitrogen metabolism" = c(0, 317.4047401),
  "Oxidative phosphorylation" = c(0, 835.2437407),
  "Primary bile acid biosynthesis" = c(0, 98.60348611),
  "Secondary bile acid biosynthesis" = c(0, 98.60348611),
  "Pyrimidine metabolism" = c(0, 244.4200904),
  "Pyruvate metabolism" = c(0, 208.0930183),
  "Riboflavin metabolism" = c(0, 24.17345561),
  "Starch and sucrose metabolism" = c(0, 82.53449229),
  "Fatty acid degradation" = c(0, 38.35518281),
  "Folate biosynthesis" = c(0, 161.0840324),
  "Type II diabetes mellitus" = c(0, 279.9895122),
  "Taurine and hypotaurine metabolism" = c(0, 13.08694069)
)

# Define the variable ranges: maximum and minimum
max_min <- data.frame(
  "Alanine, aspartate and glutamate metabolism" = c(850, 0),
  "D-Glutamine and D-glutamate metabolism" = c(850, 0),
  "Amino sugar and nucleotide sugar metabolism" = c(850, 0),
  "Butanoate metabolism" = c(850, 0),
  "Cholesterol metabolism" = c(850, 0),
  "Fatty acid biosynthesis" = c(850, 0),
  "Biosynthesis of unsaturated fatty acids" = c(850, 0),
  "Glycine, serine and threonine metabolism" = c(850, 0),
  "Glyoxylate and dicarboxylate metabolism" = c(850, 0),
  "Pentose and glucuronate interconversions" = c(850, 0),
  "Pentose phosphate pathway" = c(850, 0),
  "Purine metabolism" = c(850, 0),
  "Sphingolipid metabolism" = c(850, 0),
  "Tryptophan metabolism" = c(850, 0),
  "Valine, leucine and isoleucine degradation" = c(850, 0),
  "Valine, leucine and isoleucine biosynthesis" = c(850, 0),
  "Arginine and proline metabolism" = c(850, 0),
  "Cysteine and methionine metabolism" = c(850, 0),
  "Drug metabolism - cytochrome P450" = c(850, 0),
  "Drug metabolism - other enzymes" = c(850, 0),
  "Fructose and mannose metabolism" = c(850, 0),
  "Galactose metabolism" = c(850, 0),
  "Glutathione metabolism" = c(850, 0),
  "Glycerophospholipid metabolism" = c(850, 0),
  "Glycolysis / Gluconeogenesis" = c(850, 0),
  "Nitrogen metabolism" = c(850, 0),
  "Oxidative phosphorylation" = c(850, 0),
  "Primary bile acid biosynthesis" = c(850, 0),
  "Secondary bile acid biosynthesis" = c(850, 0),
  "Pyrimidine metabolism" = c(850, 0),
  "Pyruvate metabolism" = c(850, 0),
  "Riboflavin metabolism" = c(850, 0),
  "Starch and sucrose metabolism" = c(850, 0),
  "Fatty acid degradation" = c(850, 0),
  "Folate biosynthesis" = c(850, 0),
  "Type II diabetes mellitus" = c(850, 0),
  "Taurine and hypotaurine metabolism" = c(850, 0)
)
rownames(max_min) <- c("Max", "Min")

# Bind the variable ranges to the data
df <- rbind(max_min, pathways)
df

# Plot the data for PCA
library(fmsb)
pca_data <- df[c("Max", "Min", "PCA"), ]
radarchart(pca_data)

create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

# Reduce plot margin using par()
op <- par(mar = c(1, 2, 2, 1))
create_beautiful_radarchart(pca_data, caxislabels = c(0, 80, 160, 325, 650))
par(op)


# Plot the data for RF
library(fmsb)
rf_data <- df[c("Max", "Min", "RF"), ]
radarchart(rf_data)

create_beautiful_radarchart <- function(data, color = "#FC4E07", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

# Reduce plot margin using par()
op <- par(mar = c(1, 2, 2, 1))
create_beautiful_radarchart(rf_data, caxislabels = c(0, 105, 210, 425, 850))
par(op)

