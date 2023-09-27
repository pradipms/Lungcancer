# https://rdrr.io/cran/fmsb/man/radarchart.html
# Data must be given as the data frame, where the first cases show maximum.
# Importing library
library(fmsb) # library for radarchart
library(ggpubr) # Exporting the pdf files

# Define 
maxmin_pca <- data.frame(
  "Alanine, aspartate and glutamate metabolis"=c(0.3836, -0.0801),
  "D-Glutamine and D-glutamate metabolism"=c(0.3836, -0.0801),
  "Amino sugar and nucleotide sugar metabolism"=c(0.0529, 0.0251),
  "Butanoate metabolism"=c(-0.0775, -0.1346),
  "Cholesterol metabolism"=c(0.2009, 0.1002),
  "Fatty acid biosynthesis"=c(0.3983, -0.0216),
  "Biosynthesis of unsaturated fatty acids"=c(0.3983, -0.0216),
  "Glycine, serine and threonine metabolism"=c(0.3983, -0.0216),
  "Glyoxylate and dicarboxylate metabolism"=c(0.1970, 0.0439),
  "Pentose and glucuronate interconversions"=c(0.4129, 0.1993),
  "Pentose phosphate pathway"=c(0.4129, 0.1993),
  "Purine metabolism"=c(0.1234, 0.0074),
  "Sphingolipid metabolism"=c(0.2536, 0.2026),
  "Tryptophan metabolism"=c(0.0868, -0.1456),
  "Valine, leucine and isoleucine degradation"=c(0.4868, -0.1617),
  "Valine, leucine and isoleucine biosynthesis"=c(0.4868, -0.1617))

# data for radarchart function version 1 series, minimum value must be omitted from above.
RNGkind("Mersenne-Twister")
set.seed(123)
pca <- data.frame(
  row.names = c("Component 1", "Component 2"),
  "Alanine, aspartate and glutamate metabolis"=c(0.3836, -0.0801),
  "D-Glutamine and D-glutamate metabolism"=c(0.3836, -0.0801),
  "Amino sugar and nucleotide sugar metabolism"=c(0.0251, 0.0529),
  "Butanoate metabolism"=c(-0.0775, -0.1346),
  "Cholesterol metabolism"=c(0.2009, 0.1002),
  "Fatty acid biosynthesis"=c(0.3983, -0.0216),
  "Biosynthesis of unsaturated fatty acids"=c(0.3983, -0.0216),
  "Glycine, serine and threonine metabolism"=c(0.3983, -0.0216),
  "Glyoxylate and dicarboxylate metabolism"=c(0.0439, 0.1970),
  "Pentose and glucuronate interconversions"=c(0.1993, 0.4129),
  "Pentose phosphate pathway"=c(0.1993, 0.4129),
  "Purine metabolism"=c(0.0074, 0.1234),
  "Sphingolipid metabolism"=c(0.2026, 0.2536),
  "Tryptophan metabolism"=c(0.0868, -0.1456),
  "Valine, leucine and isoleucine degradation"=c(-0.1617, 0.4868),
  "Valine, leucine and isoleucine biosynthesis"=c(-0.1617, 0.4868)
)
pca <- rbind(maxmin_pca,pca)


# Prepare color
colors_border=c( "#00AFBB", "#E7B800"  )

# Custom the radarChart For Sum PCA Pathways
radarchart(pca, axistype=2, pcol=colors_border, plty=1, axislabcol="grey", pdensity=c(5, 10, 30), 
           pangle=c(10, 45, 120), pfcol=colors_border, 
)

# Legend
legend(x=0.85, y=1, legend = c("Component 1", "Component 2"), bty = "n", pch=20 ,  col=colors_border,  text.col = "black", cex=0.9, pt.cex=1.6)

# Export plot in Pdf
ggexport(plotlist = list(sum_pca_plot), filename = "sum_pathways_pca_rf.pdf")


