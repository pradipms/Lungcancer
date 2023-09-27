#Load libraries

library(fmsb)

# Construct the data set

# PCA data

PCA_radar_data <- data.frame("Alanine, aspartate and glutamate metabolism"=c(0.6,-0.2,0.300225603,-0.242814741),
                             "Amino sugar and nucleotide sugar metabolism"=c(0.6,-0.2,0.056394634,0.047659606),
                             "Butanoate metabolism"=c(0.6,-0.2,0.063702704,0.154918447),
                             "Cholesterol metabolism"=c(0.6,-0.2,0.380759574,0.125140988),
                             "D-Glutamine and D-glutamate metabolism"=c(0.6,-0.2,0.300225603,-0.242814741),
                             "Fatty acid biosynthesis"=c(0.6,-0.2,0.398997390,0.207062550),
                             "Biosynthesis of unsaturated fatty acids"=c(0.6,-0.2,0.398997390,0.207062550),
                             "Glycine, serine and threonine metabolism"=c(0.6,-0.2,-0.064908974,0.375766748),
                             "Glyoxylate and dicarboxylate metabolism"=c(0.6,-0.2,0.050380016,0.129077620),
                             "Pentose and glucuronate interconversions"=c(0.6,-0.2,0.298531061,0.046461876),
                             "Pentose phosphate pathway"=c(0.6,-0.2,0.298531061,0.046461876),
                             "Purine metabolism"=c(0.6,-0.2,0.073124748,0.091739794),
                             "Sphingolipid metabolism"=c(0.6,-0.2,0.288110239,0.085537016),
                             "Tryptophan metabolism"=c(0.6,-0.2,0.214409171,0.046493792),
                             "Valine, leucine and isoleucine degradation"=c(0.6,-0.2,-0.122057685,0.535876496),
                             "Valine, leucine and isoleucine biosynthesis"=c(0.6,-0.2,-0.122057685,0.535876496),
                             row.names = c("max", "min", "Component 1", "Component 2"))

# RF data

RF_radar_data <- data.frame("Alanine, aspartate and glutamate metabolism"=c(0.4,-0.4,0.198221683,-0.06668472),
                            "Amino sugar and nucleotide sugar metabolism"=c(0.4,-0.4,0.079220392,-0.035180158),
                            "Arginine and proline metabolism"=c(0.4,-0.4,0.144362933,0.118201562),
                            "Cholesterol metabolism"=c(0.4,-0.4,-0.032532042,-0.305366568),
                            "Cysteine and methionine metabolism"=c(0.4,-0.4,0.146932277,-0.096137363),
                            "D-Glutamine and D-glutamate metabolism"=c(0.4,-0.4,0.135569477,0.309533436),
                            "Drug metabolism - cytochrome P450"=c(0.4,-0.4,0.244880396,-0.018040488),
                            "Drug metabolism - other enzymes"=c(0.4,-0.4,0.244880396,-0.018040488),
                            "Fatty acid biosynthesis"=c(0.4,-0.4,0.189981091,-0.16636335),
                            "Biosynthesis of unsaturated fatty acids"=c(0.4,-0.4,0.189981091,-0.16636335),
                            "Fructose and mannose metabolism"=c(0.4,-0.4,-0.100943901,0.193819349),
                            "Galactose metabolism"=c(0.4,-0.4,-0.140160674,0.123220936),
                            "Glutathione metabolism"=c(0.4,-0.4,0.223086942,0.056615116),
                            "Glycerophospholipid metabolism"=c(0.4,-0.4,0.157720574,-0.017758468),
                            "Glycine, serine and threonine metabolism"=c(0.4,-0.4,0.164858085,-0.054687317),
                            "Glycolysis / Gluconeogenesis"=c(0.4,-0.4,0.172178893,0.149013364),
                            "Glyoxylate and dicarboxylate metabolism"=c(0.4,-0.4,0.14822735,0.003017949),
                            "Nitrogen metabolism"=c(0.4,-0.4,0.116150886,0.029320625),
                            "Oxidative phosphorylation"=c(0.4,-0.4,0.038756367,-0.064177596),
                            "Pentose and glucuronate interconversions"=c(0.4,-0.4,0.144134583,0.027254276),
                            "Pentose phosphate pathway"=c(0.4,-0.4,0.185392225,0.06751943),
                            "Primary bile acid biosynthesis"=c(0.4,-0.4,0.236897701,0.153222951),
                            "Secondary bile acid biosynthesis"=c(0.4,-0.4,0.236897701,0.153222951),
                            "Purine metabolism"=c(0.4,-0.4,0.252864883,-0.050934927),
                            "Pyrimidine metabolism"=c(0.4,-0.4,0.197450752,0.008117793),
                            "Pyruvate metabolism"=c(0.4,-0.4,0.171520262,-0.037825028),
                            "Riboflavin metabolism"=c(0.4,-0.4,0.087273463,0.187732482),
                            "Starch and sucrose metabolism"=c(0.4,-0.4,-0.128222707,0.125609144),
                            "Taurine and hypotaurine metabolism"=c(0.4,-0.4,0.209041949,0.13118677),
                            "Folate biosynthesis"=c(0.4,-0.4,0.247987776,-0.114551472),
                            "Valine, leucine and isoleucine degradation"=c(0.4,-0.4,0.003991474,-0.463676682),
                            "Valine, leucine and isoleucine biosynthesis"=c(0.4,-0.4,0.003991474,-0.463676682),
                            "Fatty acid degradation"=c(0.4,-0.4,0.217718504,0.076258173),
                            "Type II diabetes mellitus"=c(0.4,-0.4,0.160499501,-0.257053177),
                            
                            row.names = c("max", "min", "Component 1", "Component 2"))



# Define line colors

colors_line_PCA <- c(scales::alpha("green3", 0.9),
                     
                     scales::alpha("orangered", 0.9))

colors_line_RF <- c(scales::alpha("darkblue", 0.9),
                    
                    scales::alpha("red3", 0.9))

# Create plot

radarchart(PCA_radar_data, 
           
           seg = 11,  # Number of axis segments
           
           title = "Average Component Score Radar Chart",
           
           pcol = colors_line_PCA,
           
           plty = 1:1,
           
           plwd = 2,
           
           axistype = 4,
           
           caxislabels = c("-0.2","-0.08","-0.06","-0.04","-0.02","0","0.1","0.2","0.3","0.4","0.5","0.6"),
           
           cglty = 3,
           
           cglcol = "gray96",
           
           axislabcol="gray0")

# Add a legend

legend(x=1.35, y=1.25, legend = rownames(PCA_radar_data[-c(1,2),]), bty = "o", pch=20 , col=colors_line_PCA , text.col = "gray0", cex=1.2, pt.cex=3)


radarchart(RF_radar_data, 
           
           seg = 9,  # Number of axis segments
           
           title = "Average Component Score Radar Chart",
           
           pcol = colors_line_RF,
           
           plty = 1:1,
           
           plwd = 2,
           
           axistype = 4,
           
           caxislabels = c("-0.4","-0.2","-0.09","-0.06","-0.03","0","0.1","0.2","0.3","0.4"),
           
           cglty = 3,
           
           cglcol = "gray96",
           
           axislabcol="gray0")

# Add a legend

legend(x=1.35, y=1.25, legend = rownames(RF_radar_data[-c(1,2),]), bty = "o", pch=20 , col=colors_line_RF , text.col = "gray0", cex=1.2, pt.cex=3)



