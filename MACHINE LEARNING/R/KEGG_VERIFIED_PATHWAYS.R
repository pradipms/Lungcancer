library("ggradar")
library(tidyverse)
library(ggpubr)

#Data
reactions <- data.frame(
  row.names = c("PCA", "RF"),
  'Alanine, aspartate and glutamate metabolism' = c(1, 4),
  'Amino sugar and nucleotide sugar metabolism' = c(3, 16),
  'Butanoate metabolism' = c(1, 0),
  'Cholesterol metabolism' = c(1, 4),
  'D-Glutamine and D-glutamate metabolism' = c(1, 2),
  'Fatty acid biosynthesis' = c(2, 15),
  'Biosynthesis of unsaturated fatty acids' = c(2, 15),
  'Glycine, serine and threonine metabolism' = c(1, 1),
  'Glyoxylate and dicarboxylate metabolism' = c(2, 12),
  'Pentose and glucuronate interconversions' = c(1, 1),
  'Pentose phosphate pathway' = c(1, 4),
  'Purine metabolism' = c(2, 15),
  'Sphingolipid metabolism' = c(1, 0),
  'Tryptophan metabolism' = c(2, 0),
  'Valine, leucine and isoleucine degradation' = c(1, 4),
  'Valine, leucine and isoleucine biosynthesis' = c(1, 4),
  'Arginine and proline metabolism' = c(0, 5),
  'Cysteine and methionine metabolism' = c(0, 2),
  'Drug metabolism - cytochrome P450' = c(0, 9),
  'Drug metabolism - other enzymes' = c(0, 0),
  'Fructose and mannose metabolism' = c(0, 1),
  'Galactose metabolism' = c(0, 1),
  'Glutathione metabolism' = c(0, 2),
  'Glycerophospholipid metabolism' = c(0, 1),
  'Glycolysis / Gluconeogenesis' = c(0, 10),
  'Nitrogen metabolism' = c(0, 1),
  'Oxidative phosphorylation' = c(0, 2),
  'Primary bile acid biosynthesis' = c(0, 6),
  'Secondary bile acid biosynthesis' = c(0, 6),
  'Pyrimidine metabolism' = c(0, 15),
  'Pyruvate metabolism' = c(0, 2),
  'Riboflavin metabolism' = c(0, 1),
  'Starch and sucrose metabolism' = c(0, 2),
  'Taurine and hypotaurine metabolism' = c(0, 1),
  'Folate biosynthesis' = c(0, 10),
  'Fatty acid degradation' = c(0, 2),
  'Type II diabetes mellitus' = c(0, 2)
)

reactions

df2 <- t(reactions) %>%
  as.data.frame() %>%
  rownames_to_column("KEGG verified pathways")
df2

#Plot creation:
ggdotchart(
  df2, x = "KEGG verified pathways", y = "RF",
  add = "segments", 
  ylab = "reactions", title = "Feature"
)

#Data Preparation
df3 <- df2 %>%
  select("KEGG verified pathways", "PCA", "RF") %>%
  pivot_longer(
    cols = c(PCA, RF),
    names_to = "FS",
    values_to = "Number of reactions"
  )

# Plot Creation:
ggdotchart(
  df3, x = "KEGG verified pathways", y = "Number of reactions", 
  group = "FS", color = "FS", palette = "jco",
  add = "segment", position = position_dodge(0.3),
  
)


#Data Preparation
df4 <- df2 %>%
  select("KEGG verified pathways", "PCA", "RF") %>%
  pivot_longer(
    cols = c(PCA, RF),
    names_to = "FS",
    values_to = "Number of reactions"
  )
head(df4)

# Plot Creation
ggdotchart(
  df4, x = "KEGG verified pathways", y = "Number of reactions", 
  group = "FS", color = "FS", palette = "jco",
  add = "segment", position = position_dodge(0.3),
  sorting = "none", facet.by ="FS",
  rotate = TRUE, legend = "none"
)

