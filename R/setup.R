# Install Linseed 2
#devtools::load_all("../../dualsimplex")
unloadNamespace("dualsimplex")
devtools::install_github("artyomovlab/DualSimplex@knn_gene_filtering")
library("DualSimplex")