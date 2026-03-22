# Install dualsimplex package
library(devtools)
unloadNamespace("dualsimplex")
devtools::install_github("artyomovlab/DualSimplex@v1.0")
#devtools::load_all("../../dualsimplex") to install from local directory
library(DualSimplex)
options("dualsimplex-rasterize" = T)
