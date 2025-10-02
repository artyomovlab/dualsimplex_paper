# Install dualsimplex package
library(devtools)
unloadNamespace("dualsimplex")
devtools::install_github("artyomovlab/DualSimplex")
#devtools::load_all("../../dualsimplex") to install from local directory
library(DualSimplex)
options("dualsimplex-rasterize" = T)
