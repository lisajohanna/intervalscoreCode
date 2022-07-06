# =============================================================================
# Load required packages
# =============================================================================
library(ggplot2)  # for function ggprop()
library(scales)   # for function ggprop()
library(ggpubr)   # for function ggarrange() in Figure 1

# =============================================================================
# Set up parallel computing
# (functions to export need to be defined beforehand)
# =============================================================================
library(foreach)    # for parallel computing
library(doParallel) # for parallel computing
nbcores <- detectCores()
cl <- makeCluster(nbcores)
registerDoParallel(cl)
# export of a few functions needed due to function exporting issues on Windows
# (usually, fct.s/var.s exported to parallel execution workers automatically)
clusterExport(cl, varlist = c("score", "width", "coverage", "refmeanscore",
                              "wscore", "refmeanwscore"))

# # to stop parallel computing
# stopCluster(cl)