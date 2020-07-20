library(rjags)
library(boot)
library(parallel)
# Load data
setwd("/home/vigerst/MS-Thesis")
load("./data/networks/pair_data.Rdata")
load("./data/networks/cits.Rdata")