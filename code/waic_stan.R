library(bayesvl)
library(rstan)
# Working directory
setwd("C:/Users/timbv/Documents/GitHub/MS-Thesis")
# Create DAGs for all network structures
# Nodes
base <- bayesvl()
base <- bvl_addNode(base, "methyl", "norm")
base <- bvl_addNode(base, "metab", "norm")
base <- bvl_addNode(base, "t1d", "bern")
# Arcs
# 1
s1 <- bvl_addArc(base,"metab","methyl")
s1 <- bvl_addArc(s1,"t1d","metab")
s1 <- bvl_addArc(s1,"t1d","methyl")
# 2
s2 <- bvl_addArc(base,"methyl","metab")
s2 <- bvl_addArc(s2,"t1d","metab")
s2 <- bvl_addArc(s2,"t1d","methyl")
