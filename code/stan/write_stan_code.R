library(bayesvl)
# Working directory
setwd("C:/Users/timbv/Documents/GitHub/MS-Thesis")
# Create DAGs for all network structures
# Nodes
base <- bayesvl()
base <- bvl_addNode(base, "methyl", "norm")
base <- bvl_addNode(base, "metab", "norm")
base <- bvl_addNode(base, "t1d", "bern")
# Add arcs
# 1
s1 <- bvl_addArc(base,"metab","methyl")
s1 <- bvl_addArc(s1,"t1d","metab")
s1 <- bvl_addArc(s1,"t1d","methyl")
# 2
s2 <- bvl_addArc(base,"methyl","metab")
s2 <- bvl_addArc(s2,"t1d","metab")
s2 <- bvl_addArc(s2,"t1d","methyl")
# 3
s3 <- bvl_addArc(base,"methyl","metab")
s3 <- bvl_addArc(s3,"metab","t1d")
s3 <- bvl_addArc(s3,"methyl","t1d")
# 4
s4 <- bvl_addArc(base,"methyl","metab")
s4 <- bvl_addArc(s4,"methyl","t1d")
s4 <- bvl_addArc(s4,"t1d","metab")
# 5
s5 <- bvl_addArc(base,"metab","methyl")
s5 <- bvl_addArc(s5,"methyl","t1d")
s5 <- bvl_addArc(s5,"metab","t1d")
# 6
s6 <- bvl_addArc(base,"metab","methyl")
s6 <- bvl_addArc(s6,"t1d","methyl")
s6 <- bvl_addArc(s6,"metab","t1d")
# 7
s7 <- bvl_addArc(base,"t1d","methyl")
s7 <- bvl_addArc(s7,"t1d","metab")
# 8
s8 <- bvl_addArc(base,"methyl","t1d")
s8 <- bvl_addArc(s8,"t1d","metab")
# 9
s9 <- bvl_addArc(base,"t1d","methyl")
s9 <- bvl_addArc(s9,"metab","t1d")
# 10
s10 <- bvl_addArc(base,"methyl","t1d")
s10 <- bvl_addArc(s10,"metab","t1d")
# 11
s11 <- bvl_addArc(base,"methyl","t1d")
s11 <- bvl_addArc(s11,"methyl","metab")
# 12
s12 <- bvl_addArc(base,"t1d","methyl")
s12 <- bvl_addArc(s12,"methyl","metab")
# 13
s13 <- bvl_addArc(base,"methyl","t1d")
s13 <- bvl_addArc(s13,"metab","methyl")
# 14
s14 <- bvl_addArc(base,"t1d","methyl")
s14 <- bvl_addArc(s14,"metab","methyl")
# 15
s15 <- bvl_addArc(base,"metab","t1d")
s15 <- bvl_addArc(s15,"metab","methyl")
# 16
s16 <- bvl_addArc(base,"t1d","metab")
s16 <- bvl_addArc(s16,"metab","methyl")
# 17
s17 <- bvl_addArc(base,"metab","t1d")
s17 <- bvl_addArc(s17,"methyl","metab")
# 18
s18 <- bvl_addArc(base,"t1d","metab")
s18 <- bvl_addArc(s18,"methyl","metab")
# 19
s19 <- bvl_addArc(base,"t1d","methyl")
# 20
s20 <- bvl_addArc(base,"t1d","metab")
# 21
s21 <- bvl_addArc(base,"methyl","t1d")
# 22
s22 <- bvl_addArc(base,"methyl","metab")
# 23
s23 <- bvl_addArc(base,"metab","t1d")
# 24
s24 <- bvl_addArc(base,"metab","methyl")
# Write all stan code, just in case
lapply(paste0("s",1:24), function(x){
  dag <- get(x)
  write(bvl_model2Stan(dag),file = paste0("./code/stan/",x,".stan"))
  })
