setwd("/home/vigerst/EWAS")
# Import original data
load("./data/Mmatrix.platformAdj.regressOut.Rdata")
load("./data/probesFromPipeline.Rdata")
pheno = read.csv("./data/Infant_Diet_Methylation_Phenotype_Final_19MARCH2021.csv")
# Subset methylation data
M.adj = M.adj[,colnames(M.adj) %in% pheno$Array]
# Merge
M.adj = t(M.adj) 
df = cbind(pheno,M.adj[match(pheno$Array,rownames(M.adj)),])
# Write merged data
save(df,file = "./data/final_data.RData")
