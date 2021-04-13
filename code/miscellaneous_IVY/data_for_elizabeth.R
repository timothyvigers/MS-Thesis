# Import
sample_info = read.csv("~/Dropbox/School/Statistical Genomics/Final Project/data/Infant_Diet_Methylation_Phenotype_Final_19MARCH2021.csv")
load("~/Dropbox/School/MS Thesis/data/raw_data/Mmatrix.platformAdj.regressOut.Rdata")
# Filter probes
candidates = read.delim("~/Dropbox/School/Statistical Genomics/EWAS/final_candidates.txt",header = F)
missing = candidates$V1[which(!(candidates$V1 %in% rownames(M.adj)))]
candidates = candidates$V1[which((candidates$V1 %in% rownames(M.adj)))]
M.adj = M.adj[candidates,]
# Transpose and merge
M.adj = as.data.frame(t(M.adj))
M.adj$Array = rownames(M.adj)
df = merge(sample_info,M.adj,by = "Array")
# Write
write.csv(df,file = "~/Dropbox/School/Statistical Genomics/EWAS/data/data_for_elizabeth.csv",row.names = F,na = "")
write.csv(missing,file = "~/Dropbox/School/Statistical Genomics/EWAS/data/missing_probes.csv",row.names = F,na = "")