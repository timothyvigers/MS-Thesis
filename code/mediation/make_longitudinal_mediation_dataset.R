pheno = read.csv("/home/biostats_share/Norris/data/phenotype/ivyomicssample_noIdentifyingInfo.csv",stringsAsFactors = F,na.strings = "")
pheno = pheno[pheno$Visit_Type == "SV" | pheno$Visit_Type == "PSV",]
# Code from candidate selection
# Methylation
load("/home/biostats_share/Norris/data/methylation/Mmatrix.platformAdj.regressOut.Rdata")
names = colnames(M.adj)
methyl = as.data.frame(t(M.adj))
rownames(methyl) = names
# Keys
key_450k = read.csv("/home/biostats_share/Norris/data/methylation/key.450K.csv",stringsAsFactors = F)
key_450k$platform = "450K"
key_epic = read.csv("/home/biostats_share/Norris/data/methylation/key.EPIC.csv",stringsAsFactors = F)
key_epic$platform = "EPIC"
key = rbind(key_450k,key_epic)
# Make final methylation dataset
methyl$samplekey = key$samplekey[match(rownames(methyl),key$array)]
# Import metabolites
gctof = read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.bc.csv",stringsAsFactors = F)
gctof = gctof[gctof$samplekey %in% pheno$samplekey,]
hilic = read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.bc.csv",stringsAsFactors = F)
hilic = hilic[hilic$samplekey %in% pheno$samplekey,]
lipid = read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.bc.csv",stringsAsFactors = F)
lipid = lipid[lipid$samplekey %in% pheno$samplekey,]
oxylipin = read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.bc.csv",stringsAsFactors = F)
oxylipin = oxylipin[oxylipin$samplekey %in% pheno$samplekey,]
vitd = read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.bc.csv",stringsAsFactors = F)
vitd = vitd[vitd$samplekey %in% pheno$samplekey,]
# Merge
df = merge(pheno,methyl,by = "samplekey",all.x = T)
df = merge(df,gctof,by = "samplekey",all.x = T)
df = merge(df,hilic,by = "samplekey",all.x = T)
df = merge(df,lipid,by = "samplekey",all.x = T)
df = merge(df,oxylipin,by = "samplekey",all.x = T)
df = merge(df,vitd,by = "samplekey",all.x = T)
df = df[!duplicated(df$samplekey),]
# Save
all_data = data.frame(df)
save(all_data,file = "~/MS-Thesis/data/networks/longitudinal_mediation.Rdata")
