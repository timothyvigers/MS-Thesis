# Combine all necessary data into one DF
pheno = read.csv("/home/biostats_share/Norris/data/phenotype/ivyomicssample_noIdentifyingInfo.csv",stringsAsFactors = F,na.strings = "")
pheno = pheno[!is.na(pheno$T1Dgroup),]
pheno = pheno[pheno$Visit_Type == "SV",]
# Code from candidate selection
# Methylation
load("/home/biostats_share/Norris/data/methylation/Mmatrix.platformAdj.regressOut.Rdata")
methyl = as.data.frame(t(M.adj))
# Scale 
methyl = lapply(methyl, scale)
# Keys
key_450k = read.csv("/home/biostats_share/Norris/data/methylation/key.450K.csv",stringsAsFactors = F)
key_450k$platform = "450K"
key_epic = read.csv("/home/biostats_share/Norris/data/methylation/key.EPIC.csv",stringsAsFactors = F)
key_epic$platform = "EPIC"
key = rbind(key_450k,key_epic)
# Make final methylation dataset
methyl$samplekey = key$samplekey[match(rownames(methyl),key$array)]
# Import metabolites and scale
gctof = read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.bc.csv",stringsAsFactors = F)
gctof = gctof[gctof$samplekey %in% pheno$samplekey,]
#gctof = gctof[,c("samplekey",colnames(gctof)[which(colnames(gctof) %in% unique(pairs$metab))])]
gctof[,2:ncol(gctof)] = lapply(gctof[,2:ncol(gctof)],scale)
hilic = read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.bc.csv",stringsAsFactors = F)
hilic = hilic[hilic$samplekey %in% pheno$samplekey,]
#hilic = hilic[,c("samplekey","hilic_12")]
hilic[,2:ncol(hilic)] = lapply(hilic[,2:ncol(hilic)],scale)
lipid = read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.bc.csv",stringsAsFactors = F)
lipid = lipid[lipid$samplekey %in% pheno$samplekey,]
#lipid = lipid[,c("samplekey",colnames(lipid)[which(colnames(lipid) %in% unique(pairs$metab))])]
lipid[,2:ncol(lipid)] = lapply(lipid[,2:ncol(lipid)],scale)
oxylipin = read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.bc.csv",stringsAsFactors = F)
oxylipin = oxylipin[oxylipin$samplekey %in% pheno$samplekey,]
#oxylipin = oxylipin[,c("samplekey",colnames(oxylipin)[which(colnames(oxylipin) %in% unique(pairs$metab))])]
oxylipin[,2:ncol(oxylipin)] = lapply(oxylipin[,2:ncol(oxylipin)],scale)
vitd = read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.bc.csv",stringsAsFactors = F)
vitd = vitd[vitd$samplekey %in% pheno$samplekey,]
#vitd = vitd[,c("samplekey",colnames(vitd)[which(colnames(vitd) %in% unique(pairs$metab))])]
vitd[,2:ncol(vitd)] = lapply(vitd[,2:ncol(vitd)],scale)
# Merge
df = merge(pheno,methyl,by = "samplekey",all.x = T)
df = merge(df,gctof,by = "samplekey",all.x = T)
df = merge(df,hilic,by = "samplekey",all.x = T)
df = merge(df,lipid,by = "samplekey",all.x = T)
df = merge(df,oxylipin,by = "samplekey",all.x = T)
df = merge(df,vitd,by = "samplekey",all.x = T)
all_data = df
# Save
save(all_data,file = "~/MS-Thesis/data/networks/all_data_methyl_scaled.Rdata")