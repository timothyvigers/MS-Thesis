# Data import
# Phenotype
pheno <- read.csv("/home/biostats_share/Norris/data/phenotype/ivyomicssample_noIdentifyingInfo.csv",
                  stringsAsFactors = F,
                  na.strings = "")
pheno <- pheno[with(pheno,order(ID)),]
pheno <- pheno[!is.na(pheno$T1Dgroup),]
pheno <- pheno[pheno$Visit_Type == "SV",]
# Probes
load("/home/biostats_share/Norris/data/methylation/probesFromPipeline.Rdata")
# Methylation
load("/home/biostats_share/Norris/data/methylation/Mmatrix.platformAdj.Rdata")
methyl <- as.data.frame(t(M.adj))
# Keys
key_450k <- read.csv("/home/biostats_share/Norris/data/methylation/key.450K.csv",
                     stringsAsFactors = F)
key_450k$platform <- "450K"
key_epic <- read.csv("/home/biostats_share/Norris/data/methylation/key.EPIC.csv",
                     stringsAsFactors = F)
key_epic$platform <- "EPIC"
key <- rbind(key_450k,key_epic)
# Make final methylation dataset
methyl$samplekey <- key$samplekey[match(rownames(methyl),key$array)]
methyl <- methyl[match(pheno$samplekey,methyl$samplekey),]
methyl$id <- factor(pheno$ID[match(methyl$samplekey,pheno$samplekey)])
methyl <- methyl[,c(probesFromPipeline,"samplekey","id")]
# Add sex and age
methyl[,1:(ncol(methyl)-2)] <- lapply(methyl[,1:(ncol(methyl)-2)],function(x){
  x + pheno$clinage[match(methyl$samplekey,pheno$samplekey)] + 
    as.numeric(factor(pheno$SEX[match(methyl$samplekey,pheno$samplekey)]))
})
# Write CSV
out_dir = "/home/vigerst/MS-Thesis/data/step_1/sv/"
write.csv(methyl,file = paste0(out_dir,"methyl_adj.csv"),row.names = F)
# Import metabolites, scale, add sex and age
gctof <- read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.bc.csv",
                  stringsAsFactors = F)
gctof <- gctof[gctof$samplekey %in% pheno$samplekey,]
gctof[,2:ncol(gctof)] <- lapply(gctof[,2:ncol(gctof)],scale)
gctof <- lapply(gctof[,2:ncol(gctof)],function(x){
  x + pheno$clinage[match(gctof$samplekey,pheno$samplekey)] + 
    as.numeric(factor(pheno$SEX[match(gctof$samplekey,pheno$samplekey)]))
})
# Write CSV
write.csv(gctof,file = paste0(out_dir,"gctof_adj.csv"),row.names = F)

hilic <- read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.bc.csv",
                  stringsAsFactors = F)
hilic <- hilic[hilic$samplekey %in% pheno$samplekey,]
hilic[,2:ncol(hilic)] <- lapply(hilic[,2:ncol(hilic)],scale)
hilic <- lapply(hilic[,2:ncol(hilic)],function(x){
  x + pheno$clinage[match(hilic$samplekey,pheno$samplekey)] + 
    as.numeric(factor(pheno$SEX[match(hilic$samplekey,pheno$samplekey)]))
})
# Write CSV
write.csv(hilic,file = paste0(out_dir,"hilic_adj.csv"),row.names = F)

lipid <- read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.bc.csv",
                  stringsAsFactors = F)
lipid <- lipid[lipid$samplekey %in% pheno$samplekey,]
lipid[,2:ncol(lipid)] <- lapply(lipid[,2:ncol(lipid)],scale)
lipid <- lapply(lipid[,2:ncol(lipid)],function(x){
  x + pheno$clinage[match(lipid$samplekey,pheno$samplekey)] + 
    as.numeric(factor(pheno$SEX[match(lipid$samplekey,pheno$samplekey)]))
})
# Write CSV
write.csv(lipid,file = paste0(out_dir,"lipid_adj.csv"),row.names = F)

oxylipin <- read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.bc.csv",
                     stringsAsFactors = F)
oxylipin <- oxylipin[oxylipin$samplekey %in% pheno$samplekey,]
oxylipin[,2:ncol(oxylipin)] <- lapply(oxylipin[,2:ncol(oxylipin)],scale)
oxylipin <- lapply(oxylipin[,2:ncol(oxylipin)],function(x){
  x + pheno$clinage[match(oxylipin$samplekey,pheno$samplekey)] + 
    as.numeric(factor(pheno$SEX[match(oxylipin$samplekey,pheno$samplekey)]))
})
# Write CSV
write.csv(oxylipin,file = paste0(out_dir,"oxylipin_adj.csv"),row.names = F)

vitd <- read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.bc.csv",
                 stringsAsFactors = F)
vitd <- vitd[vitd$samplekey %in% pheno$samplekey,]
vitd[,2:ncol(vitd)] <- lapply(vitd[,2:ncol(vitd)],scale)
vitd <- lapply(vitd[,2:ncol(vitd)],function(x){
  x + pheno$clinage[match(vitd$samplekey,pheno$samplekey)] + 
    as.numeric(factor(pheno$SEX[match(vitd$samplekey,pheno$samplekey)]))
})
# Write CSV
write.csv(vitd,file = paste0(out_dir,"vitd_adj.csv"),row.names = F)
