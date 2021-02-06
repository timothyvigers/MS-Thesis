setwd("C:/Users/Tim Vigers/Documents/GitHub/MS-Thesis/data/raw_data")
# Data import
# Phenotype
pheno = read.csv("./ivyomicssample_noIdentifyingInfo.csv",stringsAsFactors = F,na.strings = "")
# Probes
load("./probesFromPipeline.Rdata")
# Methylation
load("./Mmatrix.platformAdj.regressOut.Rdata")
methyl = as.data.frame(t(M.adj))
# Keys
key_450k = read.csv("./key.450K.csv",
                     stringsAsFactors = F)
key_450k$platform = "450K"
key_epic = read.csv("./key.EPIC.csv",
                     stringsAsFactors = F)
key_epic$platform = "EPIC"
key = rbind(key_450k,key_epic)
# Make final methylation dataset
methyl$samplekey = key$samplekey[match(rownames(methyl),key$array)]
rm("key","key_450k","key_epic")
# Metabolites
metab_candidates = read.csv("./liz_candidates.csv",stringsAsFactors = F,na.strings = "")
# gctof
gctof = read.csv("./gctof.bc.csv")
gctof = gctof[,c("samplekey",unique(metab_candidates$gctof))]
hilic = read.csv("./hilic.bc.csv")
hilic = hilic[,c("samplekey","hilic_12")]
lipid = read.csv("./lipid.bc.csv")
lipid = lipid[,c("samplekey",unique(metab_candidates$lipid))]
oxylipin = read.csv("./oxylipin.bc.csv")
oxylipin = oxylipin[!duplicated(oxylipin$samplekey),]
vitd = read.csv("./vitD.bc.csv")
# Make big dataframe
df = merge(pheno,methyl,by = "samplekey")
df = merge(df,gctof,by = "samplekey")
df = merge(df,hilic,by = "samplekey")
df = merge(df,lipid,by = "samplekey")
df = merge(df,oxylipin,by = "samplekey")
df = merge(df,vitd,by = "samplekey")
# Pre-SV and SV only
df = df[df$Visit_Type %in% c("PSV","SV"),]
df = df[!duplicated(df$samplekey),]
rownames(df) = df$samplekey
# Save
save(df,file = "./psv_sv_dataset.Rdata")
