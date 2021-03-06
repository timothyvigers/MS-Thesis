# Import methylation association with T1D
methyl_step2 = read.csv("~/MS-Thesis/data/candidate_selection/step_2/methyl_adj.csv",stringsAsFactors = F)
methyl_step2 = methyl_step2[order(methyl_step2$p.value),]
# Import metabolite associations with T1D
gctof_step3 = read.csv("~/MS-Thesis/data/candidate_selection/step_3/gctof_adj.csv",stringsAsFactors = F)
hilic_step3 = read.csv("~/MS-Thesis/data/candidate_selection/step_3/hilic_adj.csv",stringsAsFactors = F)
lipid_step3 = read.csv("~/MS-Thesis/data/candidate_selection/step_3/lipid_adj.csv",stringsAsFactors = F)
oxylipin_step3 = read.csv("~/MS-Thesis/data/candidate_selection/step_3/oxylipin_adj.csv",stringsAsFactors = F)
vitd_step3 = read.csv("~/MS-Thesis/data/candidate_selection/step_3/vitd_adj.csv",stringsAsFactors = F)
metab_step3 = rbind(gctof_step3,hilic_step3)
metab_step3 = rbind(metab_step3,lipid_step3)
metab_step3 = rbind(metab_step3,oxylipin_step3)
metab_step3 = rbind(metab_step3,vitd_step3)
metab_step3 = metab_step3[order(metab_step3$p.value),]
# Define cutoffs
pheno_cutoff = 0.05
pair_cutoff = 0.001
# Get methyl and metabs associated with phenotype
methyl_pheno = unique(methyl_step2$methyl[methyl_step2$p.value < pheno_cutoff])
metab_pheno = unique(metab_step3$metab[metab_step3$p.value < pheno_cutoff])
# Read in methyl ~ metab results
gctof = read.csv("~/MS-Thesis/data/candidate_selection/step_1/sv/gctof_SV_unadj_scaled.csv",stringsAsFactors = F)
hilic = read.csv("~/MS-Thesis/data/candidate_selection/step_1/sv/hilic_SV_unadj_scaled.csv",stringsAsFactors = F)
lipid = read.csv("~/MS-Thesis/data/candidate_selection/step_1/sv/lipid_SV_unadj_scaled.csv",stringsAsFactors = F)
oxylipin = read.csv("~/MS-Thesis/data/candidate_selection/step_1/sv/oxylipin_SV_unadj_scaled.csv",stringsAsFactors = F)
vitd = read.csv("~/MS-Thesis/data/candidate_selection/step_1/sv/vitd_SV_unadj_scaled.csv",stringsAsFactors = F)
# Get pairs
pairs = gctof[gctof$p.value < pair_cutoff,]
pairs = rbind(pairs,hilic[hilic$p.value < pair_cutoff,])
pairs = rbind(pairs,lipid[lipid$p.value < pair_cutoff,])
pairs = rbind(pairs,oxylipin[oxylipin$p.value < pair_cutoff,])
pairs = rbind(pairs,vitd[vitd$p.value < pair_cutoff,])
# Pairs with methyl or metab associated with disease
pairs = pairs[which(pairs$methyl %in% methyl_pheno | pairs$metab %in% metab_pheno),]
# Write pairs list
save(pairs,file = "~/MS-Thesis/data/networks/pair_list.Rdata")
# Combine all necessary data into one DF
pheno = read.csv("/home/biostats_share/Norris/data/phenotype/ivyomicssample_noIdentifyingInfo.csv",stringsAsFactors = F,na.strings = "")
pheno = pheno[!is.na(pheno$T1Dgroup),]
pheno = pheno[pheno$Visit_Type == "SV",]
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
methyl = methyl[,c("samplekey",unique(pairs$methyl))]
# Import metabolites and scale
gctof = read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.bc.csv",stringsAsFactors = F)
gctof = gctof[gctof$samplekey %in% pheno$samplekey,]
gctof = gctof[,c("samplekey",colnames(gctof)[which(colnames(gctof) %in% unique(pairs$metab))])]
hilic = read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.bc.csv",stringsAsFactors = F)
hilic = hilic[hilic$samplekey %in% pheno$samplekey,]
hilic = hilic[,c("samplekey","hilic_12")]
lipid = read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.bc.csv",stringsAsFactors = F)
lipid = lipid[lipid$samplekey %in% pheno$samplekey,]
lipid = lipid[,c("samplekey",colnames(lipid)[which(colnames(lipid) %in% unique(pairs$metab))])]
oxylipin = read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.bc.csv",stringsAsFactors = F)
oxylipin = oxylipin[oxylipin$samplekey %in% pheno$samplekey,]
oxylipin = oxylipin[,c("samplekey",colnames(oxylipin)[which(colnames(oxylipin) %in% unique(pairs$metab))])]
vitd = read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.bc.csv",stringsAsFactors = F)
vitd = vitd[vitd$samplekey %in% pheno$samplekey,]
vitd = vitd[,c("samplekey",colnames(vitd)[which(colnames(vitd) %in% unique(pairs$metab))])]
# Merge
df = merge(pheno,methyl,by = "samplekey",all.x = T)
df = merge(df,gctof,by = "samplekey",all.x = T)
df = merge(df,hilic,by = "samplekey",all.x = T)
df = merge(df,lipid,by = "samplekey",all.x = T)
df = merge(df,oxylipin,by = "samplekey",all.x = T)
df = merge(df,vitd,by = "samplekey",all.x = T)
df = df[!duplicated(df$samplekey),]
# Scale
vars = min(grep("cg",colnames(df))):ncol(df)
df[,vars] = lapply(df[,vars], scale)
# Save
pair_data_methyl_scaled = as.data.frame(df)
save(pair_data_methyl_scaled,file = "~/MS-Thesis/data/networks/pair_data_methyl_scaled.Rdata")
