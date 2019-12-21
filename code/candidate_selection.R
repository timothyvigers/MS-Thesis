library(nlme)
# Data import
# Phenotype
pheno <- read.csv("/home/biostats_share/Norris/data/phenotype/ivyomicssample.csv",
                  stringsAsFactors = F,
                  na.strings = "")
# Probes
load("/home/biostats_share/Norris/data/methylation/probesFromPipeline.Rdata")
probes <- as.data.frame(probesFromPipeline)
# Methylation
# 450K
load("/home/biostats_share/Norris/data/methylation/sesame450K.batchAdj.Mmatrix.Rdata")
k450 <- as.data.frame(t(M.sesame.batch))
key_450k <- read.csv("/home/biostats_share/Norris/data/methylation/key.450K.csv",
                     stringsAsFactors = F)
key_450k$samplekey[duplicated(key_450k$samplekey)] <- 
  paste0(key_450k$samplekey[duplicated(key_450k$samplekey)],".2")
k450$samplekey <- key_450k$samplekey[match(rownames(k450),key_450k$array)]
k450$age <- pheno$clinage[match(k450$samplekey,pheno$samplekey)]
k450$sex <- factor(pheno$SEX[match(k450$samplekey,pheno$samplekey)])
# EPIC
load("/home/biostats_share/Norris/data/methylation/sesameEPIC.batchAdj.Mmatrix.Rdata")
epic <- as.data.frame(t(M.sesame.batch))
key_epic <- read.csv("/home/biostats_share/Norris/data/methylation/key.EPIC.csv",
                     stringsAsFactors = F)
epic$samplekey <- key_epic$samplekey[match(rownames(epic),key_epic$array)]
epic$age <- pheno$clinage[match(epic$samplekey,pheno$samplekey)]
epic$sex <- factor(pheno$SEX[match(epic$samplekey,pheno$samplekey)]) 
rm(M.sesame.batch,probesFromPipeline)
# Metabolites
gctof <- read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.bc.csv")
hilic <- read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.bc.csv")
lipid <- read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.bc.csv")
oxylipin <- read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.bc.csv")
vitd <- read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.bc.csv")
# Annotations
anno_gctof <- read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.featureAnno.csv")
anno_hilic <- read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.featureAnno.csv")
anno_lipid <- read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.featureAnno.csv")
anno_oxylipin <- 
  read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.featureAnno.csv")
anno_vitd <- read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.featureAnno.csv")

# gctof
## 450K
temp <- merge(gctof,k450,by = "samplekey")
# Linear models
out <- lapply(names(k450)[1:5], function(x){
  base_form <- paste0(x,"~sex+age")
  metabs <- lapply(names(gctof)[2:6], function(y) {
    form <- as.formula(paste0(base_form,"+",y))
    mod <- lme(form,random = ~1|samplekey,data = temp,na.action = na.omit)
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[4,"term"] <- paste0(x,"_",y)
    results[4,c("term","Value","p-value")]
  })
  df <- do.call(rbind,metabs)
  df
})
# Combine list of DFs into large DF
out <- do.call(rbind,out)
colnames(out) <- c("term","Value","p-value")
out <- out[order(as.numeric(as.character(out$`p-value`))),]
write.csv(out,
          file = "/home/vigerst/MS-Thesis/candidate_selection/gctof_450k_results.csv",
          row.names = F)
## Epic
temp <- merge(gctof,epic,by = "samplekey")
# Linear models
out <- lapply(names(epic)[1:5], function(x){
  base_form <- paste0(x,"~sex+age")
  metabs <- lapply(names(gctof)[2:6], function(y) {
    form <- as.formula(paste0(base_form,"+",y))
    mod <- lme(form,random = ~1|samplekey,data = temp,na.action = na.omit)
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[4,"term"] <- paste0(x,"_",y)
    results[4,c("term","Value","p-value")]
  })
  df <- do.call(rbind,metabs)
  df
})
# Combine list of DFs into large DF
out <- do.call(rbind,out)
colnames(out) <- c("term","Value","p-value")
out <- out[order(as.numeric(as.character(out$`p-value`))),]
write.csv(out,
          file = "/home/vigerst/MS-Thesis/candidate_selection/gctof_epic_results.csv",
          row.names = F)
# # hilic
# ## 450K
# temp <- merge(hilic,k450,by = "samplekey")
# # Linear models
# out <- lapply(names(k450)[1:(ncol(k450)-3)], function(x){
#   base_form <- paste0(x,"~sex+age")
#   metabs <- lapply(names(hilic)[2:ncol(hilic)], function(y) {
#     form <- as.formula(paste0(base_form,"+",y))
#     mod <- lme(form,random = ~1|samplekey,data = temp,na.action = na.omit)
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"term"] <- paste0(x,"_",y)
#     results[4,c("term","Value","p-value")]
#   })
#   df <- do.call(rbind,metabs)
#   df
# })
# # Combine list of DFs into large DF
# out <- do.call(rbind,out)
# colnames(out) <- c("term","Value","p-value")
# out <- out[order(as.numeric(as.character(out$`p-value`))),]
# write.csv(out,
#           file = "/home/vigerst/MS-Thesis/candidate_selection/hilic_450k_results.csv",
#           row.names = F)
# ## Epic
# temp <- merge(hilic,epic,by = "samplekey")
# # Linear models
# out <- lapply(names(epic)[1:(ncol(epic)-3)], function(x){
#   base_form <- paste0(x,"~sex+age")
#   metabs <- lapply(names(hilic)[2:ncol(hilic)], function(y) {
#     form <- as.formula(paste0(base_form,"+",y))
#     mod <- lme(form,random = ~1|samplekey,data = temp,na.action = na.omit)
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"term"] <- paste0(x,"_",y)
#     results[4,c("term","Value","p-value")]
#   })
#   df <- do.call(rbind,metabs)
#   df
# })
# # Combine list of DFs into large DF
# out <- do.call(rbind,out)
# colnames(out) <- c("term","Value","p-value")
# out <- out[order(as.numeric(as.character(out$`p-value`))),]
# write.csv(out,
#           file = "/home/vigerst/MS-Thesis/candidate_selection/hilic_epic_results.csv",
#           row.names = F)
# 
# # lipid
# ## 450K
# temp <- merge(lipid,k450,by = "samplekey")
# # Linear models
# out <- lapply(names(k450)[1:(ncol(k450)-3)], function(x){
#   base_form <- paste0(x,"~sex+age")
#   metabs <- lapply(names(lipid)[2:ncol(lipid)], function(y) {
#     form <- as.formula(paste0(base_form,"+",y))
#     mod <- lme(form,random = ~1|samplekey,data = temp,na.action = na.omit)
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"term"] <- paste0(x,"_",y)
#     results[4,c("term","Value","p-value")]
#   })
#   df <- do.call(rbind,metabs)
#   df
# })
# # Combine list of DFs into large DF
# out <- do.call(rbind,out)
# colnames(out) <- c("term","Value","p-value")
# out <- out[order(as.numeric(as.character(out$`p-value`))),]
# write.csv(out,
#           file = "/home/vigerst/MS-Thesis/candidate_selection/lipid_450k_results.csv",
#           row.names = F)
# ## Epic
# temp <- merge(lipid,epic,by = "samplekey")
# # Linear models
# out <- lapply(names(epic)[1:(ncol(epic)-3)], function(x){
#   base_form <- paste0(x,"~sex+age")
#   metabs <- lapply(names(lipid)[2:ncol(lipid)], function(y) {
#     form <- as.formula(paste0(base_form,"+",y))
#     mod <- lme(form,random = ~1|samplekey,data = temp,na.action = na.omit)
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"term"] <- paste0(x,"_",y)
#     results[4,c("term","Value","p-value")]
#   })
#   df <- do.call(rbind,metabs)
#   df
# })
# # Combine list of DFs into large DF
# out <- do.call(rbind,out)
# colnames(out) <- c("term","Value","p-value")
# out <- out[order(as.numeric(as.character(out$`p-value`))),]
# write.csv(out,
#           file = "/home/vigerst/MS-Thesis/candidate_selection/lipid_epic_results.csv",
#           row.names = F)
# 
# # oxylipin
# ## 450K
# temp <- merge(oxylipin,k450,by = "samplekey")
# # Linear models
# out <- lapply(names(k450)[1:(ncol(k450)-3)], function(x){
#   base_form <- paste0(x,"~sex+age")
#   metabs <- lapply(names(oxylipin)[2:ncol(oxylipin)], function(y) {
#     form <- as.formula(paste0(base_form,"+",y))
#     mod <- lme(form,random = ~1|samplekey,data = temp,na.action = na.omit)
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"term"] <- paste0(x,"_",y)
#     results[4,c("term","Value","p-value")]
#   })
#   df <- do.call(rbind,metabs)
#   df
# })
# # Combine list of DFs into large DF
# out <- do.call(rbind,out)
# colnames(out) <- c("term","Value","p-value")
# out <- out[order(as.numeric(as.character(out$`p-value`))),]
# write.csv(out,
#           file = "/home/vigerst/MS-Thesis/candidate_selection/oxylipin_450k_results.csv",
#           row.names = F)
# ## Epic
# temp <- merge(oxylipin,epic,by = "samplekey")
# # Linear models
# out <- lapply(names(epic)[1:(ncol(epic)-3)], function(x){
#   base_form <- paste0(x,"~sex+age")
#   metabs <- lapply(names(oxylipin)[2:ncol(oxylipin)], function(y) {
#     form <- as.formula(paste0(base_form,"+",y))
#     mod <- lme(form,random = ~1|samplekey,data = temp,na.action = na.omit)
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"term"] <- paste0(x,"_",y)
#     results[4,c("term","Value","p-value")]
#   })
#   df <- do.call(rbind,metabs)
#   df
# })
# # Combine list of DFs into large DF
# out <- do.call(rbind,out)
# colnames(out) <- c("term","Value","p-value")
# out <- out[order(as.numeric(as.character(out$`p-value`))),]
# write.csv(out,
#           file = "/home/vigerst/MS-Thesis/candidate_selection/oxylipin_epic_results.csv",
#           row.names = F)
# 
# # vitd
# ## 450K
# temp <- merge(vitd,k450,by = "samplekey")
# # Linear models
# out <- lapply(names(k450)[1:(ncol(k450)-3)], function(x){
#   base_form <- paste0(x,"~sex+age")
#   metabs <- lapply(names(vitd)[2:ncol(vitd)], function(y) {
#     form <- as.formula(paste0(base_form,"+",y))
#     mod <- lme(form,random = ~1|samplekey,data = temp,na.action = na.omit)
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"term"] <- paste0(x,"_",y)
#     results[4,c("term","Value","p-value")]
#   })
#   df <- do.call(rbind,metabs)
#   df
# })
# # Combine list of DFs into large DF
# out <- do.call(rbind,out)
# colnames(out) <- c("term","Value","p-value")
# out <- out[order(as.numeric(as.character(out$`p-value`))),]
# write.csv(out,
#           file = "/home/vigerst/MS-Thesis/candidate_selection/vitd_450k_results.csv",
#           row.names = F)
# ## Epic
# temp <- merge(vitd,epic,by = "samplekey")
# # Linear models
# out <- lapply(names(epic)[1:(ncol(epic)-3)], function(x){
#   base_form <- paste0(x,"~sex+age")
#   metabs <- lapply(names(vitd)[2:ncol(vitd)], function(y) {
#     form <- as.formula(paste0(base_form,"+",y))
#     mod <- lme(form,random = ~1|samplekey,data = temp,na.action = na.omit)
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"term"] <- paste0(x,"_",y)
#     results[4,c("term","Value","p-value")]
#   })
#   df <- do.call(rbind,metabs)
#   df
# })
# # Combine list of DFs into large DF
# out <- do.call(rbind,out)
# colnames(out) <- c("term","Value","p-value")
# out <- out[order(as.numeric(as.character(out$`p-value`))),]
# write.csv(out,
#           file = "/home/vigerst/MS-Thesis/candidate_selection/vitd_epic_results.csv",
#           row.names = F)
