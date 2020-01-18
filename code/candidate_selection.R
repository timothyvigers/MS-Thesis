library(parallel)
library(nlme)
# Parallel package cores
no_cores <- 60
# Output
out_dir = "/home/vigerst/MS-Thesis/candidate_selection/"
# Data import
# Phenotype
pheno <- read.csv("/home/biostats_share/Norris/data/phenotype/ivyomicssample.csv",
                  stringsAsFactors = F,
                  na.strings = "")
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
methyl$age <- pheno$clinage[match(methyl$samplekey,pheno$samplekey)]
methyl$sex <- factor(pheno$SEX[match(methyl$samplekey,pheno$samplekey)])
methyl$platform <- key$platform[match(methyl$samplekey,pheno$samplekey)]
# Metabolites
gctof <- read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.bc.csv")
hilic <- read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.bc.csv")
lipid <- read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.bc.csv")
oxylipin <- read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.bc.csv")
vitd <- read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.bc.csv")

# # gctof
# temp <- merge(gctof,methyl,by = "samplekey")
# metab <- names(gctof)[2:ncol(gctof)]
# mods <- paste0(probesFromPipeline,"~sex+age+platform")
# mods <- paste0(rep(mods,each = length(metab)),metab)
# # Make cluster
# cl <- makeCluster(no_cores,type = "FORK")
# # Linear models
# result_list <- parLapply(cl,mods,function(x){
#   form <- as.formula(x)
#   mod <- tryCatch(lme(form,random = ~1|samplekey,data = temp,na.action = na.omit),
#                   message = function(m) NULL,warning = function(m) NULL,
#                   error = function(m) NULL)
#   if (!is.null(mod)) {
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[5,"methyl"] <- strsplit(x,"~")[[1]][1]
#     results[5,"metab"] <- strsplit(x,"\\+")[[1]][4]
#     return(results[5,c("methyl","metab","Value","p-value")])
#   } else {
#     results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
#     colnames(results) <- c("methyl","metab","Value","p-value")
#     return(results)
#   }
# })
# df <- do.call(rbind,result_list)
# filename <- paste0(out_dir,"gctof","_parallel.csv")
# write.csv(df,file = filename,row.names = F)
# stopCluster(cl)
# 
# # hilic
# temp <- merge(hilic,methyl,by = "samplekey")
# metab <- names(hilic)[2:ncol(hilic)]
# mods <- paste0(probesFromPipeline,"~sex+age+platform")
# mods <- paste0(rep(mods,each = length(metab)),metab)
# # Make cluster
# cl <- makeCluster(no_cores,type = "FORK")
# # Linear models
# result_list <- parLapply(cl,mods,function(x){
#   form <- as.formula(x)
#   mod <- tryCatch(lme(form,random = ~1|samplekey,data = temp,na.action = na.omit),
#                   message = function(m) NULL,warning = function(m) NULL,
#                   error = function(m) NULL)
#   if (!is.null(mod)) {
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[5,"methyl"] <- strsplit(x,"~")[[1]][1]
#     results[5,"metab"] <- strsplit(x,"\\+")[[1]][4]
#     return(results[5,c("methyl","metab","Value","p-value")])
#   } else {
#     results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
#     colnames(results) <- c("methyl","metab","Value","p-value")
#     return(results)
#   }
# })
# df <- do.call(rbind,result_list)
# filename <- paste0(out_dir,"hilic","_parallel.csv")
# write.csv(df,file = filename,row.names = F)
# stopCluster(cl)
# 
# # lipid
# temp <- merge(lipid,methyl,by = "samplekey")
# metab <- names(lipid)[2:ncol(lipid)]
# mods <- paste0(probesFromPipeline,"~sex+age+platform")
# mods <- paste0(rep(mods,each = length(metab)),metab)
# # Make cluster
# cl <- makeCluster(no_cores,type = "FORK")
# # Linear models
# result_list <- parLapply(cl,mods,function(x){
#   form <- as.formula(x)
#   mod <- tryCatch(lme(form,random = ~1|samplekey,data = temp,na.action = na.omit),
#                   message = function(m) NULL,warning = function(m) NULL,
#                   error = function(m) NULL)
#   if (!is.null(mod)) {
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[5,"methyl"] <- strsplit(x,"~")[[1]][1]
#     results[5,"metab"] <- strsplit(x,"\\+")[[1]][4]
#     return(results[5,c("methyl","metab","Value","p-value")])
#   } else {
#     results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
#     colnames(results) <- c("methyl","metab","Value","p-value")
#     return(results)
#   }
# })
# df <- do.call(rbind,result_list)
# filename <- paste0(out_dir,"lipid","_parallel.csv")
# write.csv(df,file = filename,row.names = F)
# stopCluster(cl)

# oxylipin
temp <- merge(oxylipin,methyl,by = "samplekey")
metab <- names(oxylipin)[2:ncol(oxylipin)]
mods <- paste0(probesFromPipeline,"~sex+age+platform")
mods <- paste0(rep(mods,each = length(metab)),metab)
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods,function(x){
  form <- as.formula(x)
  mod <- tryCatch(lme(form,random = ~1|samplekey,data = temp,na.action = na.omit),
                  message = function(m) NULL,warning = function(m) NULL,
                  error = function(m) NULL)
  if (!is.null(mod)) {
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[5,"methyl"] <- strsplit(x,"~")[[1]][1]
    results[5,"metab"] <- strsplit(x,"\\+")[[1]][4]
    return(results[5,c("methyl","metab","Value","p-value")])
  } else {
    results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
    colnames(results) <- c("methyl","metab","Value","p-value")
    return(results)
  }
})
df <- do.call(rbind,result_list)
filename <- paste0(out_dir,"oxylipin","_parallel.csv")
write.csv(df,file = filename,row.names = F)
stopCluster(cl)

# vitd
temp <- merge(vitd,methyl,by = "samplekey")
metab <- names(vitd)[2:ncol(vitd)]
mods <- paste0(probesFromPipeline,"~sex+age+platform")
mods <- paste0(rep(mods,each = length(metab)),metab)
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods,function(x){
  form <- as.formula(x)
  mod <- tryCatch(lme(form,random = ~1|samplekey,data = temp,na.action = na.omit),
                  message = function(m) NULL,warning = function(m) NULL,
                  error = function(m) NULL)
  if (!is.null(mod)) {
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[5,"methyl"] <- strsplit(x,"~")[[1]][1]
    results[5,"metab"] <- strsplit(x,"\\+")[[1]][4]
    return(results[5,c("methyl","metab","Value","p-value")])
  } else {
    results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
    colnames(results) <- c("methyl","metab","Value","p-value")
    return(results)
  }
})
df <- do.call(rbind,result_list)
filename <- paste0(out_dir,"vitd","_parallel.csv")
write.csv(df,file = filename,row.names = F)
stopCluster(cl)