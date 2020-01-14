library(parallel)
library(nlme)
# Parallel package cores
no_cores <- 60
# Output
out_dir = "/home/vigerst/MS-Thesis/candidate_selection"
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
# 450k 
# gctof
temp <- merge(gctof,methyl,by = "samplekey")
metab <- names(gctof)[2:ncol(gctof)]
mods <- paste0(names(methyl)[1:(length(names(methyl))-4)],"~sex+age+platform+",metab[1])
#mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods[1:100],function(x){
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
filename <- paste0(out_dir,"/","450k","_","gctof","_parallel.csv")
write.csv(df,file = filename,row.names = F)
stopCluster(cl)

# Use Mmatrix.platform adjust, subset the probes
# Adjust for platform as a covariate?q
# Ask Theresa for candidate metabolites 


# # hilic
# temp <- merge(hilic,k450,by = "samplekey")
# methyl <- names(k450)[1:(ncol(k450)-3)]
# metab <- names(hilic)[2:ncol(hilic)]
# mods <- paste0(methyl,"~sex+age")
# mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# # Make cluster
# cl <- makeCluster(no_cores,type = "FORK")
# # Linear models
# result_list <- parLapply(cl,mods,function(x){
#   form <- as.formula(x)
#   mod <- lme(form,random = ~1|samplekey,data = temp,
#              na.action = na.omit)
#   if (!is.null(mod)) {
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
#     results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
#     return(results[4,c("methyl","metab","Value","p-value")])
#   } else {
#     results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
#     colnames(results) <- c("methyl","metab","Value","p-value")
#     return(results)
#   }
# })
# df <- do.call(rbind,result_list)
# filename <- paste0(out_dir,"/","450k","_","hilic","_parallel.csv")
# write.csv(df,file = filename,row.names = F)
# stopCluster(cl)
# # lipid
# temp <- merge(lipid,k450,by = "samplekey")
# methyl <- names(k450)[1:(ncol(k450)-3)]
# metab <- names(lipid)[2:ncol(lipid)]
# mods <- paste0(methyl,"~sex+age")
# mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# # Make cluster
# cl <- makeCluster(no_cores,type = "FORK")
# # Linear models
# result_list <- parLapply(cl,mods,function(x){
#   form <- as.formula(x)
#   mod <- lme(form,random = ~1|samplekey,data = temp,
#              na.action = na.omit)
#   if (!is.null(mod)) {
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
#     results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
#     return(results[4,c("methyl","metab","Value","p-value")])
#   } else {
#     results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
#     colnames(results) <- c("methyl","metab","Value","p-value")
#     return(results)
#   }
# })
# df <- do.call(rbind,result_list)
# filename <- paste0(out_dir,"/","450k","_","lipid","_parallel.csv")
# write.csv(df,file = filename,row.names = F)
# stopCluster(cl)
# # oxylipin
# temp <- merge(oxylipin,k450,by = "samplekey")
# methyl <- names(k450)[1:(ncol(k450)-3)]
# metab <- names(oxylipin)[2:ncol(oxylipin)]
# mods <- paste0(methyl,"~sex+age")
# mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# # Make cluster
# cl <- makeCluster(no_cores,type = "FORK")
# # Linear models
# result_list <- parLapply(cl,mods,function(x){
#   form <- as.formula(x)
#   mod <- lme(form,random = ~1|samplekey,data = temp,
#              na.action = na.omit)
#   if (!is.null(mod)) {
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
#     results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
#     return(results[4,c("methyl","metab","Value","p-value")])
#   } else {
#     results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
#     colnames(results) <- c("methyl","metab","Value","p-value")
#     return(results)
#   }
# })
# df <- do.call(rbind,result_list)
# filename <- paste0(out_dir,"/","450k","_","oxylipin","_parallel.csv")
# write.csv(df,file = filename,row.names = F)
# stopCluster(cl)
# # vitd
# temp <- merge(vitd,k450,by = "samplekey")
# methyl <- names(k450)[1:(ncol(k450)-3)]
# metab <- names(vitd)[2:ncol(vitd)]
# mods <- paste0(methyl,"~sex+age")
# mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# # Make cluster
# cl <- makeCluster(no_cores,type = "FORK")
# # Linear models
# result_list <- parLapply(cl,mods,function(x){
#   form <- as.formula(x)
#   mod <- lme(form,random = ~1|samplekey,data = temp,
#              na.action = na.omit)
#   if (!is.null(mod)) {
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
#     results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
#     return(results[4,c("methyl","metab","Value","p-value")])
#   } else {
#     results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
#     colnames(results) <- c("methyl","metab","Value","p-value")
#     return(results)
#   }
# })
# df <- do.call(rbind,result_list)
# filename <- paste0(out_dir,"/","450k","_","vitd","_parallel.csv")
# write.csv(df,file = filename,row.names = F)
# stopCluster(cl)
# 
# # EPIC
# # gctof
# temp <- merge(gctof,epic,by = "samplekey")
# methyl <- names(epic)[1:(ncol(epic)-3)]
# metab <- names(gctof)[2:ncol(gctof)]
# mods <- paste0(methyl,"~sex+age")
# mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# # Make cluster
# cl <- makeCluster(no_cores,type = "FORK")
# # Linear models
# result_list <- parLapply(cl,mods,function(x){
#   form <- as.formula(x)
#   mod <- lme(form,random = ~1|samplekey,data = temp,
#              na.action = na.omit)
#   if (!is.null(mod)) {
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
#     results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
#     return(results[4,c("methyl","metab","Value","p-value")])
#   } else {
#     results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
#     colnames(results) <- c("methyl","metab","Value","p-value")
#     return(results)
#   }
# })
# df <- do.call(rbind,result_list)
# filename <- paste0(out_dir,"/","EPIC","_","gctof","_parallel.csv")
# write.csv(df,file = filename,row.names = F)
# stopCluster(cl)
# # hilic
# temp <- merge(hilic,epic,by = "samplekey")
# methyl <- names(epic)[1:(ncol(epic)-3)]
# metab <- names(hilic)[2:ncol(hilic)]
# mods <- paste0(methyl,"~sex+age")
# mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# # Make cluster
# cl <- makeCluster(no_cores,type = "FORK")
# # Linear models
# result_list <- parLapply(cl,mods,function(x){
#   form <- as.formula(x)
#   mod <- lme(form,random = ~1|samplekey,data = temp,
#              na.action = na.omit)
#   if (!is.null(mod)) {
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
#     results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
#     return(results[4,c("methyl","metab","Value","p-value")])
#   } else {
#     results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
#     colnames(results) <- c("methyl","metab","Value","p-value")
#     return(results)
#   }
# })
# df <- do.call(rbind,result_list)
# filename <- paste0(out_dir,"/","EPIC","_","hilic","_parallel.csv")
# write.csv(df,file = filename,row.names = F)
# stopCluster(cl)
# # lipid
# temp <- merge(lipid,epic,by = "samplekey")
# methyl <- names(epic)[1:(ncol(epic)-3)]
# metab <- names(lipid)[2:ncol(lipid)]
# mods <- paste0(methyl,"~sex+age")
# mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# # Make cluster
# cl <- makeCluster(no_cores,type = "FORK")
# # Linear models
# result_list <- parLapply(cl,mods,function(x){
#   form <- as.formula(x)
#   mod <- lme(form,random = ~1|samplekey,data = temp,
#              na.action = na.omit)
#   if (!is.null(mod)) {
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
#     results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
#     return(results[4,c("methyl","metab","Value","p-value")])
#   } else {
#     results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
#     colnames(results) <- c("methyl","metab","Value","p-value")
#     return(results)
#   }
# })
# df <- do.call(rbind,result_list)
# filename <- paste0(out_dir,"/","EPIC","_","lipid","_parallel.csv")
# write.csv(df,file = filename,row.names = F)
# stopCluster(cl)
# # oxylipin
# temp <- merge(oxylipin,epic,by = "samplekey")
# methyl <- names(epic)[1:(ncol(epic)-3)]
# metab <- names(oxylipin)[2:ncol(oxylipin)]
# mods <- paste0(methyl,"~sex+age")
# mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# # Make cluster
# cl <- makeCluster(no_cores,type = "FORK")
# # Linear models
# result_list <- parLapply(cl,mods,function(x){
#   form <- as.formula(x)
#   mod <- lme(form,random = ~1|samplekey,data = temp,
#              na.action = na.omit)
#   if (!is.null(mod)) {
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
#     results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
#     return(results[4,c("methyl","metab","Value","p-value")])
#   } else {
#     results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
#     colnames(results) <- c("methyl","metab","Value","p-value")
#     return(results)
#   }
# })
# df <- do.call(rbind,result_list)
# filename <- paste0(out_dir,"/","EPIC","_","oxylipin","_parallel.csv")
# write.csv(df,file = filename,row.names = F)
# stopCluster(cl)
# # vitd
# temp <- merge(vitd,epic,by = "samplekey")
# methyl <- names(epic)[1:(ncol(epic)-3)]
# metab <- names(vitd)[2:ncol(vitd)]
# mods <- paste0(methyl,"~sex+age")
# mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# # Make cluster
# cl <- makeCluster(no_cores,type = "FORK")
# # Linear models
# result_list <- parLapply(cl,mods,function(x){
#   form <- as.formula(x)
#   mod <- lme(form,random = ~1|samplekey,data = temp,
#              na.action = na.omit)
#   if (!is.null(mod)) {
#     results <- as.data.frame(summary(mod)$tTable)
#     results$term <- rownames(results)
#     results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
#     results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
#     return(results[4,c("methyl","metab","Value","p-value")])
#   } else {
#     results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
#     colnames(results) <- c("methyl","metab","Value","p-value")
#     return(results)
#   }
# })
# df <- do.call(rbind,result_list)
# filename <- paste0(out_dir,"/","EPIC","_","vitd","_parallel.csv")
# write.csv(df,file = filename,row.names = F)
# stopCluster(cl)