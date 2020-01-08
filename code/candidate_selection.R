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
rm(M.sesame.batch)
# Metabolites
gctof <- read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.bc.csv")
hilic <- read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.bc.csv")
lipid <- read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.bc.csv")
oxylipin <- read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.bc.csv")
vitd <- read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.bc.csv")

# 450k 
# gctof
temp <- merge(gctof,k450,by = "samplekey")
methyl <- names(k450)[1:(ncol(k450)-3)]
metab <- names(gctof)[2:ncol(gctof)]
mods <- paste0(methyl,"~sex+age")
mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods,function(x){
  form <- as.formula(x)
  mod <- lme(form,random = ~1|samplekey,data = temp,
             na.action = na.omit)
  if (!is.null(mod)) {
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
    results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
    return(results[4,c("methyl","metab","Value","p-value")])
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
# hilic 
temp <- merge(hilic,k450,by = "samplekey")
methyl <- names(k450)[1:(ncol(k450)-3)]
metab <- names(hilic)[2:ncol(hilic)]
mods <- paste0(methyl,"~sex+age")
mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods,function(x){
  form <- as.formula(x)
  mod <- lme(form,random = ~1|samplekey,data = temp,
             na.action = na.omit)
  if (!is.null(mod)) {
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
    results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
    return(results[4,c("methyl","metab","Value","p-value")])
  } else {
    results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
    colnames(results) <- c("methyl","metab","Value","p-value")
    return(results)
  }
})
df <- do.call(rbind,result_list)
filename <- paste0(out_dir,"/","450k","_","hilic","_parallel.csv")
write.csv(df,file = filename,row.names = F)
stopCluster(cl)
# lipid 
temp <- merge(lipid,k450,by = "samplekey")
methyl <- names(k450)[1:(ncol(k450)-3)]
metab <- names(lipid)[2:ncol(lipid)]
mods <- paste0(methyl,"~sex+age")
mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods,function(x){
  form <- as.formula(x)
  mod <- lme(form,random = ~1|samplekey,data = temp,
             na.action = na.omit)
  if (!is.null(mod)) {
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
    results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
    return(results[4,c("methyl","metab","Value","p-value")])
  } else {
    results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
    colnames(results) <- c("methyl","metab","Value","p-value")
    return(results)
  }
})
df <- do.call(rbind,result_list)
filename <- paste0(out_dir,"/","450k","_","lipid","_parallel.csv")
write.csv(df,file = filename,row.names = F)
stopCluster(cl)
# oxylipin 
temp <- merge(oxylipin,k450,by = "samplekey")
methyl <- names(k450)[1:(ncol(k450)-3)]
metab <- names(oxylipin)[2:ncol(oxylipin)]
mods <- paste0(methyl,"~sex+age")
mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods,function(x){
  form <- as.formula(x)
  mod <- lme(form,random = ~1|samplekey,data = temp,
             na.action = na.omit)
  if (!is.null(mod)) {
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
    results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
    return(results[4,c("methyl","metab","Value","p-value")])
  } else {
    results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
    colnames(results) <- c("methyl","metab","Value","p-value")
    return(results)
  }
})
df <- do.call(rbind,result_list)
filename <- paste0(out_dir,"/","450k","_","oxylipin","_parallel.csv")
write.csv(df,file = filename,row.names = F)
stopCluster(cl)
# vitd 
temp <- merge(vitd,k450,by = "samplekey")
methyl <- names(k450)[1:(ncol(k450)-3)]
metab <- names(vitd)[2:ncol(vitd)]
mods <- paste0(methyl,"~sex+age")
mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods,function(x){
  form <- as.formula(x)
  mod <- lme(form,random = ~1|samplekey,data = temp,
             na.action = na.omit)
  if (!is.null(mod)) {
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
    results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
    return(results[4,c("methyl","metab","Value","p-value")])
  } else {
    results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
    colnames(results) <- c("methyl","metab","Value","p-value")
    return(results)
  }
})
df <- do.call(rbind,result_list)
filename <- paste0(out_dir,"/","450k","_","vitd","_parallel.csv")
write.csv(df,file = filename,row.names = F)
stopCluster(cl)

# EPIC
# gctof
temp <- merge(gctof,epic,by = "samplekey")
methyl <- names(epic)[1:(ncol(epic)-3)]
metab <- names(gctof)[2:ncol(gctof)]
mods <- paste0(methyl,"~sex+age")
mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods,function(x){
  form <- as.formula(x)
  mod <- lme(form,random = ~1|samplekey,data = temp,
             na.action = na.omit)
  if (!is.null(mod)) {
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
    results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
    return(results[4,c("methyl","metab","Value","p-value")])
  } else {
    results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
    colnames(results) <- c("methyl","metab","Value","p-value")
    return(results)
  }
})
df <- do.call(rbind,result_list)
filename <- paste0(out_dir,"/","EPIC","_","gctof","_parallel.csv")
write.csv(df,file = filename,row.names = F)
stopCluster(cl)
# hilic
temp <- merge(hilic,epic,by = "samplekey")
methyl <- names(epic)[1:(ncol(epic)-3)]
metab <- names(hilic)[2:ncol(hilic)]
mods <- paste0(methyl,"~sex+age")
mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods,function(x){
  form <- as.formula(x)
  mod <- lme(form,random = ~1|samplekey,data = temp,
             na.action = na.omit)
  if (!is.null(mod)) {
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
    results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
    return(results[4,c("methyl","metab","Value","p-value")])
  } else {
    results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
    colnames(results) <- c("methyl","metab","Value","p-value")
    return(results)
  }
})
df <- do.call(rbind,result_list)
filename <- paste0(out_dir,"/","EPIC","_","hilic","_parallel.csv")
write.csv(df,file = filename,row.names = F)
stopCluster(cl)
# lipid
temp <- merge(lipid,epic,by = "samplekey")
methyl <- names(epic)[1:(ncol(epic)-3)]
metab <- names(lipid)[2:ncol(lipid)]
mods <- paste0(methyl,"~sex+age")
mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods,function(x){
  form <- as.formula(x)
  mod <- lme(form,random = ~1|samplekey,data = temp,
             na.action = na.omit)
  if (!is.null(mod)) {
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
    results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
    return(results[4,c("methyl","metab","Value","p-value")])
  } else {
    results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
    colnames(results) <- c("methyl","metab","Value","p-value")
    return(results)
  }
})
df <- do.call(rbind,result_list)
filename <- paste0(out_dir,"/","EPIC","_","lipid","_parallel.csv")
write.csv(df,file = filename,row.names = F)
stopCluster(cl)
# oxylipin
temp <- merge(oxylipin,epic,by = "samplekey")
methyl <- names(epic)[1:(ncol(epic)-3)]
metab <- names(oxylipin)[2:ncol(oxylipin)]
mods <- paste0(methyl,"~sex+age")
mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods,function(x){
  form <- as.formula(x)
  mod <- lme(form,random = ~1|samplekey,data = temp,
             na.action = na.omit)
  if (!is.null(mod)) {
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
    results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
    return(results[4,c("methyl","metab","Value","p-value")])
  } else {
    results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
    colnames(results) <- c("methyl","metab","Value","p-value")
    return(results)
  }
})
df <- do.call(rbind,result_list)
filename <- paste0(out_dir,"/","EPIC","_","oxylipin","_parallel.csv")
write.csv(df,file = filename,row.names = F)
stopCluster(cl)
# vitd
temp <- merge(vitd,epic,by = "samplekey")
methyl <- names(epic)[1:(ncol(epic)-3)]
metab <- names(vitd)[2:ncol(vitd)]
mods <- paste0(methyl,"~sex+age")
mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
# Make cluster
cl <- makeCluster(no_cores,type = "FORK")
# Linear models
result_list <- parLapply(cl,mods,function(x){
  form <- as.formula(x)
  mod <- lme(form,random = ~1|samplekey,data = temp,
             na.action = na.omit)
  if (!is.null(mod)) {
    results <- as.data.frame(summary(mod)$tTable)
    results$term <- rownames(results)
    results[4,"methyl"] <- strsplit(x,"~")[[1]][1]
    results[4,"metab"] <- strsplit(x,"\\+")[[1]][3]
    return(results[4,c("methyl","metab","Value","p-value")])
  } else {
    results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
    colnames(results) <- c("methyl","metab","Value","p-value")
    return(results)
  }
})
df <- do.call(rbind,result_list)
filename <- paste0(out_dir,"/","EPIC","_","vitd","_parallel.csv")
write.csv(df,file = filename,row.names = F)
stopCluster(cl)