library(parallel)
# Parallel package cores
no_cores <- 80
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

# Model function with tryCatch
metab_methyl_lin_mod <- function(metabolomics,methylation,metab_name,methyl_name){
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, library(nlme))
  temp <- merge(metabolomics,methylation,by = "samplekey")
  # Linear models
  methyl <- names(methylation)[1:(ncol(methylation)-3)]
  metab <- names(metabolomics)[2:ncol(metabolomics)]
  mods <- paste0(methyl,"~sex+age")
  mods <- paste(rep(mods, each = length(metab)), metab, sep = "+")
  clusterExport(cl,c("temp","mods"))
  result_list <- parLapply(cl,mods[1:100],function(x){
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
  return(df)
  stopCluster(cl)
}

# gctof
## 450K
metab_methyl_lin_mod(gctof,k450,"gctof","450K")
# ## Epic
# metab_methyl_lin_mod(gctof,epic,"gctof","EPIC")
# 
# # hilic
# ## 450K
# metab_methyl_lin_mod(hilic,k450,"hilic","450K")
# ## Epic
# metab_methyl_lin_mod(hilic,epic,"hilic","EPIC")
# 
# # lipid
# ## 450K
# metab_methyl_lin_mod(lipid,k450,"lipid","450K")
# ## Epic
# metab_methyl_lin_mod(lipid,epic,"lipid","EPIC")
# 
# # oxylipin
# ## 450K
# metab_methyl_lin_mod(oxylipin,k450,"oxylipin","450K")
# ## Epic
# metab_methyl_lin_mod(oxylipin,epic,"oxylipin","EPIC")
# 
# # vitd
# ## 450K
# metab_methyl_lin_mod(vitd,k450,"vitd","450K")
# ## Epic
# metab_methyl_lin_mod(vitd,epic,"vitd","EPIC")
