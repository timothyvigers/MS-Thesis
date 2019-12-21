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

# Model function with tryCatch
metab_methyl_lin_mod <- function(metabolomics,methylation,metab_name,methyl_name){
  temp <- merge(metabolomics,methylation,by = "samplekey")
  # Linear models
  out <- lapply(names(methylation)[1:5], function(x){
    base_form <- paste0(x,"~sex+age")
    metabs <- lapply(names(metabolomics)[2:6], function(y) {
      form <- as.formula(paste0(base_form,"+",y))
      mod <- tryCatch(lme(form,random = ~1|samplekey,data = temp,na.action = na.omit),
                      message = function(m) NULL,warning = function(m) NULL,error = function(m) NULL)
      if (!is.null(mod)) {
        results <- as.data.frame(summary(mod)$tTable)
        results$term <- rownames(results)
        results[4,"term"] <- paste0(x,"_",y)
        return(results[4,c("term","Value","p-value")])
      } else {
        results <- as.data.frame(matrix(c(NA,NA,NA),nrow = 1))
        colnames(results) <- c("term","Value","p-value")
        return(results)
      }
    })
    df <- do.call(rbind,metabs)
    df
  })
  # Combine list of DFs into large DF
  out <- do.call(rbind,out)
  out <- out[order(as.numeric(as.character(out$`p-value`))),]
  file <- paste0("/home/vigerst/MS-Thesis/candidate_selection/",metab_name,"_",
                 methyl_name,"_results.csv")
  write.csv(out[complete.cases(out),],file = file,row.names = F)
}

# gctof
## 450K
metab_methyl_lin_mod(gctof,k450,"gctof","450K")
## Epic
metab_methyl_lin_mod(gctof,epic,"gctof","EPIC")

# hilic
## 450K
metab_methyl_lin_mod(hilic,k450,"hilic","450K")
## Epic
metab_methyl_lin_mod(hilic,epic,"hilic","EPIC")

# lipid
## 450K
metab_methyl_lin_mod(lipid,k450,"lipid","450K")
## Epic
metab_methyl_lin_mod(lipid,epic,"lipid","EPIC")

# oxylipin
## 450K
metab_methyl_lin_mod(oxylipin,k450,"oxylipin","450K")
## Epic
metab_methyl_lin_mod(oxylipin,epic,"oxylipin","EPIC")

# vitd
## 450K
metab_methyl_lin_mod(vitd,k450,"vitd","450K")
## Epic
metab_methyl_lin_mod(vitd,epic,"vitd","EPIC")
