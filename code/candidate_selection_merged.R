library(nlme)
# Data import
# Phenotype
pheno <- read.csv("/home/biostats_share/Norris/data/phenotype/ivyomicssample.csv",
                  stringsAsFactors = F,
                  na.strings = "")
# Methylation
# Combined dataset
load("/home/biostats_share/Norris/data/methylation/Mmatrix.platformAdj.Rdata")
methyl <- as.data.frame(t(M.adj))

key_450k <- read.csv("/home/biostats_share/Norris/data/methylation/key.450K.csv",
                     stringsAsFactors = F)
key_450k$platform <- "450K"
key_epic <- read.csv("/home/biostats_share/Norris/data/methylation/key.EPIC.csv",
                     stringsAsFactors = F)
key_epic$platform <- "EPIC"
key <- rbind(key_450k,key_epic)
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

# Model function with tryCatch
metab_methyl_lin_mod <- function(metabolomics,methylation,metab_name,methyl_name){
  temp <- merge(metabolomics,methylation,by = "samplekey")
  # Linear models
  out <- lapply(names(methylation)[1:5], function(x){
    base_form <- paste0(x,"~sex+age+platform")
    metabs <- lapply(names(metabolomics)[2:6], function(y) {
      form <- as.formula(paste0(base_form,"+",y))
      mod <- tryCatch(lme(form,random = ~1|samplekey,data = temp,na.action = na.omit),
                      message = function(m) NULL,warning = function(m) NULL,error = function(m) NULL)
      if (!is.null(mod)) {
        results <- as.data.frame(summary(mod)$tTable)
        results$term <- rownames(results)
        results[5,"term"] <- paste0(x,"_",y)
        return(results[5,c("term","Value","p-value")])
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
metab_methyl_lin_mod(gctof,methyl,"gctof","merged")

# hilic
metab_methyl_lin_mod(hilic,methyl,"hilic","merged")

# lipid
metab_methyl_lin_mod(lipid,methyl,"lipid","merged")

# oxylipin
metab_methyl_lin_mod(oxylipin,methyl,"oxylipin","merged")

# vitd
metab_methyl_lin_mod(vitd,methyl,"vitd","merged")
