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
methyl$samplekey = key$samplekey[match(rownames(methyl),key$array)]
methyl$platform = key$platform[match(rownames(methyl),key$array)]
methyl = methyl[match(pheno$samplekey,methyl$samplekey),]
methyl$id = factor(pheno$ID[match(methyl$samplekey,pheno$samplekey)])
methyl = methyl[,c(probesFromPipeline,"samplekey","id","platform")]
# Import metabolites and scale
gctof <- read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.bc.csv",
                  stringsAsFactors = F)
gctof <- gctof[gctof$samplekey %in% pheno$samplekey,]
gctof[,2:ncol(gctof)] <- lapply(gctof[,2:ncol(gctof)],scale)
hilic <- read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.bc.csv",
                  stringsAsFactors = F)
hilic <- hilic[hilic$samplekey %in% pheno$samplekey,]
hilic[,2:ncol(hilic)] <- lapply(hilic[,2:ncol(hilic)],scale)
lipid <- read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.bc.csv",
                  stringsAsFactors = F)
lipid <- lipid[lipid$samplekey %in% pheno$samplekey,]
lipid[,2:ncol(lipid)] <- lapply(lipid[,2:ncol(lipid)],scale)
oxylipin <- read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.bc.csv",
                     stringsAsFactors = F)
oxylipin <- oxylipin[oxylipin$samplekey %in% pheno$samplekey,]
oxylipin[,2:ncol(oxylipin)] <- lapply(oxylipin[,2:ncol(oxylipin)],scale)
vitd <- read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.bc.csv",
                 stringsAsFactors = F)
vitd <- vitd[vitd$samplekey %in% pheno$samplekey,]
vitd[,2:ncol(vitd)] <- lapply(vitd[,2:ncol(vitd)],scale)

# Model function
run_mods <- function(mods = model_list, data = temp,metabname,no_cores = 10,
                     out_dir = "/home/vigerst/MS-Thesis/data/candidate_selection/step_1/sv/") {
  require(parallel)
  # Make cluster
  cl <- makeCluster(no_cores,type = "FORK")
  # Parallel models
  result_list <- parLapply(cl,mods,function(x){
    form <- as.formula(x)
    mod <- tryCatch(lm(form,data = data),
                    message = function(m) NULL,warning = function(m) NULL,
                    error = function(m) NULL)
    if (!is.null(mod)) {
      results <- as.data.frame(summary(mod)$coefficients)
      results$term <- rownames(results)
      results[nrow(results),"methyl"] <- strsplit(x,"~")[[1]][1]
      results[nrow(results),"metab"] <- strsplit(x,"\\+")[[1]][2]
      results <- results[nrow(results),c("methyl","metab","Estimate","Pr(>|t|)")]
      colnames(results) <- c("methyl","metab","Value","p-value")
      return(results)
    } else {
      results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
      colnames(results) <- c("methyl","metab","Value","p-value")
      return(results)
    }
  })
  df <- do.call(rbind,result_list)
  filename <- paste0(out_dir,metabname,"_SV_unadj_scaled.csv")
  write.csv(df,file = filename,row.names = F)
  stopCluster(cl)
}

# gctof
temp <- merge(gctof,methyl,by = "samplekey")
metab <- names(gctof)[2:ncol(gctof)]
model_list <- paste0(rep(probesFromPipeline,each = length(metab)),"~ platform+",metab)

run_mods(model_list,metabname = "gctof")

# hilic
temp <- merge(hilic,methyl,by = "samplekey")
metab <- names(hilic)[2:ncol(hilic)]
model_list <- paste0(rep(probesFromPipeline,each = length(metab)),"~ platform+",metab)

run_mods(model_list,metabname = "hilic")

# lipid
temp <- merge(lipid,methyl,by = "samplekey")
metab <- names(lipid)[2:ncol(lipid)]
model_list <- paste0(rep(probesFromPipeline,each = length(metab)),"~ platform+",metab)

run_mods(model_list,metabname = "lipid")

# oxylipin
temp <- merge(oxylipin,methyl,by = "samplekey")
metab <- names(oxylipin)[2:ncol(oxylipin)]
model_list <- paste0(rep(probesFromPipeline,each = length(metab)),"~ platform+",metab)

run_mods(model_list,metabname = "oxylipin")

# vitd
temp <- merge(vitd,methyl,by = "samplekey")
metab <- names(vitd)[2:ncol(vitd)]
model_list <- paste0(rep(probesFromPipeline,each = length(metab)),"~ platform+",metab)

run_mods(model_list,metabname = "vitd")