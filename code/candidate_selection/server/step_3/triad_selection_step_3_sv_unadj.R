# Data import
# Phenotype
pheno <- read.csv("/home/biostats_share/Norris/data/phenotype/ivyomicssample_noIdentifyingInfo.csv",
                  stringsAsFactors = F,
                  na.strings = "")
pheno$T1Dgroup <- as.factor(pheno$T1Dgroup)
pheno <- pheno[with(pheno,order(ID,DOVISIT)),]
pheno <- pheno[!is.na(pheno$T1Dgroup),]
pheno <- pheno[pheno$Visit_Type == "SV",]
# Metabolites
gctof <- read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.bc.csv")
gctof <- gctof[gctof$samplekey %in% pheno$samplekey,]
gctof[,2:ncol(gctof)] <- lapply(gctof[,2:ncol(gctof)],scale)
hilic <- read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.bc.csv")
hilic <- hilic[hilic$samplekey %in% pheno$samplekey,]
hilic[,2:ncol(hilic)] <- lapply(hilic[,2:ncol(hilic)],scale)
lipid <- read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.bc.csv")
lipid <- lipid[lipid$samplekey %in% pheno$samplekey,]
lipid[,2:ncol(lipid)] <- lapply(lipid[,2:ncol(lipid)],scale)
oxylipin <- read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.bc.csv")
oxylipin <- oxylipin[oxylipin$samplekey %in% pheno$samplekey,]
oxylipin[,2:ncol(oxylipin)] <- lapply(oxylipin[,2:ncol(oxylipin)],scale)
vitd <- read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.bc.csv")
vitd <- vitd[vitd$samplekey %in% pheno$samplekey,]
vitd[,2:ncol(vitd)] <- lapply(vitd[,2:ncol(vitd)],scale)
# Model function
run_mods <- function(mods = model_list,no_cores = 60,metabname,data,
                     out_dir = "/home/vigerst/MS-Thesis/data/candidate_selection/step_3/") {
  require(parallel)
  # Make cluster
  cl <- makeCluster(no_cores,type = "FORK")
  # Parallel models
  result_list <- parLapply(cl,mods,function(x){
    form <- as.formula(x)
    mod <- tryCatch(glm(form,data = data,family = "binomial"),
                    message = function(m) NULL,warning = function(m) NULL,
                    error = function(m) NULL)
    if (!is.null(mod)) {
      results <- as.data.frame(summary(mod)$coefficients)
      results$term <- rownames(results)
      results[nrow(results),"metab"] <- strsplit(x,"~")[[1]][2]
      results <- results[nrow(results),c("metab","Estimate","Pr(>|z|)")]
      colnames(results) <- c("metab","Value","p-value")
      return(results)
    } else {
      results <- as.data.frame(matrix(c(NA,NA,NA),nrow = 1))
      colnames(results) <- c("methyl","Value","p-value")
      return(results)
    }
  })
  df <- do.call(rbind,result_list)
  filename <- paste0(out_dir,metabname,"_unadj.csv")
  write.csv(df,file = filename,row.names = F)
  stopCluster(cl)
}

# gctof
model_list <- paste0("T1Dgroup~",names(gctof)[2:ncol(gctof)])
run_mods(model_list,metabname = "gctof",data = gctof)

# hilic
model_list <- paste0("T1Dgroup~",names(hilic)[2:ncol(hilic)])
run_mods(model_list,metabname = "hilic",data = hilic)

# lipid
model_list <- paste0("T1Dgroup~",names(lipid)[2:ncol(lipid)])
run_mods(model_list,metabname = "lipid",data = lipid)

# oxylipin
model_list <- paste0("T1Dgroup~",names(oxylipin)[2:ncol(oxylipin)])
run_mods(model_list,metabname = "oxylipin",data = oxylipin)

# vitd
model_list <- paste0("T1Dgroup~",names(vitd)[2:ncol(vitd)])
run_mods(model_list,metabname = "vitd",data = vitd)