# Parallel package cores
no_cores <- 60
# Data import
# Phenotype
pheno <- read.csv("/home/biostats_share/Norris/data/phenotype/ivyomicssample.csv",
                  stringsAsFactors = F,
                  na.strings = "")
pheno <- pheno[with(pheno,order(ID,DOVISIT)),]
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
methyl$id <- factor(pheno$ID[match(methyl$samplekey,pheno$samplekey)])
methyl$visit <- factor(pheno$Visit_Type[match(methyl$samplekey,pheno$samplekey)])
methyl$platform <- key$platform[match(methyl$samplekey,pheno$samplekey)]
# Metabolites
gctof <- read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.bc.csv")
hilic <- read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.bc.csv")
lipid <- read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.bc.csv")
oxylipin <- read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.bc.csv")
vitd <- read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.bc.csv")
# Liz's candidates
candidates <- read.csv("/home/vigerst/MS-Thesis/data/metabolomics/liz_candidates.csv",
                       stringsAsFactors = F,na.strings = "")

# Model function
run_mods <- function(mods = model_list, data = temp,metabname,
                     out_dir = "/home/vigerst/MS-Thesis/candidate_selection/step_1") {
  require(parallel)
  require(nlme)
  # Make cluster
  cl <- makeCluster(no_cores,type = "FORK")
  # Parallel models
  result_list <- parLapply(cl,mods,function(x){
    form <- as.formula(x)
    mod <- tryCatch(lme(form,random = ~1|id,data = data,na.action = na.omit),
                    message = function(m) NULL,warning = function(m) NULL,
                    error = function(m) NULL)
    if (!is.null(mod)) {
      results <- as.data.frame(summary(mod)$tTable)
      results$term <- rownames(results)
      results[nrow(results),"methyl"] <- strsplit(x,"~")[[1]][1]
      results[nrow(results),"metab"] <- strsplit(x,"\\+")[[1]][5]
      return(results[nrow(results),c("methyl","metab","Value","p-value")])
    } else {
      results <- as.data.frame(matrix(c(NA,NA,NA,NA),nrow = 1))
      colnames(results) <- c("methyl","metab","Value","p-value")
      return(results)
    }
  })
  df <- do.call(rbind,result_list)
  filename <- paste0(out_dir,metabname,"_longitudinal_adj.csv")
  write.csv(df,file = filename,row.names = F)
  stopCluster(cl)
}

# gctof
temp <- merge(gctof,methyl,by = "samplekey")
metab <- unique(candidates$gctof[!is.na(candidates$gctof)])
probes <- paste0(probesFromPipeline,"~sex+age+platform+visit+")
model_list <- paste0(rep(probes,each = length(metab)),metab)

run_mods(metabname = "gctof")

# hilic
temp <- merge(hilic,methyl,by = "samplekey")
metab <- unique(candidates$hilic[!is.na(candidates$hilic)])
model_list <- paste0(rep(probes,each = length(metab)),metab)

run_mods(metabname = "hilic")

# lipid
temp <- merge(lipid,methyl,by = "samplekey")
metab <- unique(candidates$lipid[!is.na(candidates$lipid)])
model_list <- paste0(rep(probes,each = length(metab)),metab)

run_mods(metabname = "lipid")

# oxylipin
temp <- merge(oxylipin,methyl,by = "samplekey")
metab <- names(oxylipin)[2:ncol(oxylipin)]
model_list <- paste0(rep(probes,each = length(metab)),metab)

run_mods(metabname = "oxylipin")

# vitd
temp <- merge(vitd,methyl,by = "samplekey")
metab <- names(vitd)[2:ncol(vitd)]
model_list <- paste0(rep(probes,each = length(metab)),metab)

run_mods(metabname = "vitd")
