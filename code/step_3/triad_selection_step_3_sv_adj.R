# Data import
# Phenotype
pheno <- read.csv("/home/biostats_share/Norris/data/phenotype/ivyomicssample.csv",
                  stringsAsFactors = F,
                  na.strings = "")
pheno$T1Dgroup <- as.factor(pheno$T1Dgroup)
pheno$sex <- as.factor(pheno$SEX)
pheno$age <- as.numeric(pheno$clinage)
pheno <- pheno[with(pheno,order(ID,DOVISIT)),]
pheno <- pheno[!is.na(pheno$T1Dgroup),]
pheno <- pheno[pheno$Visit_Type == "SV",]
# Metabolites
gctof <- read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.bc.csv")
gctof <- merge(pheno,gctof,by = "samplekey")
hilic <- read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.bc.csv")
hilic <- merge(pheno,hilic,by = "samplekey")
lipid <- read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.bc.csv")
lipid <- merge(pheno,lipid,by = "samplekey")
oxylipin <- read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.bc.csv")
oxylipin <- merge(pheno,oxylipin,by = "samplekey")
vitd <- read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.bc.csv")
vitd <- merge(pheno,vitd,by = "samplekey")

# Model function
run_mods <- function(mods = model_list,no_cores = 60,metabname,data,
                     out_dir = "/home/vigerst/MS-Thesis/candidate_selection/step_3/") {
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
  filename <- paste0(out_dir,metabname,"_adj.csv")
  write.csv(df,file = filename,row.names = F)
  stopCluster(cl)
}

# gctof
model_list <- paste0("T1Dgroup~sex+age+",names(gctof)[20:ncol(gctof)])
run_mods(model_list[1:50],metabname = "gctof",data = gctof)

# hilic
model_list <- paste0("T1Dgroup~sex+age+",names(hilic)[20:ncol(hilic)])
run_mods(model_list[1:50],metabname = "hilic",data = hilic)

# lipid
model_list <- paste0("T1Dgroup~sex+age+",names(lipid)[20:ncol(lipid)])
run_mods(model_list[1:50],metabname = "lipid",data = lipid)

# oxylipin
model_list <- paste0("T1Dgroup~sex+age+",names(oxylipin)[20:ncol(oxylipin)])
run_mods(model_list,metabname = "oxylipin",data = oxylipin)

# vitd
model_list <- paste0("T1Dgroup~sex+age+",names(vitd)[20:ncol(vitd)])
run_mods(model_list,metabname = "vitd",data = vitd)