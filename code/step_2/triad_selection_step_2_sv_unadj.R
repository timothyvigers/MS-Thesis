# Data import
# Phenotype
pheno <- read.csv("/home/biostats_share/Norris/data/phenotype/ivyomicssample.csv",
                  stringsAsFactors = F,
                  na.strings = "")
pheno <- pheno[with(pheno,order(ID,DOVISIT)),]
pheno <- pheno[!is.na(pheno$T1Dgroup),]
pheno <- pheno[pheno$Visit_Type == "SV",]
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
methyl <- methyl[match(pheno$samplekey,methyl$samplekey),]
methyl$id <- factor(pheno$ID[match(methyl$samplekey,pheno$samplekey)])
methyl$T1Dgroup <- factor(pheno$T1Dgroup[match(methyl$samplekey,pheno$samplekey)])

# Model function
run_mods <- function(mods = model_list, data = methyl,no_cores = 60,
                     out_dir = "/home/vigerst/MS-Thesis/candidate_selection/step_2/") {
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
      results[nrow(results),"methyl"] <- strsplit(x,"~")[[1]][2]
      results <- results[nrow(results),c("methyl","Estimate","Pr(>|z|)")]
      colnames(results) <- c("methyl","Value","p-value")
      return(results)
    } else {
      results <- as.data.frame(matrix(c(NA,NA,NA),nrow = 1))
      colnames(results) <- c("methyl","Value","p-value")
      return(results)
    }
  })
  df <- do.call(rbind,result_list)
  filename <- paste0(out_dir,"methyl_unadj.csv")
  write.csv(df,file = filename,row.names = F)
  stopCluster(cl)
}

# models
model_list <- paste0("T1Dgroup~",names(methyl)[1:(ncol(methyl)-3)])

run_mods(model_list)
