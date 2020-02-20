# Data import
methyl <- read.csv("/home/vigerst/MS-Thesis/data/step_1/sv/methyl_adj.csv",stringsAsFactors = F)

# Liz's candidates
candidates <- read.csv("/home/vigerst/MS-Thesis/data/metabolomics/liz_candidates.csv",
                       stringsAsFactors = F,na.strings = "")

# Model function
run_mods <- function(mods = model_list, data = temp,metabname,no_cores = 60,
                     out_dir = "/home/vigerst/MS-Thesis/candidate_selection/step_1/sv/") {
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
      results[nrow(results),"metab"] <- strsplit(x,"~")[[1]][2]
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
  filename <- paste0(out_dir,metabname,"_SV_adj_scaled.csv")
  write.csv(df,file = filename,row.names = F)
  stopCluster(cl)
}

# gctof
temp <- merge(gctof,methyl,by = "samplekey")
metab <- unique(candidates$gctof[!is.na(candidates$gctof)])
model_list <- paste0(rep(probesFromPipeline,each = length(metab)),"~",metab)

run_mods(model_list[1:100],metabname = "gctof")

# hilic
temp <- merge(hilic,methyl,by = "samplekey")
metab <- unique(candidates$hilic[!is.na(candidates$hilic)])
model_list <- paste0(rep(probesFromPipeline,each = length(metab)),"~",metab)

run_mods(model_list,metabname = "hilic")

# lipid
temp <- merge(lipid,methyl,by = "samplekey")
metab <- unique(candidates$lipid[!is.na(candidates$lipid)])
model_list <- paste0(rep(probesFromPipeline,each = length(metab)),"~",metab)

run_mods(model_list,metabname = "lipid")

# oxylipin
temp <- merge(oxylipin,methyl,by = "samplekey")
metab <- names(oxylipin)[2:ncol(oxylipin)]
model_list <- paste0(rep(probesFromPipeline,each = length(metab)),"~",metab)

run_mods(model_list,metabname = "oxylipin")

# vitd
temp <- merge(vitd,methyl,by = "samplekey")
metab <- names(vitd)[2:ncol(vitd)]
model_list <- paste0(rep(probesFromPipeline,each = length(metab)),"~",metab)

run_mods(model_list,metabname = "vitd")