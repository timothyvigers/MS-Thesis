setwd("C:/Users/Tim Vigers/Documents/GitHub/MS-Thesis/data/raw_data")
load("./psv_sv_dataset.Rdata")
load("./probesFromPipeline.Rdata")
metab_candidates = read.csv("./liz_candidates.csv",stringsAsFactors = F,na.strings = "")
metab_candidates = unlist(metab_candidates)
metab_candidates = unique(metab_candidates[!is.na(metab_candidates)])
# Model functions
methyl_psv_mods = function(df,no_cores = 6) {
  require(parallel)
  # Make cluster
  # cl = makeCluster(no_cores,type = "FORK")
  # Models
  candidates = lapply(probesFromPipeline, function(p){
    
    t = df[]
    methyl = df[,p]
  })
  
  
  
}
