library(parallel)
setwd("C:/Users/Tim Vigers/Dropbox/School/MS Thesis/data/raw_data")
load("./psv_sv_dataset.Rdata")
load("./probesFromPipeline.Rdata")
metab_candidates = read.csv("./liz_candidates.csv",stringsAsFactors = F,na.strings = "")
metab_candidates = unlist(metab_candidates)
metab_candidates = unique(metab_candidates[!is.na(metab_candidates)])
# Parallel
n_cores = 8
cl = makeCluster(n_cores)
clusterExport(cl,list('metab_candidates','psv','sv','probesFromPipeline'))
# Methylation at PSV, metabolite at SV
methyl_psv_candidates = parLapply(cl,probesFromPipeline, function(p){
  candidates = lapply(metab_candidates, function(m){
    ia = factor(psv$IAgroup2)
    meth = psv[,p]
    meta = sv[,m]
    a = summary(lm(meta ~ meth))$coefficients[2,4]
    b = summary(glm(ia ~ meta,family = binomial("logit")))$coefficients[2,4]
    c = summary(glm(ia ~ meth,family = binomial("logit")))$coefficients[2,4]
    if(a < 0.05 & b < 0.05 & c < 0.05){
      return(c(p,m))
    } else {return(c(NA,NA))}
  })
  candidates = do.call(rbind,candidates)
  candidates[complete.cases(candidates),]
})
stopCluster(cl)
methyl_psv_candidates = do.call(rbind,methyl_psv_candidates)
save(methyl_psv_candidates,file = "./methyl_psv_candidates_p_05.Rdata")
rm("methyl_psv_candidates")
# Restart cluster
cl = makeCluster(n_cores)
clusterExport(cl,list('metab_candidates','psv','sv','probesFromPipeline'))
# Methylation at SV, metabolite at PSV
metab_psv_candidates = parLapply(cl,probesFromPipeline, function(p){
  candidates = lapply(metab_candidates, function(m){
    ia = factor(psv$IAgroup2)
    meth = sv[,p]
    meta = psv[,m]
    a = summary(lm(meta ~ meth))$coefficients[2,4]
    b = summary(glm(ia ~ meta,family = binomial("logit")))$coefficients[2,4]
    c = summary(glm(ia ~ meth,family = binomial("logit")))$coefficients[2,4]
    if(a < 0.05 & b < 0.05 & c < 0.05){
      return(c(p,m))
    } else {return(c(NA,NA))}
  })
  candidates = do.call(rbind,candidates)
  candidates[complete.cases(candidates),]
})
stopCluster(cl)
metab_psv_candidates = do.call(rbind,metab_psv_candidates)
save(metab_psv_candidates,file = "./metab_psv_candidates_p_05.Rdata")