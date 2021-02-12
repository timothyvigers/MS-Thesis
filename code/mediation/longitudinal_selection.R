library(parallel)
setwd("/home/vigerst/MS-Thesis/data/raw_data")
load("./psv_sv_dataset.Rdata")
load("./probesFromPipeline.Rdata")
ia = factor(psv$IAgroup2)
# Parallel
n_cores = 16
cl = makeCluster(n_cores,type = "FORK")
# Methylation at PSV, metabolite at SV
methyl_psv_candidates = parLapply(cl,probesFromPipeline, function(p){
  candidates = lapply(metabolites, function(m){
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
cl = makeCluster(n_cores,type = "FORK")
# Methylation at SV, metabolite at PSV
metab_psv_candidates = parLapply(cl,probesFromPipeline, function(p){
  candidates = lapply(metabolites, function(m){
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
