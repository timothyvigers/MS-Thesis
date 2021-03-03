import multiprocessing as mp
import pandas as pd
import statsmodels.api as sm
# Home directory
wd = "/Users/timvigers/Dropbox/School/MS Thesis/data/"
# Import data
psv = pd.read_csv(wd+"raw_data/psv.csv")
psv
#sv = pd.read_csv(wd+"raw_data/sv.csv")
# Get outcome
ia = pd.read_csv(wd+"raw_data/psv.csv",usecols=["IAgroup2"])
ia = ia.replace({'case': 1,'control':0})
# Get lists of probes and metabolites to test
probes = pd.read_csv(wd+"raw_data/probes.csv",header=None)
probes = probes.astype(str).values.tolist()
metabolites = pd.read_csv(wd+"raw_data/metabolites.csv",header=None)
metabolites = metabolites.astype(str).values.tolist()
# For each probe, test the association with IA status
probe_candidates = [] # To store results
for i in range(0,len(probes)+1):
    c = probes[i][0]
    logit_model=sm.Logit(ia,sm.add_constant(psv[[c]]))
    result=logit_model.fit(disp=0)
    if result.pvalues[[1]][0] < 0.5:
        probe_candidates.append(c)
# setwd("/home/vigerst/MS-Thesis/data")
# #setwd("/Users/timvigers/Dropbox/School/MS Thesis/data")
# load("./raw_data/psv_sv_dataset.Rdata")
# load("./raw_data/probesFromPipeline.Rdata")
# ia = factor(psv$IAgroup2)
# # Parallel
# n_cores = 8
# # Methylation at PSV, metabolite at SV
# # Probes associated with outcome
# cl = makeCluster(n_cores,type = "FORK")
# probe_candidates = parLapply(cl,probesFromPipeline, function(p){
#   meth = psv[,p]
#   c = summary(glm(ia ~ meth,family = binomial("logit")))$coefficients[2,4]
#   if(c < 0.05){return(p)}else{return(NA)}
# })
# stopCluster(cl)
# probe_candidates = unlist(probe_candidates[!is.na(probe_candidates)])
# # Metabolites associated with outcome
# cl = makeCluster(n_cores,type = "FORK")
# metab_candidates = parLapply(cl,metabolites, function(m){
#   meta = sv[,m]
#   c = summary(glm(ia ~ meta,family = binomial("logit")))$coefficients[2,4]
#   if(c < 0.05){return(m)}else{return(NA)}
# })
# stopCluster(cl)
# metab_candidates = unlist(metab_candidates[!is.na(metab_candidates)])
# # Check if these are related
# cl = makeCluster(n_cores,type = "FORK")
# methyl_psv_candidates = parLapply(cl,probe_candidates, function(p){
#   candidates = lapply(metab_candidates, function(m){
#     meth = psv[,p]
#     meta = sv[,m]
#     a = summary(lm(meta ~ meth))$coefficients[2,4]
#     if(a < 0.05){return(c(p,m))}else{return(c(NA,NA))}
#   })
#   candidates = do.call(rbind,candidates)
#   candidates[complete.cases(candidates),]
# })
# stopCluster(cl)
# methyl_psv_candidates = do.call(rbind,methyl_psv_candidates)
# save(methyl_psv_candidates,file = "./mediation/methyl_psv_candidates_p_05.Rdata")
# rm("methyl_psv_candidates")
# # Methylation at SV, metabolite at PSV
# # Probes associated with outcome
# cl = makeCluster(n_cores,type = "FORK")
# probe_candidates = parLapply(cl,probesFromPipeline, function(p){
#   meth = sv[,p]
#   c = summary(glm(ia ~ meth,family = binomial("logit")))$coefficients[2,4]
#   if(c < 0.05){return(p)}else{return(NA)}
# })
# stopCluster(cl)
# probe_candidates = unlist(probe_candidates[!is.na(probe_candidates)])
# # Metabolites associated with outcome
# cl = makeCluster(n_cores,type = "FORK")
# metab_candidates = parLapply(cl,metabolites, function(m){
#   meta = psv[,m]
#   c = summary(glm(ia ~ meta,family = binomial("logit")))$coefficients[2,4]
#   if(c < 0.05){return(m)}else{return(NA)}
# })
# stopCluster(cl)
# metab_candidates = unlist(metab_candidates[!is.na(metab_candidates)])
# # Check if these are related
# cl = makeCluster(n_cores,type = "FORK")
# metab_psv_candidates = parLapply(cl,probe_candidates, function(p){
#   candidates = lapply(metab_candidates, function(m){
#     meth = sv[,p]
#     meta = psv[,m]
#     a = summary(lm(meta ~ meth))$coefficients[2,4]
#     if(a < 0.05){return(c(p,m))}else{return(c(NA,NA))}
#   })
#   candidates = do.call(rbind,candidates)
#   candidates[complete.cases(candidates),]
# })
# stopCluster(cl)
# metab_psv_candidates = do.call(rbind,metab_psv_candidates)
# save(metab_psv_candidates,file = "./mediation/metab_psv_candidates_p_05.Rdata")
# rm("metab_psv_candidates")
