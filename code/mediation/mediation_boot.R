library(mediation)
library(parallel)
set.seed(1017)
# Load data
setwd("/Users/timvigers/GitHub/MS-Thesis/")
load("./data/raw_data/methyl_psv_candidates.Rdata")
load("./data/raw_data/psv_sv_dataset.Rdata")
# Outcome and adjustment variables
ia = factor(psv$IAgroup2)
clinage = psv$clinage
SEX = factor(psv$SEX)
dr34 = factor(psv$dr34,labels = c("No","Yes"))
# Cluster
cl = makeCluster(8,type = "FORK")
# Iterate through all
mediation_mods = parApply(cl,methyl_psv_candidates,1,function(r){
  methyl = psv[,r[1]]
  metab = sv[,r[2]]
  # Mediation models
  m = lm(metab~methyl+clinage+SEX+dr34)
  c = glm(ia ~ methyl+metab+clinage+SEX+dr34,family = binomial("logit"))
  med = mediate(m,c,treat="methyl",mediator="metab",boot = T,sims = 100000,long = F)
}) 
stopCluster(cl)
# Save
save(mediation_mods,file = "./data/mediation/longitudinal_mediation_mods_100k.Rdata")
