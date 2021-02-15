library(mediation)
library(parallel)
set.seed(1017)
# Load data
setwd("/Users/timvigers/GitHub/MS-Thesis/")
load("./data/mediation/methyl_psv_candidates_p_05.Rdata")
load("./data/raw_data/psv_sv_dataset.Rdata")
# Outcome and adjustment variables
ia = factor(psv$IAgroup2)
SEX = factor(psv$SEX)
dr34 = factor(psv$dr34,labels = c("No","Yes"))
# Change in age from PSV to SV
age = sv$clinage - psv$clinage
# Cluster
cl = makeCluster(8,type = "FORK")
# Iterate through all
mediation_mods = parApply(cl,methyl_psv_candidates,1,function(r){
  methyl = psv[,r[1]]
  metab = sv[,r[2]]
  # Mediation models
  m = lm(metab~methyl+age+SEX+dr34)
  c = glm(ia ~ methyl+metab+age+SEX+dr34,family = binomial("logit"))
  med = mediate(m,c,treat="methyl",mediator="metab",boot = T,sims = 100000)
}) 
stopCluster(cl)
# Save
save(mediation_mods,file = "./data/mediation/longitudinal_mediation_p_05_100k_sims.Rdata")
