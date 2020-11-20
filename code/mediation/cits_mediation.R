library(mediation)
# Load data
setwd("/Users/timvigers/GitHub/MS-Thesis")
load("./data/networks/all_data_methyl_scaled.Rdata")
load("./data/networks/cits.Rdata")
set.seed(1017)
# Mediation model lists
out_list = unique(paste0("factor(T1Dgroup)~",cits$methyl,"+",cits$metab))
# Iterate through all
cits_mediation = lapply(out_list[1], function(t){
  split = strsplit(t,"\\~|\\+")
  y1 = glm(as.formula(t),data = all_data_methyl_scaled,family = binomial("probit"))
  # methyl mediator
  m1_form = as.formula(paste0(split[[1]][2],"~",split[[1]][3]))
  m = lm(m1_form,data = all_data_methyl_scaled)
  med1 = mediate(m,y,mediator = split[[1]][2],treat = split[[1]][3])
  # metabolite mediator
  m2_form = as.formula(paste0(split[[1]][3],"~",split[[1]][2]))
  m = lm(m2_form,data = all_data_methyl_scaled)
  med2 = mediate(m,y,mediator = split[[1]][3],treat = split[[1]][2])
  # Results - what to report here?
  list(summary(med1),summary(med2))
}) 
save(cits_mediation,file = "./data/mediation/cits_mediation.Rdata")