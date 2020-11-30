library(mediation)
library(parallel)
# Load data
setwd("/Users/timvigers/GitHub/MS-Thesis")
load("./data/networks/all_data_methyl_scaled.Rdata")
load("./data/networks/cits.Rdata")
set.seed(1017)
# Mediation model lists
out_list = unique(paste0("factor(T1Dgroup)~",cits$methyl,"+",cits$metab))
out_list = out_list[1:2]
# Cluster
cl = makeCluster(8,type = "FORK")
# Iterate through all
cits_mediation = parLapply(cl,out_list, function(t){
  split = strsplit(t,"\\~|\\+")
  y_mod = glm(as.formula(t),data = all_data_methyl_scaled,family = "binomial")
  # methyl mediator
  m = split[[1]][2]
  x = split[[1]][3]
  m1_form = as.formula(paste0(m,"~",x))
  m_mod = lm(m1_form,data = all_data_methyl_scaled)
  med1 = mediate(m_mod,y_mod,mediator = m,treat = x)
  r1 = c(x,m,med1$tau.coef,med1$tau.p,med1$d.avg,med1$d.avg.p,med1$z.avg,med1$z.avg.p,med1$n.avg,med1$n.avg.p)
  # metabolite mediator
  m = split[[1]][3]
  x = split[[1]][2]
  m2_form = as.formula(paste0(m,"~",x))
  m_mod = lm(m2_form,data = all_data_methyl_scaled)
  med2 = mediate(m_mod,y_mod,mediator = m,treat = x)
  r2 = c(x,m,med2$tau.coef,med2$tau.p,med2$d.avg,med2$d.avg.p,med2$z.avg,med2$z.avg.p,med2$n.avg,med2$n.avg.p)
  # Results
  r = rbind(r1,r2)
}) 
cits_mediation = do.call(rbind,cits_mediation)
rownames(cits_mediation) = 1:nrow(cits_mediation)
colnames(cits_mediation) = c("X","M","Total Effect","Total Effect p",
                             "Average ACME","Average ACME p",
                             "Average ADE","Average ADE p",
                             "Average Prop. Med.","Average Prop. Med. p")
# Sensitivity analysis for all
cits_mediation_sensitivity = parLapply(cl,out_list, function(t){
  split = strsplit(t,"\\~|\\+")
  y = glm(as.formula(t),data = all_data_methyl_scaled,family = binomial("probit"))
  # methyl mediator
  m1_form = as.formula(paste0(split[[1]][2],"~",split[[1]][3]))
  m = lm(m1_form,data = all_data_methyl_scaled)
  med1 = mediate(m,y,mediator = split[[1]][2],treat = split[[1]][3])
  sens1 = medsens(med1)
  # metabolite mediator
  m2_form = as.formula(paste0(split[[1]][3],"~",split[[1]][2]))
  m = lm(m2_form,data = all_data_methyl_scaled)
  med2 = mediate(m,y,mediator = split[[1]][3],treat = split[[1]][2])
  sens2 = medsens(med2)
  # Results
  list(sens1,sens2)
}) 
save(cits_mediation,file = "./data/mediation/cits_mediation.Rdata")
save(cits_mediation_sensitivity,file = "./data/mediation/cits_mediation_sensitivity.Rdata")
stopCluster(cl)
