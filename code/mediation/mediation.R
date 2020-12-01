library(mediation)
library(parallel)
set.seed(1017)
# Load data
setwd("/home/vigerst/MS-Thesis")
load("./data/networks/all_data_methyl_scaled.Rdata")
load("./data/networks/pair_list.Rdata")
# Mediation model lists
out_list = unique(paste0("factor(T1Dgroup)~",pairs$methyl,"+",pairs$metab))
# Cluster
cl = makeCluster(8,type = "FORK")
# Iterate through all
mediation = parLapply(cl,out_list, function(t){
  try({split = strsplit(t,"\\~|\\+")
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
  r = rbind(r1,r2)})
}) 
mediation = do.call(rbind,mediation)
rownames(mediation) = 1:nrow(mediation)
colnames(mediation) = c("X","M","Total Effect","Total Effect p",
                        "Average ACME","Average ACME p",
                        "Average ADE","Average ADE p",
                        "Average Prop. Med.","Average Prop. Med. p")
# Sensitivity analysis for all
mediation_sensitivity = parLapply(cl,out_list, function(t){
  try({split = strsplit(t,"\\~|\\+")
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
  list(sens1,sens2)})
}) 
save(mediation,file = "./data/mediation/mediation.Rdata")
save(mediation_sensitivity,file = "./data/mediation/mediation_sensitivity.Rdata")
stopCluster(cl)
