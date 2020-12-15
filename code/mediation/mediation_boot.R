library(mediation)
library(parallel)
set.seed(1017)
# Load data
setwd("/Users/timvigers/GitHub/MS-Thesis")
load("./data/networks/all_data.Rdata")
load("./data/networks/pair_list.Rdata")
all_data$T1Dgroup = factor(all_data$T1Dgroup,)
# Mediation model lists
out_list = unique(paste0("T1Dgroup~",pairs$methyl,"+",pairs$metab))
# Cluster
cl = makeCluster(8,type = "FORK")
# Iterate through all
mediation = parLapply(cl,out_list, function(t){
  try({
    split = strsplit(t,"\\~|\\+")
    y_form = as.formula(paste0(t,"+clinage+SEX+dr34"))
    y_mod = glm(y_form,data = all_data,family = binomial("logit"))
    # methyl mediator
    m = split[[1]][2]
    x = split[[1]][3]
    m1_form = as.formula(paste0(m,"~",x,"+clinage+SEX+dr34"))
    m_mod = lm(m1_form,data = all_data)
    med1 = mediate(m_mod,y_mod,mediator = m,treat = x)
    r1 = c(x,m,med1$tau.coef,med1$tau.p,med1$d.avg,med1$d.avg.p,med1$z.avg,med1$z.avg.p,med1$n.avg,med1$n.avg.p)
    # metabolite mediator
    m = split[[1]][3]
    x = split[[1]][2]
    m2_form = as.formula(paste0(m,"~",x,"+clinage+SEX"))
    m_mod = lm(m2_form,data = all_data)
    med2 = mediate(m_mod,y_mod,mediator = m,treat = x)
    r2 = c(x,m,med2$tau.coef,med2$tau.p,med2$d.avg,med2$d.avg.p,med2$z.avg,med2$z.avg.p,med2$n.avg,med2$n.avg.p)
    # Results
    r = rbind(r1,r2)
    })
}) 
stopCluster(cl)
# Format
mediation = do.call(rbind,mediation)
rownames(mediation) = 1:nrow(mediation)
colnames(mediation) = c("X","M","Total Effect","Total Effect p",
                        "Average ACME","Average ACME p",
                        "Average ADE","Average ADE p",
                        "Average Prop. Med.","Average Prop. Med. p")
# Save results
mediation_quasi = mediation
save(mediation_quasi,file = "./data/mediation/mediation_quasi.Rdata")
