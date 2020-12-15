library(mediation)
library(parallel)
set.seed(1017)
# Load data
setwd("/Users/timvigers/GitHub/MS-Thesis")
load("./data/networks/all_data.Rdata")
load("./data/networks/cits.Rdata")
all_data$T1Dgroup = factor(all_data$T1Dgroup)
# Mediation model list
out_list = unique(paste0("T1Dgroup~",cits$methyl,"+",cits$metab))
# Cluster
cl = makeCluster(4,type = "FORK")
# Iterate through all
cits_mediation = parLapply(cl,out_list, function(t){
  try({
    # methyl mediator
    split = strsplit(t,"\\~|\\+")
    df = all_data[complete.cases(all_data[,c(split[[1]][2],split[[1]][3])]),]
    t1d = c(df$T1Dgroup) - 1
    m = c(df[,split[[1]][2]])
    x = c(df[,split[[1]][3]])
    clinage = c(df$clinage)
    SEX = c(df$SEX)
    dr34 = c(df$dr34)
    # Y
    y_mod = glm(t1d ~ m+x+clinage+SEX+dr34,family = binomial("logit"))
    # M
    m_mod = lm(m ~ x+clinage+SEX+dr34)
    # Mediate
    med1 = mediate(m_mod,y_mod,mediator = "m",treat = "x",boot = T)
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
cits_mediation = do.call(rbind,cits_mediation)
rownames(cits_mediation) = 1:nrow(cits_mediation)
colnames(cits_mediation) = c("X","M","Total Effect","Total Effect p",
                             "Average ACME","Average ACME p",
                             "Average ADE","Average ADE p",
                             "Average Prop. Med.","Average Prop. Med. p")
# Save
cits_mediation_boot = cits_mediation
save(cits_mediation_boot,file = "./data/mediation/cits_mediation_boot.Rdata")
