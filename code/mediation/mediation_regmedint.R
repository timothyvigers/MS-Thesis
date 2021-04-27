library(regmedint)
library(parallel)
set.seed(1017)
# Load data
setwd("/home/vigerst/MS-Thesis/")
#setwd("~/Dropbox/School/MS Thesis")
load("./data/raw_data/probesFromPipeline.Rdata")
load("./data/raw_data/psv_sv_dataset.Rdata")
# List pairs
all = expand.grid(probesFromPipeline,metabolites,stringsAsFactors = F)
# Outcome and adjustment variables
ia = as.numeric(factor(psv$IAgroup2)) - 1
SEX = as.numeric(factor(psv$SEX)) - 1
dr34 = psv$dr34
age_delta = as.numeric(sv$clinage) - as.numeric(psv$clinage)
age = as.numeric(psv$clinage)
covariates = as.data.frame(cbind(ia,SEX,dr34,age,age_delta))
# Mediation function
regmed = function(d){
  regmedint_obj = 
    regmedint(data = d,
              ## Variables
              yvar = "ia",
              avar = "methyl",
              mvar = "metab",
              cvar = c("SEX","dr34","age","age_delta"),
              ## Values at which effects are evaluated
              a0 = 0,
              a1 = 1,
              m_cde = 1,
              c_cond = c(1,1,1,1),
              ## Model types
              mreg = "linear",
              yreg = "logistic",
              ## Additional specification
              interaction = T,
              casecontrol = T)
  m = summary(regmedint_obj)$summary_myreg
  m["tnie","p"]
}
# Iterate through all - parallel
cores = 8
cl = makeCluster(cores,type = "FORK")
methyl_psv_pvalues = parApply(cl,all,1,function(r){
  methyl = as.character(r["Var1"])
  methyl = as.numeric(scale(psv[,methyl]))
  metab = as.character(r["Var2"])
  metab = as.numeric(scale(sv[,metab]))
  # Dataframe 
  df = as.data.frame(cbind(methyl,metab,covariates))
  df = df[complete.cases(df),]
  # Mediation
  med = regmed(df)
  return(med)
})
methyl_psv_pvalues = as.data.frame(methyl_psv_pvalues)
rownames(methyl_psv_pvalues) = apply(all,1,paste,collapse = " & ")
# Save
save(methyl_psv_pvalues,file = "./data/mediation/methyl_psv_all_pvalues.Rdata")
stopCluster(cl)
# Same again for metab at PSV
# Iterate through all - parallel
cl = makeCluster(cores,type = "FORK")
metab_psv_pvalues = apply(all,1,function(r){
  metab = as.character(r["Var2"])
  metab = as.numeric(scale(psv[,metab]))
  methyl = as.character(r["Var1"])
  methyl = as.numeric(scale(sv[,methyl]))
  # Dataframe 
  df = as.data.frame(cbind(methyl,metab,age,covariates))
  df = df[complete.cases(df),]
  # Mediation
  med = regmed(df)
  return(med)
})
metab_psv_pvalues = as.data.frame(metab_psv_pvalues)
rownames(metab_psv_pvalues) = apply(all,1,paste,collapse = " & ")
# Save
save(metab_psv_pvalues,file = "./data/mediation/metab_psv_all_pvalues.Rdata")
stopCluster(cl)