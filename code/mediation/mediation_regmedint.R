library(regmedint)
library(parallel)
set.seed(1017)
# Load data
setwd("/home/vigerst/MS-Thesis/")
#setwd("~/Dropbox/School/MS Thesis")
load("./data/raw_data/probesFromPipeline.Rdata")
load("./data/raw_data/psv_sv_dataset.Rdata")
# List pairs
metabolites = read.csv("./data/metabolomics/liz_candidates.csv",na.strings = "")
metabolites = as.character(unlist(metabolites))
all = expand.grid(probesFromPipeline,metabolites[!is.na(metabolites)],stringsAsFactors = F)
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
  m = as.data.frame(summary(regmedint_obj)$summary_myreg)
  return(m)
}
# Iterate through all - parallel
cores = 12
cl = makeCluster(cores,type = "FORK")
methyl_psv_results = parApply(cl,all,1,function(r){
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
stopCluster(cl)
# Save
names(methyl_psv_results) = apply(all,1,paste,collapse = " & ")
save(methyl_psv_results,file = "./data/mediation/methyl_psv_liz_results.Rdata")
rm(methyl_psv_results)
# Same again for metab at PSV
# Iterate through all - parallel
cl = makeCluster(cores,type = "FORK")
metab_psv_results = apply(all,1,function(r){
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
stopCluster(cl)
# Save
names(metab_psv_results) = apply(all,1,paste,collapse = " & ")
save(metab_psv_results,file = "./data/mediation/metab_psv_liz_results.Rdata")