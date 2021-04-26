library(regmedint)
library(boot)
library(broom)
set.seed(1017)
# Load data
#setwd("/home/vigerst/MS-Thesis/")
setwd("~/Dropbox/School/MS Thesis")
load("./data/raw_data/probesFromPipeline.Rdata")
load("./data/raw_data/psv_sv_dataset.Rdata")
# Outcome and adjustment variables
ia = as.numeric(factor(psv$IAgroup2)) - 1
SEX = as.numeric(factor(psv$SEX)) - 1
dr34 = psv$dr34
age_delta = as.numeric(sv$clinage) - as.numeric(psv$clinage)
age = as.numeric(psv$clinage)
covariates = as.data.frame(cbind(ia,SEX,dr34,age,age_delta))
# Bootstrap function
regmed_boot = function(d,i){
  regmedint_obj = 
    regmedint(data = d[i,],
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
              casecontrol = T,
              na_omit = T)
  summary(regmedint_obj)$summary_myreg[,1]
}
# Bootstrap options
boot_cores = 16
boots = 10000
conf.m = "perc"
# Iterate through all
all = expand.grid(probesFromPipeline,metabolites,stringsAsFactors = F)
methyl_psv_results = apply(all,1,function(r){
  methyl = as.character(r["Var1"])
  methyl = as.numeric(scale(psv[,methyl]))
  metab = as.character(r["Var2"])
  metab = as.numeric(scale(sv[,metab]))
  # Dataframe 
  df = as.data.frame(cbind(methyl,metab,covariates))
  # Bootstrap
  b = boot(data = df, statistic = regmed_boot, R = boots,parallel = "multicore",
           ncpus = boot_cores)
  b = tidy(b,conf.int = T,conf.method = conf.m)
  return(b)
})
names(methyl_psv_results) = apply(all,1,paste,collapse = " & ")
# Save
save(methyl_psv_results,file = "./data/mediation/methyl_psv_all_results.Rdata")
# Same again for metab at PSV
# Iterate through all
metab_psv_results = apply(all,1,function(r){
  metab = as.character(r["Var2"])
  metab = as.numeric(scale(psv[,metab]))
  methyl = as.character(r["Var1"])
  methyl = as.numeric(scale(sv[,methyl]))
  # Dataframe 
  df = as.data.frame(cbind(methyl,metab,age,covariates))
  # Mediation
  b = boot(data = df, statistic = regmed_boot, R = boots,parallel = "multicore",
           ncpus = boot_cores)
  b = tidy(b,conf.int = T,conf.method = conf.m)
  return(b)
})
names(metab_psv_results) = apply(all,1,paste,collapse = " & ")
# Save
save(metab_psv_results,file = "./data/mediation/metab_psv_all_results.Rdata")