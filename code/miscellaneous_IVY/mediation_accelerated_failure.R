library(regmedint)
library(boot)
set.seed(1017)
# Load data
#setwd("/home/vigerst/MS-Thesis/")
setwd("~/Documents/School/MS Thesis")
df = read.csv("./data/raw_data/Progression_Sugar_Mediation_Dataset_Expanded.csv",
              na.strings = "")
# Outcome and adjustment variables. Per Patrick:
# "When we ran the mediation models, we included the following variables as covariates/confounders: dr34,  NHW  Age_SV,  Female_YN,  calor (calories), FFQ_YN (type of questionnaire used). 
# Most importantly,  T1D_Time and T1D_Event are the variables used for the progression analysis."
dr34 = df$dr34
nhw = df$NHW
age = df$Age_SV
sex = df$Female_YN
cal = df$calor
ffq = factor(df$FFQ_type)
t1d_time = df$T1D_time
t1d_event = df$T1D_Event
# Bootstrap function
regmed_boot = function(d,i){
  regmedint_obj = 
    regmedint(data = d[i,],
              ## Variables
              yvar = "ia",
              avar = "methyl",
              mvar = "metab",
              cvar = c("dr34","nhw","age","sex","cal","ffq"),
              eventvar = "t1d_event",
              ## Values at which effects are evaluated
              a0 = 0,
              a1 = 1,
              m_cde = 1,
              c_cond = c(1,1,1,1),
              ## Model types
              mreg = "linear",
              yreg = "survAFT_exp",
              ## Additional specification
              interaction = T,
              casecontrol = T,
              na_omit = T)
  summary(regmedint_obj)$summary_myreg[,1]
}
# Bootstrap options
boot_cores = 12
boots = 10000
# Iterate through all
methyl_psv_results = apply(methyl_psv_candidates,1,function(r){
  methyl = as.numeric(scale(psv[,r[1]]))
  metab = as.numeric(scale(sv[,r[2]]))
  # Dataframe 
  df = as.data.frame(cbind(methyl,metab,age,covariates))
  # Bootstrap
  b = boot(data = df, statistic = regmed_boot, R = boots,parallel = "multicore",
           ncpus = boot_cores)
  return(b)
})
names(methyl_psv_results) = apply(methyl_psv_candidates,1,paste,collapse = " & ")
# Save
save(methyl_psv_results,file = "./data/mediation/methyl_psv_both_scaled_results.Rdata")
# Same again for metab at PSV
# Iterate through all
metab_psv_results = apply(metab_psv_candidates,1,function(r){
  metab = as.numeric(scale(psv[,r[2]]))
  methyl = as.numeric(scale(sv[,r[1]]))
  # Dataframe 
  df = as.data.frame(cbind(methyl,metab,age,covariates))
  # Mediation
  b = boot(data = df, statistic = regmed_boot, R = boots,parallel = "multicore",
           ncpus = boot_cores)
  return(b)
})
names(metab_psv_results) = apply(metab_psv_candidates,1,paste,collapse = " & ")
# Save
save(metab_psv_results,file = "./data/mediation/metab_psv_both_scaled_results.Rdata")