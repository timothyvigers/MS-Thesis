library(regmedint)
library(boot)
set.seed(1017)
# Load data
setwd("/home/vigerst/EWAS/")
df = read.csv("./data/Progression_Sugar_Mediation_Dataset_Expanded.csv",
              na.strings = "")
# Outcome and adjustment variables. Per Patrick:
# "When we ran the mediation models, we included the following variables as covariates/confounders: dr34,  NHW  Age_SV,  Female_YN,  calor (calories), FFQ_YN (type of questionnaire used). 
# Most importantly,  T1D_Time and T1D_Event are the variables used for the progression analysis."
df$FFQ_type = as.numeric(as.factor(df$FFQ_type))-1
# Bootstrap function
regmed_boot = function(d,i,probe){
  regmedint_obj = 
    regmedint(data = d[i,],
              ## Variables
              yvar = "T1D_time",
              avar = "Sugar.Z",
              mvar = probe,
              cvar = c("dr34","NHW","Age_SV","Female_YN","calor","FFQ_type"),
              eventvar = "T1D_Event",
              ## Values at which effects are evaluated
              a0 = 0,
              a1 = 1,
              m_cde = 1,
              c_cond = c(1,1,1,1,1,1),
              ## Model types
              mreg = "linear",
              yreg = "survAFT_exp",
              ## Additional specification
              na_omit = T)
  summary(regmedint_obj)$summary_myreg[,1]
}
# Bootstrap options
boot_cores = 16
boots = 10000
# Iterate through all
probes = colnames(df)[grep("^cg\\d",colnames(df))]
results = lapply(probes,function(p){
  # Bootstrap
  b = suppressMessages(boot(data = df, statistic = regmed_boot, R = boots, parallel = "multicore",
           ncpus = boot_cores, probe = p))
  return(b)
})
names(results) = probes
# Save
results_fruc_adjusted = results
save(results_fruc_adjusted,file = "./results/aft_results_fruc_adjusted.RData")
