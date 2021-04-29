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
# Bootstrap options
boot_cores = 16
boots = 10000
# Bootstrap function
regmed_boot = function(d,i,probe,int = F,cond,surv_type){
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
              c_cond = cond,
              ## Model types
              mreg = "linear",
              yreg = surv_type,
              ## Additional specification
              na_omit = T,
              interaction = int)
  summary(regmedint_obj)$summary_myreg[,1]
}
# Get probes from column names
probes = colnames(df)[grep("^cg\\d",colnames(df))]
# Exp survival, 1 for conditional
results = lapply(probes,function(p){
  # Bootstrap
  b = suppressMessages(boot(data = df, statistic = regmed_boot, R = boots, parallel = "multicore",
           ncpus = boot_cores, probe = p,
           cond = c(1,1,1,1,1,1),surv_type = "survAFT_exp"))
  return(b)
})
names(results) = probes
# Save
save(results,file = "./results/aft_results_cond1_exp.RData")
# Weibull survival, 1 for conditional
results = lapply(probes,function(p){
  # Bootstrap
  b = suppressMessages(boot(data = df, statistic = regmed_boot, R = boots, parallel = "multicore",
                            ncpus = boot_cores, probe = p,
                            cond = c(1,1,1,1,1,1),surv_type = "survAFT_weibull"))
  return(b)
})
names(results) = probes
# Save
save(results,file = "./results/aft_results_cond1_weibull.RData")
# Exp survival, means for conditional
results = lapply(probes,function(p){
  # Bootstrap
  b = suppressMessages(boot(data = df, statistic = regmed_boot, R = boots, parallel = "multicore",
                            ncpus = boot_cores, probe = p,
                            cond = c(0.5,0.5,mean(df$Age_SV,na.rm=T),0.5,mean(df$calor,na.rm=T),0.5),surv_type = "survAFT_exp"))
  return(b)
})
names(results) = probes
# Save
save(results,file = "./results/aft_results_condmean_exp.RData")
# Weibull survival, 1 for conditional
results = lapply(probes,function(p){
  # Bootstrap
  b = suppressMessages(boot(data = df, statistic = regmed_boot, R = boots, parallel = "multicore",
                            ncpus = boot_cores, probe = p,
                            cond = c(0.5,0.5,mean(df$Age_SV,na.rm=T),0.5,mean(df$calor,na.rm=T),0.5),surv_type = "survAFT_weibull"))
  return(b)
})
names(results) = probes
# Save
save(results,file = "./results/aft_results_condmean_weibull.RData")
