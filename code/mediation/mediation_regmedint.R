library(regmedint)
library(boot)
set.seed(1017)
# Load data
#setwd("/home/vigerst/MS-Thesis/")
setwd("/Users/timvigers/Dropbox/School/MS Thesis")
load("./data/raw_data/psv_sv_dataset.Rdata")
load("./data/mediation/methyl_psv_candidates_p_01.Rdata")
load("./data/mediation/metab_psv_candidates_p_01.Rdata")
# Outcome and adjustment variables
ia = psv$IAgroup2
SEX = psv$SEX
dr34 = psv$dr34
age_delta = psv$clinage - sv$clinage
age = psv$clinage
covariates = as.data.frame(cbind(ia,SEX,dr34,age,age_delta))
covariates = as.data.frame(lapply(covariates,function(x){as.numeric(as.factor(x))-1}))
# Iterate through all
boot_cores = 4
boots = 10000
methyl_psv_results = apply(methyl_psv_candidates,1,function(r){
  methyl = psv[,r[1]]
  metab = sv[,r[2]]
  # Dataframe 
  df = as.data.frame(cbind(methyl,metab,age,covariates))
  # Mediation
  # Bootstrap
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
                casecontrol = T)
    coef(regmedint_obj)
  }
  b = boot(data = df, statistic = regmed_boot, R = boots,parallel = "multicore",
           ncpus = boot_cores)
  bci = boot.ci(b,type = "perc",index = 1)
})
names(methyl_psv_results) = apply(methyl_psv_candidates,1,paste,collapse = " & ")
# Save
save(methyl_psv_results,file = "./data/mediation/methyl_psv_results.Rdata")
# Same again for metab at PSV
cl = makeCluster(n_cores,type = "FORK")
# Iterate through all
metab_psv_results = apply(metab_psv_candidates,1,function(r){
  methyl = sv[,r[1]]
  metab = psv[,r[2]]
  # Dataframe 
  df = as.data.frame(cbind(methyl,metab,age,covariates))
  # Mediation
  regmedint_obj <- 
    regmedint(data = df,
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
  regmedint_obj
})
stopCluster(cl)
names(metab_psv_results) = apply(metab_psv_candidates,1,paste,collapse = " & ")
# Save
save(metab_psv_results,file = "./data/mediation/metab_psv_results.Rdata")