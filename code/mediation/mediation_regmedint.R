library(regmedint)
library(parallel)
set.seed(1017)
# Load data
#setwd("/home/vigerst/MS-Thesis/")
setwd("/Users/timvigers/Dropbox/School/MS Thesis")
load("./data/mediation/psv_age_package_res.Rdata")
load("./data/mediation/methyl_psv_candidates_p_05.Rdata")
methyl_psv_candidates = methyl_psv_candidates[methyl_psv_candidates[,1] %in% psv_age[,"Probe"],]
load("./data/raw_data/psv_sv_dataset.Rdata")
# Outcome and adjustment variables
ia = psv$IAgroup2
SEX = psv$SEX
dr34 = psv$dr34
covariates = as.data.frame(cbind(ia,SEX,dr34))
covariates = as.data.frame(lapply(covariates,function(x){as.numeric(as.factor(x))-1}))
# Model function
med_mods = function(age,out_name,n_cores = 8){
  # Cluster
  cl = makeCluster(n_cores,type = "FORK")
  # Iterate through all
  mediation_results = apply(methyl_psv_candidates,1,function(r){
    methyl = psv[,r[1]]
    metab = sv[,r[2]]
    # Dataframe 
    df = as.data.frame(cbind(methyl,metab,age,covariates))
    # Mediation
    regmedint_obj <- 
      regmedint(data = df,
                ## Variables
                yvar = "ia",
                avar = "methyl",
                mvar = "metab",
                cvar = c("age","SEX","dr34"),
                ## Values at which effects are evaluated
                a0 = 0,
                a1 = 1,
                m_cde = 1,
                c_cond = c(1,1,1),
                ## Model types
                mreg = "linear",
                yreg = "logistic",
                ## Additional specification
                interaction = F,
                casecontrol = F)
    s = as.data.frame(summary(regmedint_obj)$summary_myreg)
    ret = c(s["tnde","est"],s["tnde","est"] - (1.96*s["tnde","se"]),s["tnde","est"] + (1.96*s["tnde","se"]),
            s["tnie","est"],s["tnie","est"] - (1.96*s["tnie","se"]),s["tnie","est"] + (1.96*s["tnie","se"]),
            s["pm","est"],s["pm","est"] - (1.96*s["pm","se"]),s["pm","est"] + (1.96*s["pm","se"]))
    ret = c(r[1],r[2],round(ret,5))
    ret
  })
  stopCluster(cl)
  # Format and save
  mediation_results = t(mediation_results)
  colnames(mediation_results) = c("Probe","Metabolite",
                                  "DE","DE Lower","DE Upper",
                                  "IE","IE Lower","IE Upper",
                                  "PM","PM Lower","PM Upper")
  # Save
  save(mediation_results,file = paste0("./data/mediation/",out_name,"_regmedint.Rdata"))
}
# Age at PSV
med_mods(age = psv$clinage,out_name = "long_med_p_05_psv_age")
