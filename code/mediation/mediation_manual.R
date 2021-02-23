library(mediation)
library(parallel)
set.seed(1017)
# Load data
setwd("/Users/timvigers/Documents/School/MS Thesis/")
load("./data/mediation/methyl_psv_candidates_p_05.Rdata")
load("./data/raw_data/psv_sv_dataset.Rdata")
# Outcome and adjustment variables
ia = factor(psv$IAgroup2)
SEX = factor(psv$SEX)
dr34 = factor(psv$dr34,labels = c("No","Yes"))
controls = which(psv$IAgroup2 == "control")
# Model function
med_mods_tim = function(age,out_name,n_cores = 8,n_sims = 10000,long = F){
  # Iterate through all
  mediation_results = apply(methyl_psv_candidates,1,function(r){
    methyl = psv[,r[1]]
    metab = sv[,r[2]]
    # Mediation models
    # Because we have a case-control design, the mediator model should be fit to 
    # only the controls to account for oversampling a rare outcome (VanderWeele pg. 28)
    m = lm(metab[controls]~methyl[controls]+age[controls]+SEX[controls]+dr34[controls])
    c = glm(ia ~ methyl*metab+age+SEX+dr34,family = binomial("logit"))
    # Calculate direct and indirect effect (starting on VandwerWeele pg. 27)
    theta_1 = summary(c)$coefficients[2,1]
    theta_2 = summary(c)$coefficients[3,1]
    
    med = mediate(m,c,treat="methyl",mediator="metab",boot = T,sims = n_sims,long = long)
  })
  stopCluster(cl)
  # Save
  save(mediation_mods,file = paste0("./data/mediation/",out_name,"_",n_sims / 1000,"k_sims.Rdata"))
}
# Age at PSV
med_mods(age = psv$clinage,out_name = "long_med_p_05_psv_age")
# Time from PSV to SV
med_mods(age = sv$clinage - psv$clinage,out_name = "long_med_p_05_delta_age")
# Time from PSV to SV
med_mods(age = psv$clinage - psv$clinage,out_name = "long_med_p_05_no_age")
