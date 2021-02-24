library(boot)
set.seed(1017)
# Load data
setwd("/home/vigerst/MS-Thesis/")
# setwd("/Users/timvigers/Documents/School/MS Thesis")
load("./data/mediation/methyl_psv_candidates_p_05.Rdata")
load("./data/raw_data/psv_sv_dataset.Rdata")
# Outcome and adjustment variables
ia = psv$IAgroup2
SEX = psv$SEX
dr34 = psv$dr34
covariates = as.data.frame(cbind(ia,SEX,dr34))
covariates = as.data.frame(lapply(covariates,as.factor))
# Bootstrap function
boot_med = function(data,i){
  d = data[i,] # allows boot to select sample
  # Mediation models
  # Because we have a case-control design, the mediator model should be fit to 
  # only the controls to account for oversampling a rare outcome (VanderWeele pg. 28)
  m = lm(metab~methyl+age+SEX+dr34,data = d[d$ia=="control",])
  # Refit without interaction term
  c = glm(ia ~ methyl+metab+age+SEX+dr34,data = d,family = binomial("logit"))
  # Calculate direct and indirect effect (starting on VandwerWeele pg. 27)
  de = summary(c)$coefficients[2,1]
  ie = summary(m)$coefficients[2,1] * summary(c)$coefficients[3,1]
  # Proportion mediated (VandwerWeele pg.48 and per Patrick Carry)
  pmed = (exp(de)*(exp(ie)-1))/(exp(de)*exp(ie) - 1)
  return(c(de,ie,pmed))
}
# Model function
med_mods_manual = function(age,out_name,n_cores = 16,n_sims = 10000,ci.type = "bca"){
  # Iterate through all
  mediation_results = list()
  for(row in 1:nrow(methyl_psv_candidates)){
    r = methyl_psv_candidates[row,]
    methyl = psv[,r[1]]
    metab = sv[,r[2]]*1000 # Metabolite estimates are difficult to interpret on current scale
    # Dataframe for bootstrap
    df = as.data.frame(cbind(methyl,metab,age,covariates))
    # Check interaction effect - skip any with interaction for now
    c = glm(ia ~ methyl*metab+age+SEX+dr34,data = df,family = binomial("logit"))
    if (summary(c)$coefficients[7,4] < 0.05){next}
    # Check Baron and Kenny steps
    m = lm(metab~methyl+age+SEX+dr34,data = df[df$ia=="control",])
    c = glm(ia ~ methyl+metab+age+SEX+dr34,data = df,family = binomial("logit"))
    if (summary(c)$coefficients[2,4] > 0.05 | summary(c)$coefficients[3,4] > 0.05 |
        summary(m)$coefficients[2,4] > 0.05) {next}
    # Run bootstrap
    res = boot(data = df,statistic = boot_med,R = n_sims,ncpus = n_cores)
    # Confidence intervals
    de_ci = boot.ci(res,type = ci.type,index = 1)
    de_ci = de_ci[[4]][c(4,5)]
    ie_ci = boot.ci(res,type = ci.type,index = 2)
    ie_ci = ie_ci[[4]][c(4,5)]
    pmed_ci = boot.ci(res,type = ci.type,index = 3)
    pmed_ci = pmed_ci[[4]][c(4,5)]
    # Return
    ret = c(res$t0[1],de_ci,res$t0[2],ie_ci,res$t0[3],pmed_ci)
    ret = c(r[1],r[2],round(ret,5))
    mediation_results[[row]] = ret
    }
  # Format and save
  mediation_results = do.call(rbind,mediation_results)
  colnames(mediation_results) = c("Probe","Metabolite",
                                  "DE","DE Lower","DE Upper",
                                  "IE","IE Lower","IE Upper",
                                  "PM","PM Lower","PM Upper")
  save(mediation_results,file = paste0("./data/mediation/",out_name,"_manual.Rdata"))
}
# Age at PSV
med_mods_manual(age = psv$clinage,out_name = "long_med_p_05_psv_age")
# Time from PSV to SV
med_mods_manual(age = sv$clinage - psv$clinage,out_name = "long_med_p_05_delta_age")
