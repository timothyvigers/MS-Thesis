setwd("/Users/timvigers/GitHub/MS-Thesis")
# Load cit results and data
load("./data/networks/cits.Rdata")
load("./data/networks/pair_data.Rdata")
pair_data$T1Dgroup = ifelse(pair_data$T1Dgroup == "T1D control",0,1)
# Weights for case-control study (because outcome is not rare in this cohort)
p = sum(pair_data$T1Dgroup) / nrow(pair_data)
# Pi is T1D prevalence
pi = 0.003
pair_data$weight = ifelse(pair_data$T1Dgroup == 0,(1-pi)/(1-p),pi/p)
# For each pair significant by cit, use mediation package to estimate ACME, etc.
med = apply(cits,1,function(x){
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  if (x["direction"] == ">"){
    exposure = methyl
    mediator = metab
  } else {
    exposure = metab
    mediator = methyl
  }
  outform = as.formula(paste("T1Dgroup","~",exposure,"+",mediator))
  out.mod = glm(outform,family = "binomial",data = pair_data)
  medform = as.formula(paste(mediator,"~",exposure))
  med.mod = lm(medform,data = pair_data,weights = pair_data$weight)
  # Check for interaction between exposure and mediator
  
  # Model estimates
  theta0 = summary(out.mod)$coefficients[1,1]
  theta1 = summary(out.mod)$coefficients[2,1]
  theta2 = summary(out.mod)$coefficients[3,1]
  #theta3 = summary(out.mod)$coefficients[4,1]
  beta0 = summary(med.mod)$coefficients[1,1]
  beta1 = summary(med.mod)$coefficients[2,1]
  # Conditional direct and indirect effects (VanderWeele & Vansteelandt, 2010)
  cde = theta1
  cie = theta2*beta1
  # Proportion mediated
  pmed = (exp(cde)*(exp(cie)-1))/((exp(cde)*exp(cie))-1)
  # Results
  # Check associations - if methyl or metab not associated with outcome in 
  # multivariate model, or mediator not associated with exposure return NA
  if(any(summary(out.mod)$coefficients[2:3,4]>0.05)|
     summary(med.mod)$coefficients[2,4]>0.05){
    return(rep(NA,5))
  } else {
    return(c(methyl,metab,round(c(cde,cie,pmed),3)))
  }
})
med = t(med)
med = med[complete.cases(med),]
colnames(med) = c("CpG","Metabolite","DE","IE","Prop. Mediated")
# Questions
# 1. To extend this when there are interaction effects, do we need to pick a level 
# of the mediator (m) to evaluate at? Would the mean be best?
# 2. What about treatment? The paper uses a binary treatment (a = therapy vs. 
# a* = no therapy), but it seems like the approach can be generalized. 
# 3. Why am I getting negative and > 1 proportions? Is it because the outcome 
# isn't rare in this cohort?
