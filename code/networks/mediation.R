setwd("/Users/timvigers/GitHub/MS-Thesis")
# Load cit results and data
load("./data/networks/cits.Rdata")
load("./data/networks/pair_data.Rdata")
pair_data$T1Dgroup = ifelse(pair_data$T1Dgroup == "T1D control",0,1)
# Weights for case-control study (because outcome is not rare in this cohort)
p = sum(pair_data$T1Dgroup) / nrow(pair_data)
# Pi is T1D prevalence
pi = 0.03847
pair_data$weight = ifelse(pair_data$T1Dgroup == 0,(1-pi)/(1-p),pi/p)
# Mediation function
mediate <- function(data,i){
  data = data[i,]
  outform = as.formula(paste("T1Dgroup","~",exposure,"+",mediator))
  out.mod = glm(outform,family = "binomial",data = data)
  medform = as.formula(paste(mediator,"~",exposure))
  med.mod = lm(medform,data = data,weights = data$weight)
  # Model estimates
  theta0 = summary(out.mod)$coefficients[1,1]
  theta1 = summary(out.mod)$coefficients[2,1]
  theta2 = summary(out.mod)$coefficients[3,1]
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
  return(c(cde,cie,pmed))
}
# For each pair significant by cit, use mediation package to estimate ACME, etc.
med = data.frame(matrix(ncol = 12,nrow = 0))
for (r in 1:nrow(cits)) {
  x = cits[r,]
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
  out.mod = glm(outform,family = "binomial",data = data)
  medform = as.formula(paste(mediator,"~",exposure))
  med.mod = lm(medform,data = data,weights = data$weight)
  intform = as.formula(paste("T1Dgroup","~",exposure,"*",mediator))
  int.mod = glm(intform,family = "binomial",data = data)
  if(any(summary(out.mod)$coefficients[2:3,4]>0.05)|
     summary(med.mod)$coefficients[2,4]>0.05|
     summary(int.mod)$coefficients[4,4]<0.05){
    next()
  }
  b <- boot(pair_data,mediate,R=10)
  # Bootstrap CIs
  c(cde,cie,pmed)
  cde.ci = boot.ci(b,conf = 0.95, type = c("norm"), index=1)
  cie.ci = boot.ci(b,conf = 0.95, type = c("norm"), index=2)
  pmed.ci = boot.ci(b,conf = 0.95, type = c("norm"), index=3)
  out = c(methyl,metab,
          cde.ci$t0,cde.ci$normal[2],cde.ci$normal[3],pnorm(abs((b$t0[1] - mean(b$t[,1]) )) / sqrt(var(b$t[, 1])), lower.tail=F)*2,
          cie.ci$t0,cie.ci$normal[2],cie.ci$normal[3],
          pmed.ci$t0,pmed.ci$normal[2],pmed.ci$normal[3])
  names(out) = colnames(med)
  med = rbind(med,out)
}
colnames(med) = 
  c("CpG","Metabolite","DE","DE.LL","DE.UL","p","IE","IE.LL","IE.UL",
    "Prop. Med.","Prop. Med. LL","Prop. Med. UL")
# Questions
# 1. To extend this when there are interaction effects, do we need to pick a level 
# of the mediator (m) to evaluate at? Would the mean be best?
