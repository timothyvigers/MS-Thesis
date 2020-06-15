library(mediation)
setwd("/Users/timvigers/GitHub/MS-Thesis")
# Load cit results
load("./data/networks/cits.Rdata")
load("./data/networks/pair_data.Rdata")
# For each pair significant by cit, use mediation package to estimate ACME, etc.
med = apply(cits,1,function(x){
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  temp = pair_data[,c("T1Dgroup",methyl,metab)]
  temp$T1Dgroup = ifelse(temp$T1Dgroup == "T1D control",0,1)
  if (x[,"direction"] == ">"){
    exposure = methyl
    mediator = metab
  } else {
    exposure = metab
    mediator = methyl
  }
  outform = as.formula(paste("T1Dgroup","~",exposure,"+",mediator))
  out.mod = glm(outform,family = "binomial",data = temp)
  medform = as.formula(paste(mediator,"~",exposure))
  med.mod = lm(medform,data = temp)
  m = mediate(med.mod,out.mod,treat = exposure,mediator = mediator,
              treat.value = 10,control.value = -10)
})
