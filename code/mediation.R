library(boot)
set.seed(1017)
setwd("/home/vigerst/MS-Thesis/")
# Load cit results, data, and annotations
load("./data/networks/cits.Rdata")
load("./data/networks/pair_data.Rdata")
load("./data/candidate_selection/annotation.850K.Rdata")
gctof_anno = read.csv("./data/metabolomics/gctof.featureAnno.csv",
                      stringsAsFactors = F,na.strings = "")
lipid_anno = read.csv("./data/metabolomics/lipid.featureAnno.csv",
                      stringsAsFactors = F,na.strings = "")
# 0 and 1 outcome
pair_data$T1Dgroup = ifelse(pair_data$T1Dgroup == "T1D control",0,1)
# Weights for case-control study (because outcome more rare than in this cohort)
p = sum(pair_data$T1Dgroup) / nrow(pair_data)
# Pi is T1D prevalence
pi = 0.03847
pair_data$weight = ifelse(pair_data$T1Dgroup == 0,(1-pi)/(1-p),pi/p)
# Mediation function
mediate <- function(data,i){
  data = data[i,]
  # Regress outcome on exposure and mediator
  outform = as.formula(paste("T1Dgroup","~",exposure,"+",mediator))
  out.mod = glm(outform,family = "binomial",data = data)
  # Regress mediator on exposure
  medform = as.formula(paste(mediator,"~",exposure))
  med.mod = lm(medform,data = data,weights = data$weight)
  # Get model estimates
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
  return(c(cde,cie,pmed))
}
# For each pair significant by cit, bootstrap mediation estimates
# Results dataframe
med = list()
result_names = 
  c("CpG","Direction","Metabolite","DE","DE.LL","DE.UL","DE p value",
    "IE","IE.LL","IE.UL","IE p value",
    "Prop. Med.","Prop. Med. LL","Prop. Med. UL","Prop. Med. p value")
for (r in 1:nrow(cits)) {
  x = cits[r,]
  d = cits[r,"direction"]
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  if (x["direction"] == ">"){
    exposure = methyl
    mediator = metab
  } else {
    exposure = metab
    mediator = methyl
  }
  # Check mediation assumptions. If violated, next.
  outform = as.formula(paste("T1Dgroup","~",exposure,"+",mediator))
  out.mod = glm(outform,family = "binomial",data = pair_data)
  medform = as.formula(paste(mediator,"~",exposure))
  med.mod = lm(medform,data = pair_data,weights = pair_data$weight)
  intform = as.formula(paste("T1Dgroup","~",exposure,"*",mediator))
  int.mod = glm(intform,family = "binomial",data = pair_data)
  if(any(summary(out.mod)$coefficients[2:3,4]>0.05)|
     summary(med.mod)$coefficients[2,4]>0.05|
     summary(int.mod)$coefficients[4,4]<0.05){
    next()
  }
  # Bootstrap data
  b <- boot(pair_data,mediate,R=1000)
  # Bootstrap CIs
  cde.ci = boot.ci(b,conf = 0.95, type = c("norm"), index=1)
  cie.ci = boot.ci(b,conf = 0.95, type = c("norm"), index=2)
  pmed.ci = boot.ci(b,conf = 0.95, type = c("norm"), index=3)
  # Results
  out = c(methyl,d,metab,
          cde.ci$t0,cde.ci$normal[2],cde.ci$normal[3],
          pnorm(abs((2*b$t0[1] - mean(b$t[,1]) )) / sqrt(var(b$t[, 1])), 
                lower.tail=F)*2,
          cie.ci$t0,cie.ci$normal[2],cie.ci$normal[3],
          pnorm(abs((2*b$t0[2] - mean(b$t[,2]) )) / sqrt(var(b$t[, 2])), 
                lower.tail=F)*2,
          pmed.ci$t0,pmed.ci$normal[2],pmed.ci$normal[3],
          pnorm(abs((2*b$t0[3] - mean(b$t[,3]) )) / sqrt(var(b$t[, 3])), 
                lower.tail=F)*2)
  names(out) = result_names
  # Add to results dataframe
  med[[r]] = out
}
med = as.data.frame(do.call(rbind,med))
# Format results
med[,4:ncol(med)] = 
  lapply(med[,4:ncol(med)],function(x){
    round(as.numeric(as.character(x)),4)
  })
# Annotate
anno$CpG = rownames(anno)

gctof_anno = gctof_anno[,c("feature_name","BinBase.name","quant.mz","ret.index")]
colnames(gctof_anno) = c("Metabolite","Name","Mass","Ret.")

lipid_anno = lipid_anno[,c("feature_name","Metabolite.name","row.m.z","row.retention.time")]
colnames(lipid_anno) = c("Metabolite","Name","Mass","Ret.")

metabs = rbind(gctof_anno,lipid_anno)

med = merge(med,anno[,c("CpG","UCSC_RefGene_Name","chr","pos")],
            by = "CpG",sort = F)

med = merge(med,metabs,by = "Metabolite",sort = F)

med = med[,c(2,3,1,4:ncol(med))]

# Save
save(med,file = "./data/mediation.Rdata")
