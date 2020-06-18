set.seed(1017)
setwd("/Users/timvigers/GitHub/MS-Thesis")
# Load cit results, data, and annotations
load("./data/networks/cits.Rdata")
load("./data/networks/pair_data.Rdata")
load("./data/candidate_selection/annotation.850K.Rdata")
gctof_anno = read.csv("./data/metabolomics/gctof.featureAnno.csv",
                      stringsAsFactors = F)
hilic_anno = read.csv("./data/metabolomics/hilic.featureAnno.csv",
                      stringsAsFactors = F)
lipid_anno = read.csv("./data/metabolomics/lipid.featureAnno.csv",
                      stringsAsFactors = F)
oxylipin_anno = read.csv("./data/metabolomics/oxylipin.featureAnno.csv",
                         stringsAsFactors = F)
vitd_anno = read.csv("./data/metabolomics/vitd.featureAnno.csv",
                     stringsAsFactors = F)
# 0 and 1 outcome
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
med = data.frame(matrix(ncol = 15,nrow = 0))
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
  names(out) = colnames(med)
  med = rbind(med,out)
}
colnames(med) = 
  c("CpG","Direction","Metabolite","DE","DE.LL","DE.UL","DE p value",
    "IE","IE.LL","IE.UL","IE p value",
    "Prop. Med.","Prop. Med. LL","Prop. Med. UL","Prop. Med. p value")
med[,4:ncol(med)] = lapply(med[,4:ncol(med)],function(x){round(as.numeric(x),4)})
# Annotate
anno$CpG = rownames(anno)

gctof_anno = gctof_anno[,c("feature_name","BinBase.name","quant.mz","ret.index")]
colnames(gctof_anno) = c("Metabolite","Name","Mass","Ret.")

lipid_anno = lipid_anno[,c("feature_name","Metabolite.name","row.m.z","row.retention.time")]
colnames(lipid_anno) = c("Metabolite","Name","Mass","Ret.")

metabs = rbind(gctof_anno,lipid_anno)

med = merge(med,anno[,c("CpG","UCSC_RefGene_Name","chr","pos")],by = "CpG")

med = merge(med,metabs,by = "Metabolite")
