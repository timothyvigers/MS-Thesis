library(parallel)
library(nlme)
setwd("/home/vigerst/EWAS")
load("./data/final_data.RData")
load("./data/probesFromPipeline.Rdata")
set.seed(1017)
# Late infancy
late_infancy = df[df$Visit == "Late Infancy",]
rm(df)
# Variables
Sex = factor(late_infancy$SEX)
Age = late_infancy$clinage
CD8T	= late_infancy$CD8T
CD4T = late_infancy$CD4T
NK = late_infancy$NK
Bcell	= late_infancy$Bcell
Mono = late_infancy$Mono
Gran = late_infancy$Gran # drop to make independent
ID = late_infancy$ID
Platform = factor(late_infancy$Data)
# Fix continuous variables
month_vars = c("frstdairy","frstwbr","frstro","frstsolidfruit","frstveg","frstmeat","frstwheat","frstbarley")
late_infancy[,month_vars] = lapply(late_infancy[,month_vars],function(c){
  c[which(c==-999)] = 17
  c = c - 1
})
# Make vectors for models
labs = c("<4 months","4-5 months","6+ months")
exbf = late_infancy$exbf
bfdur = late_infancy$bfdur
bfwhbar = factor(late_infancy$bfwhbar,levels = c("n","y"),labels = c("N","Y"))
frstdairy = late_infancy$frstdairy
id_solidfood = factor(late_infancy$id_solidfood,labels = labs)
id_cereal = factor(late_infancy$id_cereal,labels = labs)
id_wbr = factor(late_infancy$id_wbr,labels = labs)
id_riceoat = factor(late_infancy$id_riceoat,labels = labs)
id_solidfruit = factor(late_infancy$id_solidfruit,labels = labs)
id_veg = factor(late_infancy$id_veg,labels = labs)
id_meat = factor(late_infancy$id_meat,labels = labs)
id_meat6mon = cut(late_infancy$frstmeat,c(-Inf,6,Inf),right = F,
                  labels = c("<6 months",">=6 months"))
# Remove unnecessary columns
late_infancy = late_infancy[,probesFromPipeline]
# List of variables
analysis_vars = c("exbf","bfdur","bfwhbar","frstdairy","id_solidfood","id_cereal",
                  "id_wbr","id_riceoat","id_solidfruit","id_veg","id_meat","id_meat6mon")
# Cluster variables
cores = 16
# Model function
mod_fun = function(m,var){
  mod = try(lm(m ~ var + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
  ifelse(class(mod) == "try-error",NA,return(mod))
}
# Loop through all variables
for(v in analysis_vars){
  iv = get(v)
  save_obj = paste0(v,"_mods")
  save_path = paste0("./results/",save_obj,".RData")
  cl = makeCluster(cores,type = "FORK")
  mods = parLapply(cl,late_infancy[,1:10],function(m){mod_fun(m,iv)})
  stopCluster(cl)
  save(mods,file = save_path)
}
# Childhood
rm(list=ls())
library(parallel)
library(nlme)
setwd("/home/tim/.local/share/Cryptomator/mnt/Vault/School/Statistical Genomics/Final Project")
load("./data/final_data.RData")
rm(probes)
load("./data/probesFromPipeline.Rdata")
set.seed(1017)
childhood = df[df$Visit == "Childhood",]
rm(df)
# Variables
Sex = factor(childhood$SEX)
Age = childhood$clinage
CD8T	= childhood$CD8T
CD4T = childhood$CD4T
NK = childhood$NK
Bcell	= childhood$Bcell
Mono = childhood$Mono
Gran = childhood$Gran # drop to make independent
ID = childhood$ID
Platform = factor(childhood$Data)
labs = c("<4 months","4-5 months","6+ months")
bfdur = childhood$bfdur
dairy = cut(childhood$frstdairy,c(0,4,6,Inf),right = F,labels = labs)
egg = cut(childhood$frstegg,c(0,median(childhood$frstegg,na.rm = T),Inf),
          labels = c("0-11 months","> 11 months"))
meat = cut(childhood$frstmeat,c(0,4,6,Inf),right = F,labels = labs)
veg = cut(childhood$frstveg,c(0,4,6,Inf),right = F,labels = labs)
fruit = cut(childhood$frstfruit,c(0,4,6,Inf),right = F,labels = labs)
gluten = cut(childhood$frstwheat,c(0,4,6,Inf),right = F,labels = labs)
cereal = factor(childhood$id_cerealn,labels = labs)
# Cluster variables
cores = 4
# Remove unnecessary columns
childhood = childhood[,probesFromPipeline]
# Breastfeeding duration
cl = makeCluster(cores,type = "FORK")
bfdur_childhood_mods = parLapply(cl,childhood, function(m){
  mod = try(lm(m ~ bfdur + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
  ifelse(class(mod) == "try-error",NA,return(summary(mod)$coefficients))
})
stopCluster(cl)
save(bfdur_childhood_mods,file = "./data/ewas/bfdur_childhood_mods.RData")
rm(bfdur_childhood_mods)
# Dairy
cl = makeCluster(cores,type = "FORK")
dairy_childhood_mods = parLapply(cl,childhood, function(m){
  mod = try(lm(m ~ dairy + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
  ifelse(class(mod) == "try-error",NA,return(summary(mod)$coefficients))
})
stopCluster(cl)
save(dairy_childhood_mods,file = "./data/ewas/dairy_childhood_mods.RData")
rm(dairy_childhood_mods)
# Egg
cl = makeCluster(cores,type = "FORK")
egg_childhood_mods = parLapply(cl,childhood, function(m){
  mod = try(lm(m ~ egg + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
  ifelse(class(mod) == "try-error",NA,return(summary(mod)$coefficients))
})
stopCluster(cl)
save(egg_childhood_mods,file = "./data/ewas/egg_childhood_mods.RData")
rm(egg_childhood_mods)
# Meat
cl = makeCluster(cores,type = "FORK")
meat_childhood_mods = parLapply(cl,childhood, function(m){
  mod = try(lm(m ~ meat + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
  ifelse(class(mod) == "try-error",NA,return(summary(mod)$coefficients))
})
stopCluster(cl)
save(meat_childhood_mods,file = "./data/ewas/meat_childhood_mods.RData")
rm(meat_childhood_mods)
# Veg
cl = makeCluster(cores,type = "FORK")
veg_childhood_mods = parLapply(cl,childhood, function(m){
  mod = try(lm(m ~ veg + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
  ifelse(class(mod) == "try-error",NA,return(summary(mod)$coefficients))
})
stopCluster(cl)
save(veg_childhood_mods,file = "./data/ewas/veg_childhood_mods.RData")
rm(veg_childhood_mods)
# Fruit
cl = makeCluster(cores,type = "FORK")
fruit_childhood_mods = parLapply(cl,childhood, function(m){
  mod = try(lm(m ~ fruit + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
  ifelse(class(mod) == "try-error",NA,return(summary(mod)$coefficients))
})
stopCluster(cl)
save(fruit_childhood_mods,file = "./data/ewas/fruit_childhood_mods.RData")
rm(fruit_childhood_mods)
# Gluten
cl = makeCluster(cores,type = "FORK")
gluten_childhood_mods = parLapply(cl,childhood, function(m){
  mod = try(lm(m ~ gluten + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
  ifelse(class(mod) == "try-error",NA,return(summary(mod)$coefficients))
})
stopCluster(cl)
save(gluten_childhood_mods,file = "./data/ewas/gluten_childhood_mods.RData")
rm(gluten_childhood_mods)
# Cereal
cl = makeCluster(cores,type = "FORK")
cereal_childhood_mods = parLapply(cl,childhood, function(m){
  mod = try(lm(m ~ cereal + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
  ifelse(class(mod) == "try-error",NA,return(summary(mod)$coefficients))
})
stopCluster(cl)
save(cereal_childhood_mods,file = "./data/ewas/cereal_childhood_mods.RData")
rm(cereal_childhood_mods)
