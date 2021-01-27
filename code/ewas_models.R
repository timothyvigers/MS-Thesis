library(parallel)
library(nlme)
setwd("/home/tim/.local/share/Cryptomator/mnt/Vault/School/Statistical Genomics/Final Project")
load("./data/final_data.RData")
rm(probes)
load("./data/probesFromPipeline.Rdata")
set.seed(1017)
# Covariates
Sex = factor(df$SEX)
Visit = factor(df$Visit,levels = c("CWB","Late Infancy","Childhood"))
CD8T	= df$CD8T
CD4T = df$CD4T
NK = df$NK
Bcell	= df$Bcell
Mono = df$Mono
Gran = df$Gran # drop to make independent
ID = df$ID
Platform = factor(df$Data)
# Diet variables - breast-feeding duration, age at introduction to milk, eggs,
# meat, vegetables, fruit, gluten, and cereals.
# Age at introduction should be categorized (0=<4 months, 1= 4-5 months, and 2= ≥6 months)
labs = c("<4 months","4-5 months","6+ months")
bfdur = df$bfdur
dairy = cut(df$frstdairy,c(0,4,6,Inf),right = F,labels = labs)
egg = cut(df$frstegg,c(0,median(df$frstegg,na.rm = T),Inf),
          labels = c("0-11 months","> 11 months"))
meat = cut(df$frstmeat,c(0,4,6,Inf),right = F,labels = labs)
veg = cut(df$frstveg,c(0,4,6,Inf),right = F,labels = labs)
fruit = cut(df$frstfruit,c(0,4,6,Inf),right = F,labels = labs)
gluten = cut(df$frstwheat,c(0,4,6,Inf),right = F,labels = labs)
cereal = factor(df$id_cerealn,labels = labs)
# Remove unnecessary columns
df = df[,probesFromPipeline]
# RI only models for all probes
cores = 4
# Breastfeeding duration
cl = makeCluster(cores,type = "FORK")
bfdur_mods = parLapply(cl,df, function(m){
  mod = try(lme(m ~ Visit*bfdur + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform,
                random=~1|ID,na.action = na.omit))
  ifelse(class(mod) == "try-error",NA,return(anova(mod)))
})
stopCluster(cl)
save(bfdur_mods,file = "./data/ewas/bfdur_mods.RData")
rm(bfdur_mods)
# Dairy
cl = makeCluster(cores,type = "FORK")
dairy_mods = parLapply(cl,df, function(m){
  mod = try(lme(m ~ Visit*dairy + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform,
                random=~1|ID,na.action = na.omit))
  ifelse(class(mod) == "try-error",NA,return(anova(mod)))
})
stopCluster(cl)
save(dairy_mods,file = "./data/ewas/dairy_mods.RData")
rm(dairy_mods)
# Egg
cl = makeCluster(cores,type = "FORK")
egg_mods = parLapply(cl,df, function(m){
  mod = try(lme(m ~ Visit*egg + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform,
                random=~1|ID,na.action = na.omit))
  ifelse(class(mod) == "try-error",NA,return(anova(mod)))
})
stopCluster(cl)
save(egg_mods,file = "./data/ewas/egg_mods.RData")
rm(egg_mods)
# Meat
cl = makeCluster(cores,type = "FORK")
meat_mods = parLapply(cl,df, function(m){
  mod = try(lme(m ~ Visit*meat + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform,
                random=~1|ID,na.action = na.omit))
  ifelse(class(mod) == "try-error",NA,return(anova(mod)))
})
stopCluster(cl)
save(meat_mods,file = "./data/ewas/meat_mods.RData")
rm(meat_mods)
# Veg
cl = makeCluster(cores,type = "FORK")
veg_mods = parLapply(cl,df, function(m){
  mod = try(lme(m ~ Visit*veg + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform,
                random=~1|ID,na.action = na.omit))
  ifelse(class(mod) == "try-error",NA,return(anova(mod)))
})
stopCluster(cl)
save(veg_mods,file = "./data/ewas/veg_mods.RData")
rm(veg_mods)
# Fruit
cl = makeCluster(cores,type = "FORK")
fruit_mods = parLapply(cl,df, function(m){
  mod = try(lme(m ~ Visit*fruit + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform,
                random=~1|ID,na.action = na.omit))
  ifelse(class(mod) == "try-error",NA,return(anova(mod)))
})
stopCluster(cl)
save(fruit_mods,file = "./data/ewas/fruit_mods.RData")
rm(fruit_mods)
# Gluten
cl = makeCluster(cores,type = "FORK")
gluten_mods = parLapply(cl,df, function(m){
  mod = try(lme(m ~ Visit*gluten + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform,
                random=~1|ID,na.action = na.omit))
  ifelse(class(mod) == "try-error",NA,return(anova(mod)))
})
stopCluster(cl)
save(gluten_mods,file = "./data/ewas/gluten_mods.RData")
rm(gluten_mods)
# Cereal
cl = makeCluster(cores,type = "FORK")
cereal_mods = parLapply(cl,df, function(m){
  mod = try(lme(m ~ Visit*cereal + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform,
                random=~1|ID,na.action = na.omit))
  ifelse(class(mod) == "try-error",NA,return(anova(mod)))
})
stopCluster(cl)
save(cereal_mods,file = "./data/ewas/cereal_mods.RData")
rm(cereal_mods)

# # Simple models
# # rm(list=ls())
# library(parallel)
# library(nlme)
# setwd("/home/tim/.local/share/Cryptomator/mnt/Vault/School/Statistical Genomics/Final Project")
# load("./data/final_data.RData")
# rm(probes)
# load("./data/probesFromPipeline.Rdata")
# set.seed(1017)
# # Late infancy
# late_infancy = df[df$Visit == "Late Infancy",]
# rm(df)
# # Variables
# Sex = factor(late_infancy$SEX)
# Age = late_infancy$clinage
# CD8T	= late_infancy$CD8T
# CD4T = late_infancy$CD4T
# NK = late_infancy$NK
# Bcell	= late_infancy$Bcell
# Mono = late_infancy$Mono
# Gran = late_infancy$Gran # drop to make independent
# ID = late_infancy$ID
# Platform = factor(late_infancy$Data)
# labs = c("<4 months","4-5 months","6+ months")
# bfdur = late_infancy$bfdur
# dairy = cut(late_infancy$frstdairy,c(0,4,6,Inf),right = F,labels = labs)
# egg = cut(late_infancy$frstegg,c(0,median(late_infancy$frstegg,na.rm = T),Inf),
#           labels = c("0-11 months","> 11 months"))
# meat = cut(late_infancy$frstmeat,c(0,4,6,Inf),right = F,labels = labs)
# veg = cut(late_infancy$frstveg,c(0,4,6,Inf),right = F,labels = labs)
# fruit = cut(late_infancy$frstfruit,c(0,4,6,Inf),right = F,labels = labs)
# gluten = cut(late_infancy$frstwheat,c(0,4,6,Inf),right = F,labels = labs)
# cereal = factor(late_infancy$id_cerealn,labels = labs)
# # Cluster variables
# cores = 4
# # Remove unnecessary columns
# late_infancy = late_infancy[,probesFromPipeline]
# # # Breastfeeding duration
# # cl = makeCluster(cores,type = "FORK")
# # bfdur_infancy_mods = parLapply(cl,late_infancy, function(m){
# #   mod = try(lm(m ~ bfdur + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
# #   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# # })
# # stopCluster(cl)
# # save(bfdur_infancy_mods,file = "./data/ewas/bfdur_infancy_mods.RData")
# # rm(bfdur_infancy_mods)
# # # Dairy
# # cl = makeCluster(cores,type = "FORK")
# # dairy_infancy_mods = parLapply(cl,late_infancy, function(m){
# #   mod = try(lm(m ~ dairy + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
# #   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# # })
# # stopCluster(cl)
# # save(dairy_infancy_mods,file = "./data/ewas/dairy_infancy_mods.RData")
# # rm(dairy_infancy_mods)
# # # Egg
# # cl = makeCluster(cores,type = "FORK")
# # egg_infancy_mods = parLapply(cl,late_infancy, function(m){
# #   mod = try(lm(m ~ egg + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
# #   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# # })
# # stopCluster(cl)
# # save(egg_infancy_mods,file = "./data/ewas/egg_infancy_mods.RData")
# # rm(egg_infancy_mods)
# # # Meat
# # cl = makeCluster(cores,type = "FORK")
# # meat_infancy_mods = parLapply(cl,late_infancy, function(m){
# #   mod = try(lm(m ~ meat + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
# #   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# # })
# # stopCluster(cl)
# # save(meat_infancy_mods,file = "./data/ewas/meat_infancy_mods.RData")
# # rm(meat_infancy_mods)
# # # Veg
# # cl = makeCluster(cores,type = "FORK")
# # veg_infancy_mods = parLapply(cl,late_infancy, function(m){
# #   mod = try(lm(m ~ veg + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
# #   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# # })
# # stopCluster(cl)
# # save(veg_infancy_mods,file = "./data/ewas/veg_infancy_mods.RData")
# # rm(veg_infancy_mods)
# # Fruit
# cl = makeCluster(cores,type = "FORK")
# fruit_infancy_mods = parLapply(cl,late_infancy, function(m){
#   mod = try(lm(m ~ fruit + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
#   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# })
# stopCluster(cl)
# save(fruit_infancy_mods,file = "./data/ewas/fruit_infancy_mods.RData")
# rm(fruit_infancy_mods)
# # Gluten
# cl = makeCluster(cores,type = "FORK")
# gluten_infancy_mods = parLapply(cl,late_infancy, function(m){
#   mod = try(lm(m ~ gluten + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
#   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# })
# stopCluster(cl)
# save(gluten_infancy_mods,file = "./data/ewas/gluten_infancy_mods.RData")
# rm(gluten_infancy_mods)
# # Cereal
# cl = makeCluster(cores,type = "FORK")
# cereal_infancy_mods = parLapply(cl,late_infancy, function(m){
#   mod = try(lm(m ~ cereal + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
#   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# })
# stopCluster(cl)
# save(cereal_infancy_mods,file = "./data/ewas/cereal_infancy_mods.RData")
# rm(cereal_infancy_mods)
# 
# # Childhood
# rm(list=ls())
# library(parallel)
# library(nlme)
# setwd("/home/tim/.local/share/Cryptomator/mnt/Vault/School/Statistical Genomics/Final Project")
# load("./data/final_data.RData")
# rm(probes)
# load("./data/probesFromPipeline.Rdata")
# set.seed(1017)
# childhood = df[df$Visit == "Childhood",]
# rm(df)
# # Variables
# Sex = factor(childhood$SEX)
# Age = childhood$clinage
# CD8T	= childhood$CD8T
# CD4T = childhood$CD4T
# NK = childhood$NK
# Bcell	= childhood$Bcell
# Mono = childhood$Mono
# Gran = childhood$Gran # drop to make independent
# ID = childhood$ID
# Platform = factor(childhood$Data)
# labs = c("<4 months","4-5 months","6+ months")
# bfdur = childhood$bfdur
# dairy = cut(childhood$frstdairy,c(0,4,6,Inf),right = F,labels = labs)
# egg = cut(childhood$frstegg,c(0,median(childhood$frstegg,na.rm = T),Inf),
#           labels = c("0-11 months","> 11 months"))
# meat = cut(childhood$frstmeat,c(0,4,6,Inf),right = F,labels = labs)
# veg = cut(childhood$frstveg,c(0,4,6,Inf),right = F,labels = labs)
# fruit = cut(childhood$frstfruit,c(0,4,6,Inf),right = F,labels = labs)
# gluten = cut(childhood$frstwheat,c(0,4,6,Inf),right = F,labels = labs)
# cereal = factor(childhood$id_cerealn,labels = labs)
# # Cluster variables
# cores = 4
# # Remove unnecessary columns
# childhood = childhood[,probesFromPipeline]
# # Breastfeeding duration
# cl = makeCluster(cores,type = "FORK")
# bfdur_childhood_mods = parLapply(cl,childhood, function(m){
#   mod = try(lm(m ~ bfdur + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
#   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# })
# stopCluster(cl)
# save(bfdur_childhood_mods,file = "./data/ewas/bfdur_childhood_mods.RData")
# rm(bfdur_childhood_mods)
# # Dairy
# cl = makeCluster(cores,type = "FORK")
# dairy_childhood_mods = parLapply(cl,childhood, function(m){
#   mod = try(lm(m ~ dairy + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
#   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# })
# stopCluster(cl)
# save(dairy_childhood_mods,file = "./data/ewas/dairy_childhood_mods.RData")
# rm(dairy_childhood_mods)
# # Egg
# cl = makeCluster(cores,type = "FORK")
# egg_childhood_mods = parLapply(cl,childhood, function(m){
#   mod = try(lm(m ~ egg + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
#   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# })
# stopCluster(cl)
# save(egg_childhood_mods,file = "./data/ewas/egg_childhood_mods.RData")
# rm(egg_childhood_mods)
# # Meat
# cl = makeCluster(cores,type = "FORK")
# meat_childhood_mods = parLapply(cl,childhood, function(m){
#   mod = try(lm(m ~ meat + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
#   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# })
# stopCluster(cl)
# save(meat_childhood_mods,file = "./data/ewas/meat_childhood_mods.RData")
# rm(meat_childhood_mods)
# # Veg
# cl = makeCluster(cores,type = "FORK")
# veg_childhood_mods = parLapply(cl,childhood, function(m){
#   mod = try(lm(m ~ veg + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
#   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# })
# stopCluster(cl)
# save(veg_childhood_mods,file = "./data/ewas/veg_childhood_mods.RData")
# rm(veg_childhood_mods)
# # Fruit
# cl = makeCluster(cores,type = "FORK")
# fruit_childhood_mods = parLapply(cl,childhood, function(m){
#   mod = try(lm(m ~ fruit + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
#   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# })
# stopCluster(cl)
# save(fruit_childhood_mods,file = "./data/ewas/fruit_childhood_mods.RData")
# rm(fruit_childhood_mods)
# # Gluten
# cl = makeCluster(cores,type = "FORK")
# gluten_childhood_mods = parLapply(cl,childhood, function(m){
#   mod = try(lm(m ~ gluten + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
#   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# })
# stopCluster(cl)
# save(gluten_childhood_mods,file = "./data/ewas/gluten_childhood_mods.RData")
# rm(gluten_childhood_mods)
# # Cereal
# cl = makeCluster(cores,type = "FORK")
# cereal_childhood_mods = parLapply(cl,childhood, function(m){
#   mod = try(lm(m ~ cereal + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
#   ifelse(class(mod) == "try-error",NA,return(anova(mod)))
# })
# stopCluster(cl)
# save(cereal_childhood_mods,file = "./data/ewas/cereal_childhood_mods.RData")
# rm(cereal_childhood_mods)
