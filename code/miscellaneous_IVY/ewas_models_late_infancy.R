library(parallel)
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
NHW = late_infancy$NHW
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
# For 3 group factors, make 4-5 months the reference group
id_solidfood = factor(late_infancy$id_solidfood,labels = labs)
id_solidfood = relevel(id_solidfood,ref = "4-5 months")
id_cereal = factor(late_infancy$id_cereal,labels = labs)
id_cereal = relevel(id_cereal,ref = "4-5 months")
id_wbr = factor(late_infancy$id_wbr,labels = labs)
id_wbr = relevel(id_wbr,ref = "4-5 months")
id_wbr6mon = cut(late_infancy$frstwbr,c(-Inf,6,Inf),right = F,
                 labels = c("<6 months",">=6 months"))
id_riceoat = factor(late_infancy$id_riceoat,labels = labs)
id_riceoat = relevel(id_riceoat,ref = "4-5 months")
id_solidfruit = factor(late_infancy$id_solidfruit,labels = labs)
id_solidfruit = relevel(id_solidfruit,ref = "4-5 months")
id_veg = factor(late_infancy$id_veg,labels = labs)
id_veg = relevel(id_veg,ref = "4-5 months")
id_meat = factor(late_infancy$id_meat,labels = labs)
id_meat = relevel(id_meat,ref = "4-5 months")
id_meat6mon = cut(late_infancy$frstmeat,c(-Inf,6,Inf),right = F,
                  labels = c("<6 months",">=6 months"))
# List of variables
analysis_vars = c("exbf","bfdur","bfwhbar","frstdairy","id_solidfood","id_cereal",
                  "id_wbr","id_wbr6mon","id_riceoat","id_solidfruit","id_veg",
                  "id_meat","id_meat6mon")
# Cluster variables
cores = 16
# Model function
mod_fun = function(m,var){
  mod = try(lm(m ~ var + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform + NHW))
  if(class(mod)=="lm"){
    return(summary(mod)$coefficients)
  } else {return(NA)}
}
# Loop through all variables
for(v in analysis_vars){
  iv = get(v)
  save_obj = paste0(v,"_mods")
  save_path = paste0("./results/late_infancy/",save_obj,".RData")
  cl = makeCluster(cores,type = "FORK")
  mods = parLapply(cl,late_infancy[,probesFromPipeline],function(p){mod_fun(p,iv)})
  stopCluster(cl)
  save(mods,file = save_path)
  rm(mods)
}