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
id_wbr6mon = cut(late_infancy$frstwbr,c(-Inf,6,Inf),right = F,
                 labels = c("<6 months",">=6 months"))
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
                  "id_wbr","id_wbr6mon","id_riceoat","id_solidfruit","id_veg",
                  "id_meat","id_meat6mon")
# Cluster variables
cores = 16
# Model function
mod_fun = function(m,var){
  mod = try(lm(m ~ var + Age + Sex + CD8T +	CD4T +	NK +	Bcell +	Mono + Platform))
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
  mods = parLapply(cl,late_infancy[,1:100],function(m){mod_fun(m,iv)})
  stopCluster(cl)
  save(mods,file = save_path)
  rm(mods)
}