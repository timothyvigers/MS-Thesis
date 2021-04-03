library(parallel)
setwd("/home/vigerst/EWAS")
load("./data/final_data.RData")
load("./data/probesFromPipeline.Rdata")
set.seed(1017)
# Childhood data
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
# Fix continuous variables
month_vars = c("frstdairy","frstwbr","frstro","frstsolidfruit","frstveg","frstmeat","frstwheat","frstbarley")
childhood[,month_vars] = lapply(childhood[,month_vars],function(c){
  c[which(c==-999)] = 17
  c = c - 1
})
# Make vectors for models
labs = c("<4 months","4-5 months","6+ months")
exbf = childhood$exbf
bfdur = childhood$bfdur
bfwhbar = factor(childhood$bfwhbar,levels = c("n","y"),labels = c("N","Y"))
frstdairy = childhood$frstdairy
id_solidfood = factor(childhood$id_solidfood,labels = labs)
id_cereal = factor(childhood$id_cereal,labels = labs)
id_wbr = factor(childhood$id_wbr,labels = labs)
id_wbr6mon = cut(childhood$frstwbr,c(-Inf,6,Inf),right = F,
                 labels = c("<6 months",">=6 months"))
id_riceoat = factor(childhood$id_riceoat,labels = labs)
id_solidfruit = factor(childhood$id_solidfruit,labels = labs)
id_veg = factor(childhood$id_veg,labels = labs)
id_meat = factor(childhood$id_meat,levels = c(0,1,2),labels = labs)
id_meat6mon = cut(childhood$frstmeat,c(-Inf,6,Inf),right = F,
                  labels = c("<6 months",">=6 months"))
# Remove unnecessary columns
childhood = childhood[,probesFromPipeline]
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
  save_path = paste0("./results/childhood/",save_obj,".RData")
  cl = makeCluster(cores,type = "FORK")
  mods = parLapply(cl,childhood[,1:100],function(p){mod_fun(p,iv)})
  stopCluster(cl)
  save(mods,file = save_path)
  rm(mods)
}
