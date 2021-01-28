setwd("/home/tim/.local/share/Cryptomator/mnt/Vault/School/Statistical Genomics/Final Project/data")
# Annotation  
annotate = c("UCSC_RefGene_Name","chr","pos","Probe_rs","Islands_Name","Relation_to_Island")
load("annotation.450k.Rdata")
k450 = as.data.frame(anno)
rm(anno)
load("annotation.850K.Rdata")
k850 = as.data.frame(anno)
rm(anno)
anno = rbind(k450[,annotate],k850[,annotate])
rm("k450","k850")
# Function
format_candidate = function(diet_var,mods,probe){
  t = mods[[probe]]
  t = t[,c("Estimate","Pr(>|t|)")]
  n = rownames(t)
  n = n[grep(diet_var,n)]
  t = t[grep(diet_var,rownames(t)),]
  t = c(t(t))
  names(t) = paste(rep(n,each = length(n)),c("Beta","P"))
  t
}
# Infancy
load("./ewas/bfdur_infancy_mods.RData")
load("./ewas/dairy_infancy_mods.RData")
load("./ewas/egg_infancy_mods.RData")
load("./ewas/meat_infancy_mods.RData")
load("./ewas/veg_infancy_mods.RData")
load("./ewas/fruit_infancy_mods.RData")
load("./ewas/gluten_infancy_mods.RData")
load("./ewas/cereal_infancy_mods.RData")
candidates = read.csv("./ewas/late_infancy_candidates.csv",na.strings = "")
candidates = as.character(unlist(candidates))
candidates = candidates[!is.na(candidates)]
infancy_candidates = lapply(candidates, function(p){
  b = format_candidate("bfdur",bfdur_infancy_mods,p)
  d = format_candidate("dairy",dairy_infancy_mods,p)
  e = format_candidate("egg",egg_infancy_mods,p)
  m = format_candidate("meat",meat_infancy_mods,p)
  v = format_candidate("veg",veg_infancy_mods,p)
  f = format_candidate("fruit",fruit_infancy_mods,p)
  g = format_candidate("gluten",gluten_infancy_mods,p)
  c = format_candidate("cereal",cereal_infancy_mods,p)
  r = c(b,d,e,m,v,f,g,c)
})
infancy_candidates = do.call(rbind,infancy_candidates)
rownames(infancy_candidates) = candidates
# Annotate
infancy_anno = anno[match(rownames(infancy_candidates),rownames(anno)),]
infancy_candidates = cbind(infancy_candidates,infancy_anno)
write.csv(infancy_candidates,"./ewas/late_infancy_model_results.csv",na = "")
# Childhood
load("./ewas/bfdur_childhood_mods.RData")
load("./ewas/dairy_childhood_mods.RData")
load("./ewas/egg_childhood_mods.RData")
load("./ewas/meat_childhood_mods.RData")
load("./ewas/veg_childhood_mods.RData")
load("./ewas/fruit_childhood_mods.RData")
load("./ewas/gluten_childhood_mods.RData")
load("./ewas/cereal_childhood_mods.RData")
candidates = read.csv("./ewas/childhood_candidates.csv",na.strings = "")
candidates = as.character(unlist(candidates))
candidates = candidates[!is.na(candidates)]
childhood_candidates = lapply(candidates, function(p){
  b = format_candidate("bfdur",bfdur_childhood_mods,p)
  d = format_candidate("dairy",dairy_childhood_mods,p)
  e = format_candidate("egg",egg_childhood_mods,p)
  m = format_candidate("meat",meat_childhood_mods,p)
  v = format_candidate("veg",veg_childhood_mods,p)
  f = format_candidate("fruit",fruit_childhood_mods,p)
  g = format_candidate("gluten",gluten_childhood_mods,p)
  c = format_candidate("cereal",cereal_childhood_mods,p)
  r = c(b,d,e,m,v,f,g,c)
})
childhood_candidates = do.call(rbind,childhood_candidates)
rownames(childhood_candidates) = candidates
# Annotate
childhood_anno = anno[match(rownames(childhood_candidates),rownames(anno)),]
childhood_candidates = cbind(childhood_candidates,childhood_anno)
write.csv(childhood_candidates,"./ewas/childhood_model_results.csv",na = "")
