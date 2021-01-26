library(dplyr)
setwd("/home/tim/.local/share/Cryptomator/mnt/Vault/School/Statistical Genomics/Final Project/data/ewas")
# Simple models
# Function
get_candidates = function(mods,cutoff = 0.1){
  mods = mods[!is.na(mods)]
  t = lapply(1:length(mods), function(x){
    m = mods[[x]]
    r = as.data.frame(m[1,])
    r["probe"] = names(mods)[x]
    return(r)
  })
  t = as.data.frame(do.call(rbind,t))
  t$p_adj = p.adjust(as.numeric(t$`Pr(>F)`),"fdr")
  t = t[t$p_adj < cutoff,]
  return(t$probe[order(t$p_adj)])
}
# Infancy - FDR < 0.05
# BF duration
load("./bfdur_infancy_mods.RData")
bfdur_candidates = get_candidates(bfdur_infancy_mods,0.05)
# Dairy
load("./dairy_infancy_mods.RData")
dairy_candidates = get_candidates(dairy_infancy_mods,0.05)
# Egg
load("./egg_infancy_mods.RData")
egg_candidates = get_candidates(egg_infancy_mods,0.05)
# Meat
load("./meat_infancy_mods.RData")
meat_candidates = get_candidates(meat_infancy_mods,0.05)
# Veg
load("./veg_infancy_mods.RData")
veg_candidates = get_candidates(veg_infancy_mods,0.05)
# Fruit
load("./fruit_infancy_mods.RData")
fruit_candidates = get_candidates(fruit_infancy_mods,0.05)
# Gluten
load("./gluten_infancy_mods.RData")
gluten_candidates = get_candidates(gluten_infancy_mods,0.05)
# Cereal
load("./cereal_infancy_mods.RData")
cereal_candidates = get_candidates(cereal_infancy_mods,0.05)
# Candidate table
n = max(length(bfdur_candidates),length(dairy_candidates),length(egg_candidates),
        length(meat_candidates),length(veg_candidates),length(fruit_candidates),
        length(gluten_candidates),length(cereal_candidates))
length(bfdur_candidates) = n
length(dairy_candidates) = n
length(egg_candidates) = n
length(meat_candidates) = n
length(veg_candidates) = n
length(fruit_candidates) = n
length(gluten_candidates) = n
length(cereal_candidates) = n
late_infancy_candidates = 
  bind_cols(bfdur_candidates,dairy_candidates,egg_candidates,meat_candidates,
            veg_candidates,fruit_candidates,gluten_candidates,cereal_candidates)
colnames(late_infancy_candidates) = c("bfdur","dairy","egg","meat","veg",
                                      "fruit","gluten","cereal")
write.csv(late_infancy_candidates,"late_infancy_candidates.csv",row.names = F,
          na = "")
# Childhood
# BF duration
load("./bfdur_childhood_mods.RData")
bfdur_candidates = get_candidates(bfdur_childhood_mods,0.05)
# Dairy
load("./dairy_childhood_mods.RData")
dairy_candidates = get_candidates(dairy_childhood_mods,0.05)
# Egg
load("./egg_childhood_mods.RData")
egg_candidates = get_candidates(egg_childhood_mods,0.05)
# Meat
load("./meat_childhood_mods.RData")
meat_candidates = get_candidates(meat_childhood_mods,0.05)
# Veg
load("./veg_childhood_mods.RData")
veg_candidates = get_candidates(veg_childhood_mods,0.05)
# Fruit
load("./fruit_childhood_mods.RData")
fruit_candidates = get_candidates(fruit_childhood_mods,0.05)
# Gluten
load("./gluten_childhood_mods.RData")
gluten_candidates = get_candidates(gluten_childhood_mods,0.05)
# Cereal
load("./cereal_childhood_mods.RData")
cereal_candidates = get_candidates(cereal_childhood_mods,0.05)
# Candidate table
n = max(length(bfdur_candidates),length(dairy_candidates),length(egg_candidates),
        length(meat_candidates),length(veg_candidates),length(fruit_candidates),
        length(gluten_candidates),length(cereal_candidates))
length(bfdur_candidates) = n
length(dairy_candidates) = n
length(egg_candidates) = n
length(meat_candidates) = n
length(veg_candidates) = n
length(fruit_candidates) = n
length(gluten_candidates) = n
length(cereal_candidates) = n
late_childhood_candidates = 
  bind_cols(bfdur_candidates,dairy_candidates,egg_candidates,meat_candidates,
            veg_candidates,fruit_candidates,gluten_candidates,cereal_candidates)
colnames(late_childhood_candidates) = c("bfdur","dairy","egg","meat","veg",
                                      "fruit","gluten","cereal")
write.csv(late_childhood_candidates,"late_childhood_candidates.csv",row.names = F,
          na = "")