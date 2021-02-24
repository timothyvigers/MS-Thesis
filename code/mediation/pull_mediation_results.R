library(mediation)
setwd("/Users/timvigers/Documents/School/MS Thesis")
load("./data/mediation/methyl_psv_candidates_p_05.Rdata")

load("./data/mediation/long_med_p_05_psv_age_10k_sims.Rdata")
names(mediation_mods) = methyl_psv_candidates[,1]
psv_age_package_res = lapply(names(mediation_mods), function(n){
  s = summary(mediation_mods[[n]])
  ret = c(s$d.avg,s$d.avg.ci[1],s$d.avg.ci[2],s$z.avg,s$z.avg.ci[1],s$z.avg.ci[2],
          s$n.avg,s$n.avg.ci[1],s$n.avg.ci[2])
  return(c(n,"lipid_582",round(ret,5)))
})
psv_age_package_res = do.call(rbind,psv_age_package_res)
colnames(psv_age_package_res) = c("Probe","Metabolite","DE","DE Lower","DE Upper",
                                  "IE","IE Lower","IE Upper","PM","PM Lower","PM Upper")
save(psv_age_package_res,file = "./data/mediation/psv_age_package_res.Rdata")
rm(mediation_mods,psv_age_package_res)

load("./data/mediation/long_med_p_05_delta_age_10k_sims.Rdata")
names(mediation_mods) = methyl_psv_candidates[,1]
delta_age_package_res = lapply(names(mediation_mods), function(n){
  s = summary(mediation_mods[[n]])
  ret = c(s$d.avg,s$d.avg.ci[1],s$d.avg.ci[2],s$z.avg,s$z.avg.ci[1],s$z.avg.ci[2],
          s$n.avg,s$n.avg.ci[1],s$n.avg.ci[2])
  return(c(n,"lipid_582",round(ret,5)))
})
delta_age_package_res = do.call(rbind,delta_age_package_res)
colnames(delta_age_package_res) = c("Probe","Metabolite","DE","DE Lower","DE Upper",
                                    "IE","IE Lower","IE Upper","PM","PM Lower","PM Upper")
save(delta_age_package_res,file = "./data/mediation/delta_age_package_res.Rdata")
rm(mediation_mods,delta_age_package_res)