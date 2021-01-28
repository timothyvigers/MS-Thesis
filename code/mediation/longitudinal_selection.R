# Data import
# Phenotype
pheno <- read.csv("/home/biostats_share/Norris/data/phenotype/ivyomicssample_noIdentifyingInfo.csv",
                  stringsAsFactors = F,
                  na.strings = "")
# Probes
load("/home/biostats_share/Norris/data/methylation/probesFromPipeline.Rdata")
# Methylation
load("/home/biostats_share/Norris/data/methylation/Mmatrix.platformAdj.regressOut.Rdata")
methyl <- as.data.frame(t(M.adj))
# Keys
key_450k <- read.csv("/home/biostats_share/Norris/data/methylation/key.450K.csv",
                     stringsAsFactors = F)
key_450k$platform <- "450K"
key_epic <- read.csv("/home/biostats_share/Norris/data/methylation/key.EPIC.csv",
                     stringsAsFactors = F)
key_epic$platform <- "EPIC"
key <- rbind(key_450k,key_epic)
# Make final methylation dataset
methyl_samplekey <- key$samplekey[match(rownames(methyl),key$array)]
# PSV and SV samples
psv_samples = pheno$samplekey[which(pheno$Visit_Type == "PSV")]
sv_samples = pheno$samplekey[which(pheno$Visit_Type == "SV")]
# Metabolites
gctof <- read.csv("/home/biostats_share/Norris/data/metabolomics/gctof.bc.csv")
rownames(gctof) = gctof$samplekey
gctof$samplekey = NULL
hilic <- read.csv("/home/biostats_share/Norris/data/metabolomics/hilic.bc.csv")
rownames(hilic) = hilic$samplekey
hilic$samplekey = NULL
lipid <- read.csv("/home/biostats_share/Norris/data/metabolomics/lipid.bc.csv")
rownames(lipid) = lipid$samplekey
lipid$samplekey = NULL
oxylipin <- read.csv("/home/biostats_share/Norris/data/metabolomics/oxylipin.bc.csv")
rownames(oxylipin) = oxylipin$samplekey
oxylipin$samplekey = NULL
vitd <- read.csv("/home/biostats_share/Norris/data/metabolomics/vitD.bc.csv")
rownames(vitd) = vitd$samplekey
vitd$samplekey = NULL
# Model function
run_mods <- function(methyl,metab,no_cores = 20) {
  require(parallel)
  # Match indices
  psv_methyl_match = match(psv_samples,methyl_samplekey)
  psv_metab_match = match(psv_samples,rownames(metab))
  sv_methyl_match = match(sv_samples,methyl_samplekey)
  sv_metab_match = match(sv_samples,rownames(metab))
  ia = factor(pheno$IAgroup2)
  # Make cluster
  cl <- makeCluster(no_cores,type = "FORK")
  # PSV
  psv_candidates = parLapply(cl,names(methyl), function(c){
    # Check association between methylation and IA
    meth = methyl[,c]
    meth_mod = lm(meth[psv_methyl_match] ~ ia[match(psv_samples,pheno$samplekey)])
    metab_mods = lapply(names(metab), function(m){
      # Check association between metabolite and IA
      meta = metab[,m]
      meta_mod = lm(meta[psv_metab_match] ~ ia[match(psv_samples,pheno$samplekey)])
      # Check association between metabolite and methylation
      mm_mod = lm(meth[psv_methyl_match] ~ meta[psv_metab_match])
      if(summary(meta_mod)$coefficients[2,4] > 0.5 | 
         summary(meth_mod)$coefficients[2,4] > 0.5) {
        return(NA)
      } else if(summary(mod)$coefficients[2,4] < 0.4) {
        return(paste(c,m))
      } else {
          return(NA)
        }
    })
  })
  psv_candidates = unlist(psv_candidates)[!is.na(unlist(psv_candidates))]
  # SV
  sv_candidates = parLapply(cl,names(methyl), function(c){
    # Check association between methylation and IA
    meth = methyl[,c]
    meth_mod = lm(meth[sv_methyl_match] ~ ia[match(sv_samples,pheno$samplekey)])
    metab_mods = lapply(names(metab), function(m){
      # Check association between metabolite and IA
      meta = metab[,m]
      meta_mod = lm(meta[sv_metab_match] ~ ia[match(sv_samples,pheno$samplekey)])
      # Check association between metabolite and methylation
      mm_mod = lm(meth[sv_methyl_match] ~ meta[sv_metab_match])
      if(summary(meta_mod)$coefficients[2,4] > 0.5 | 
         summary(meth_mod)$coefficients[2,4] > 0.5) {
        return(NA)
      } else if(summary(mod)$coefficients[2,4] < 0.4) {
        return(paste(c,m))
      } else {
        return(NA)
      }
    })
  })
  sv_candidates = unlist(sv_candidates)[!is.na(unlist(sv_candidates))]
  stopCluster(cl)
  # Results
  return(intersect(psv_candidates,sv_candidates))
}
# gctof
gctof_candidates = run_mods(methyl = methyl[,1:20],metab = gctof)
save(gctof_candidates,file = "/home/vigerst/MS-Thesis/data/mediation/gctof_candidates.RData")
# hilic
hilic_candidates = run_mods(methyl = methyl[,1:20],metab = hilic)
save(hilic_candidates,file = "/home/vigerst/MS-Thesis/data/mediation/hilic_candidates.RData")
# lipid
lipid_candidates = run_mods(methyl = methyl[,1:20],metab = lipid)
save(lipid_candidates,file = "/home/vigerst/MS-Thesis/data/mediation/lipid_candidates.RData")
# oxylipin
oxylipin_candidates = run_mods(methyl = methyl[,1:20],metab = oxylipin)
save(oxylipin_candidates,file = "/home/vigerst/MS-Thesis/data/mediation/oxylipin_candidates.RData")
# vitd
vitd_candidates = run_mods(methyl = methyl[,1:20],metab = vitd)
save(vitd_candidates,file = "/home/vigerst/MS-Thesis/data/mediation/vitd_candidates.RData")