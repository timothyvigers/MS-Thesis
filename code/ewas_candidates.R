library(tibble)
library(dplyr)
setwd("/home/vigerst/EWAS")
# Annotation  
annotate = c("UCSC_RefGene_Name","chr","pos","Probe_rs","Islands_Name","Relation_to_Island")
load("/home/vigerst/MS-Thesis/data/raw_data/annotation.450k.Rdata")
k450 = as.data.frame(anno)
rm(anno)
load("/home/vigerst/MS-Thesis/data/raw_data/annotation.850K.Rdata")
k850 = as.data.frame(anno)
rm(anno)
anno = rbind(k450[,annotate],k850[,annotate])
rm("k450","k850")
load("./data/probesFromPipeline.Rdata")
# List of variables
analysis_vars = c("exbf","bfdur","bfwhbar","frstdairy","id_solidfood","id_cereal",
                  "id_wbr","id_wbr6mon","id_riceoat","id_solidfruit","id_veg",
                  "id_meat","id_meat6mon")
# Format function
format_mods = function(t){
  t = t[,c("Estimate","Pr(>|t|)")]
  n = rownames(t)
  n = n[grep("var",n)]
  t = t[grep("var",rownames(t)),]
  t = c(t(t))
  names(t) = paste(rep(n,each = length(n)),c("Beta","P"))
  names(t) = gsub("var",v,names(t))
  return(t)
}
# Candidate function
format_candidates = function(timepoint){
  for(v in analysis_vars){
    load(paste0("./results/",timepoint,"/",v,"_mods.RData"))
    mods = mods[!is.na(mods)]
    candidates = lapply(mods,format_candidate)
    candidates = as.data.frame(do.call(rbind,candidates))
    candidates = candidates %>% rownames_to_column("probe")
    # Add adjusted p value columns
    candidates = candidates %>% 
      mutate(across(ends_with("P"), ~p.adjust(.,"fdr"), .names = "{col}.FDR"))
    # Annotate
    candidate_anno = anno[match(candidates$probe,rownames(anno)),]
    candidates = cbind(candidates,candidate_anno)
    # Save
    save_obj = paste0(v,"_candidates")
    save_path = paste0("./results/",timepoint,"/",save_obj,".csv")
    write.csv(candidates,file = save_path,row.names = F,na="")
  }
}
# Late infancy
format_candidates("late_infancy")
# Childhood
format_candidates("childhood")