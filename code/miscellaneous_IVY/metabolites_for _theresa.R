metabs = c("linolenic acid","linoleic acid","arachidonic acid","eicosapentaenoic acid",
           "docosahexaenoic acid","Docosapentaenoic acid")
lipid_anno = read.csv("/Users/timvigers/Dropbox/School/MS Thesis/data/raw_data/lipid.featureAnno.csv")
matches = lipid_anno$feature_name[unique(grep(paste(metabs,collapse="|"),lipid_anno$Metabolite.name))]
lipids = read.csv("/Users/timvigers/Dropbox/School/MS Thesis/data/raw_data/lipid.bc.csv")
lipids = lipids[,c("samplekey",matches)]
write.csv(lipids,file = "~/Desktop/lipids_for_theresa.csv",na="",row.names = F)
