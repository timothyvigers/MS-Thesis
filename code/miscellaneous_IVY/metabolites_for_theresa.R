metabs = c("linolenic acid","linoleic acid","arachidonic acid","eicosapentaenoic acid",
           "docosahexaenoic acid","Docosapentaenoic acid")
lipid_anno = read.csv("~/Dropbox/School/MS Thesis/data/raw_data/lipid.featureAnno.csv")
matches = lipid_anno$feature_name[unique(grep(paste(metabs,collapse="|"),lipid_anno$Metabolite.name))]
lipids = read.csv("~/Dropbox/School/MS Thesis/data/raw_data/lipid.bc.csv")
sample = read.csv("~/Dropbox/School/MS Thesis/data/raw_data/ivyomicssample_noIdentifyingInfo.csv")
lipids = merge(lipids,sample[,c("samplekey","ID","Visit_Type")])
lipids = lipids[,c("samplekey","ID","Visit_Type",matches,"lipid_593")]

write.csv(lipids,file = "~/Desktop/lipids_for_theresa.csv",na="",row.names = F)
