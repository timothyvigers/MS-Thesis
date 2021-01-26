load("/home/tim/Documents/GitHub/MS-Thesis/data/mediation/longitudinal_mediation.Rdata")
load("/home/tim/Documents/GitHub/MS-Thesis/data/networks/pair_list.Rdata")

metabolites = unique(pairs$metab)
probes = colnames(all_data)[grepl("cg.*",colnames(all_data))]

# Metabolites PSV, methylation SV
metab_psv = all_data[all_data$Visit_Type == "PSV",
                     c("samplekey","IAgroup2","T1Dgroup",metabolites)]
metab_psv = metab_psv[!(rowSums(is.na(metab_psv[,metabolites]))==length(metabolites)),]
methyl_sv = all_data[all_data$samplekey %in% metab_psv$samplekey,
                     c("samplekey","IAgroup2","T1Dgroup",probes)]
methyl_sv = methyl_sv[!(rowSums(is.na(methyl_sv[,probes]))==length(probes)),]
# Metabolites SV, methylation PSV
metab_sv = all_data[all_data$Visit_Type == "SV",
                     c("samplekey","IAgroup2","T1Dgroup","Visit_Type",metabolites)]
metab_sv = metab_sv[!(rowSums(is.na(metab_sv[,metabolites]))==length(metabolites)),]

methyl_psv = all_data[all_data$samplekey %in% metab_sv$samplekey,
                     c("samplekey","IAgroup2","T1Dgroup",probes)]
methyl_psv = methyl_psv[!(rowSums(is.na(methyl_psv[,probes]))==length(probes)),]
