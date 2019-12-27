# Probes
load("/home/biostats_share/Norris/data/methylation/probesFromPipeline.Rdata")
probes <- as.data.frame(probesFromPipeline)
write.csv(probes,file = "/home/vigerst/MS-Thesis/data/methylation/probes.csv",row.names = F)
# Methylation
# 450K
load("/home/biostats_share/Norris/data/methylation/sesame450K.batchAdj.Mmatrix.Rdata")
k450 <- as.data.frame(M.sesame.batch)
write.csv(k450,file = "/home/vigerst/MS-Thesis/data/methylation/450k.csv",row.names = F)
# EPIC
load("/home/biostats_share/Norris/data/methylation/sesameEPIC.batchAdj.Mmatrix.Rdata")
epic <- as.data.frame(M.sesame.batch)
write.csv(epic,file = "/home/vigerst/MS-Thesis/data/methylation/epic.csv",row.names = F)
# Combined dataset
load("/home/biostats_share/Norris/data/methylation/Mmatrix.platformAdj.Rdata")
methyl <- as.data.frame(M.adj)
write.csv(methyl,file = "/home/vigerst/MS-Thesis/data/methylation/merged.csv",row.names = F)