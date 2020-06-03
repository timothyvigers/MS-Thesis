library(cit)
# Load list and data
load("/Users/timvigers/Documents/GitHub/MS-Thesis/data/networks/pair_list.Rdata")
load("/Users/timvigers/Documents/GitHub/MS-Thesis/data/networks/pair_data.Rdata")
pair_data$T1Dgroup = factor(pair_data$T1Dgroup)
# Run CIT package
cits = apply(pairs, 1, function(x){
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  temp = pair_data[,c("T1Dgroup",methyl,metab)]
  temp = temp[complete.cases(temp),]
  y = as.numeric(temp$T1Dgroup) - 1
  X = temp[,c(methyl,metab)]
  cinf1 = cit.bp(L = X[,1],G = X[,2],T = y)
  cinf2 = cit.bp(L = X[,2],G = X[,1],T = y)
  p = min(cinf1$p_cit,cinf2$p_cit)
  d = which.min(c(cinf1$p_cit,cinf2$p_cit))
  d = ifelse(d == 1,">","<")
  if (p < 0.05) {return(paste(methyl,d,metab,p))} else {NA}
})
# Convert to dataframe
cits = as.data.frame(cits[!is.na(cits)])
cits = separate(cits,"cits[!is.na(cits)]",c("methyl","direction","metab","cit p"),
                sep = " ")
cits$`cit p` = round(as.numeric(cits$`cit p`),3)
cits = cits[order(cits$`cit p`),]
rownames(cits) = 1:nrow(cits)
# Write
save(cits,file = "./data/networks/cits.Rdata")