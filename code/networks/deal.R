library(bnlearn)
library(deal)
# Load data
setwd("/Users/timvigers/Documents/GitHub/MS-Thesis")
load("./data/networks/pair_data.Rdata")
load("./data/networks/cits.Rdata")
load("./data/networks/pair_list.Rdata")
# Iterate through pairs
best <- apply(pairs,1,function(x){
  # Prepare data with no missing
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  pair = pair_data[,c("T1Dgroup",methyl,metab)]
  pair = pair[complete.cases(pair),]
  pair$T1Dgroup <- factor(pair$T1Dgroup,levels = c("T1D control","T1D case"))
  pair[,2:3] <- lapply(pair[,2:3],as.numeric)
  # Deal package
  # Make an empty network
  net <- network(pair)
  # Learn parameters to get score for empty network
  learned <- learn(net,pair)$nw
  # Greedy search - based on Bayes factor
  search <- autosearch(learned,pair)
  # Calculate BIC from score (score = loglikelihood)
  search_t <- as.data.frame(search$table)
  search_t$score <- as.numeric(as.character(search_t$score))
  search_t$model <- as.character(search_t$model)
  # Get number of parameters
  arcs <- lapply(search_t$model, function(m){nrow(arcs(bnlearn::model2network(m)))})
  search_t$d <- unlist(arcs)
  d <- nrow(arcs)
  search_t$BIC <- search_t$score - (search_t$d/2) * log(nrow(pair))
  # Same for perturbation check
  check <- heuristic(search$nw, pair, restart = 24)
  check_t <- as.data.frame(search$table)
  check_t$score <- as.numeric(as.character(check_t$score))
  check_t$model <- as.character(check_t$model)
  # Get number of parameters
  arcs <- lapply(check_t$model, function(m){nrow(arcs(bnlearn::model2network(m)))})
  check_t$d <- unlist(arcs)
  d <- nrow(arcs)
  check_t$BIC <- check_t$score - (check_t$d/2) * log(nrow(pair))
  # Return best model and score
  if (all(search_t[which.max(search_t$BIC),] == check_t[which.max(check_t$BIC),])){
    r <- search_t[which.max(search_t$BIC),]
    r$model <- gsub(methyl,"methyl",r$model)
    r$model <- gsub(metab,"metab",r$model)
    return(r)
  } else {
    return(rep(NA,4))
  }
})
# Format output and write
deal_best_models <- do.call(rbind,best)
save(deal_best_models,file = "./data/networks/deal_best_models.Rdata")
