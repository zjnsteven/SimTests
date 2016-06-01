
# #Propensity Tree
# dbb = trans_dta@data
# crxvdata = dbb
# 
# idx = sample(sample(1:2, nrow(crxvdata), replace = TRUE))
# trainingset <- subset(crxvdata, idx %in% 1)
# testset <- subset(crxvdata, idx %in% 2)
# # half data for the propensity tree
# fit1 = rpart(treatment.status ~ modelVar + coord1 + coord2, method="class", data=crxvdata)
# #prp(fit1)
# id = c(1:dim(fit1$frame)[1])
# leafs = id[which(as.character(fit1$frame$var) == "<leaf>")]
# 
# # half for estimation
# fit = as.party(fit1)
# #Dan added spdf@data here - it's wrong, but testing.
# testset = spdf@data
# node_id =  predict(fit,testset,type="node")
# res = 0
# 
# res = rep(0,length(leafs))
# id = c(1:dim(testset)[1])
# 
# for(i in 1:length(leafs)){
#   leaf = leafs[i]
#   nodes = id[node_id==leaf]
#   df = testset[nodes,]
#   trtcount = length(df$trueOutcome[df$treatment.status==1])
#   untrtcount = length(df$trueOutcome[df$treatment.status==0])
#   treated = sum(df$trueOutcome[df$treatment.status==1])
#   untreated = sum(df$trueOutcome[df$treatment.status==0])
#   res[i] = treated/trtcount - untreated/untrtcount
# }
# 
# treatment.predictions@data$propensity.tree = unlist(lapply(as.numeric(node_id),function(x) res[match(x,leafs)]))

