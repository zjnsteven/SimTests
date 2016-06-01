
tot.fit.spill <- rpart(transOutcome ~ modelVar + coord1 + coord2 +trans_dist,
                       data = trans_dta@data,
                       control=rpart.control(cp=0, minsplit=tree_split_lim),
                       method="anova")
cpbest <- tot.fit.spill$cptable[which.min(tot.fit.spill$cptable[,"xerror"]),"CP"]

tot.fit.spillB <- rpart(transOutcome ~ modelVar + coord1 + coord2,
                        data = trans_dta@data,
                        control=rpart.control(cp=cpbest, minsplit=tree_split_lim),
                        method="anova")

spdf@data["coord1"] <- coordinates(spdf)[,1]
spdf@data["coord2"] <- coordinates(spdf)[,2]

res <- predict(tot.fit.spillB, newdata=spdf@data)

# Add by Jianing
db = spdf@data
tree2 = tot.fit.spillB
tree2$frame$yval = as.numeric(rownames(tree2$frame))
res2 = predict(tree2,newdata=spdf@data)

leaf = unique(tot.fit.spillB$where)
res3 = 0
for(i in 1:length(leaf)){
  
  treatcount = 0
  untreatcount = 0
  count = 0
  temp  = 0
  
  for(j in 1:nrow(db)){
    if(res2[j] == leaf[i]){
      temp = c(temp,j)
      count = count + 1
      if(db$treatment.status[j] == 1){
        treatcount = treatcount + 1
      }
      if(db$treatment.status[j] == 0){
        untreatcount = untreatcount + 1
      }
    }
    
  }
  
  if(count == treatcount | count == untreatcount ){
    temp = temp[-1]
    res3 = c(res3,temp)
  }
}

res_leaf2prune = res3[-1]