ctev <- function(y, wt,parms) {
  out = node_evaluate(y)
  list(label= out[1], deviance=out[2])
}

ctsplit <- function(y, wt, x, parms, continuous) {
  if (continuous) {
    n = nrow(y)
    res = splitc(y)
    list(goodness=res[1:(n-1)], direction=res[n:(2*(n-1))])
  }
  else{
    res = splitnc(y,x)
    n=(length(res)+1)/2
    list(goodness=res[1:(n-1)], direction=res[n:(2*n-1)])
  }
}


ctinit <- function(y, offset, parms, wt) {
  sfun <- function(yval, dev, wt, ylevel, digits ) {

    paste("events=", round(yval[,1]),
          ", coef= ", format(signif(yval[,2], digits)),
          ", deviance=" , format(signif(dev, digits)),
          sep = '')}
  environment(sfun) <- .GlobalEnv
  list(y =y, parms = 0, numresp = 1, numy = 4,
       summary = sfun)
}

alist <- list(eval=ctev, split=ctsplit, init=ctinit)
dbb = trans_dta@data
k = 10
n = dim(dbb)[1]
crxvdata = dbb
crxvdata$id <- sample(1:k, nrow(crxvdata), replace = TRUE)
list = 1:k
fit1 = rpart(cbind(modelOutcome,treatment.status,m1.pscore,transOutcome) ~ modelVar + coord1 + coord2 + trans_dist,
             crxvdata,
             control = rpart.control(cp = 0,minsplit = tree_split_lim),
             method=alist)
fit = data.matrix(fit1$frame)
index = as.numeric(rownames(fit1$frame))
tsize = dim(fit1$frame[which(fit1$frame$var=="<leaf>"),])[1]

alpha = 0
alphalist = 0
alphalist = cross_validate(fit, index,alphalist)
res = rep(0,length(alphalist)-1)
if(length(alphalist) <= 2){
  res = alphalist
}else{
  for(j in 2:(length(alphalist)-1)){
    res[j] = sqrt(alphalist[j]*alphalist[j+1])
  }
}

alphacandidate = res
alphaset = rep(0,length(alphacandidate))
errset = rep(0,length(alphacandidate))
tsize = 0

for(l in 1:length(alphacandidate)){
  alpha = alphacandidate[l]
  error = 0
  treesize = 0
  loginfo("alpha index %d value %f", l, alpha, logger=paste("simtest.", iteration, ".", "CT", sep="") )
  for (i in 1:k){
    trainingset <- subset(crxvdata, id %in% list[-i])
    testset <- subset(crxvdata, id %in% c(i))
    fit1 = rpart (cbind(modelOutcome,treatment.status,m1.pscore,transOutcome)  ~ modelVar + coord1 + coord2 +trans_dist,
                  trainingset,
                  control = rpart.control(cp = alpha,minsplit = tree_split_lim),
                  method=alist)
    temperr = 0
    treesize = treesize + dim(fit1$frame[which(fit1$frame$var=="<leaf>"),])[1]
    pt = predict(fit1,testset,type = "matrix")
    y = data.frame(pt)
    val = data.matrix(y)
    idx = as.numeric(rownames(testset))
    dbidx = as.numeric(rownames(crxvdata))
    for(pid in 1:(dim(testset)[1])){
      id = match(idx[pid],dbidx) 
      if(any(is.na(id))){
        logerror("id has na", logger=paste("simtest.", iteration, ".", "CT", sep=""))
      }     
      if(any(is.na(pid))){
        logerror("pid has na", logger=paste("simtest.", iteration, ".", "CT", sep=""))
      } 
      temperr = temperr + (crxvdata$transOutcome[id] - val[pid])^2
      
    }
    if(is.numeric(temperr)){
      loginfo("alpha id %d fold %d error %f", l, i ,temperr, logger=paste("simtest.", iteration, ".", "CT", sep=""))
      error = error + temperr
    }
    else{
      logerror("temperr is not numeric", logger=paste("simtest.", iteration, ".", "CT", sep=""))
    }
    
  }
  
  tsize = c(tsize,treesize/k)
  if(error == 0){
    errset[l] = 1000000
  }
  else{
    errset[l] = error/k
  }
  loginfo("alpha id %d  error avg %f", l, errset[l], logger=paste("simtest.", iteration, ".", "CT", sep=""))
}

tsize = tsize[-1]
alpha_res = alphacandidate[which.min(errset)]

loginfo("best alpha %f", alpha_res, logger=paste("simtest.", iteration, ".", "CT", sep=""))
fit_ctpred <- rpart(cbind(modelOutcome,treatment.status,m1.pscore,transOutcome) ~ modelVar + coord1 + coord2,
                    crxvdata, control=rpart.control(minsplit=tree_split_lim,cp=alpha_res),
                    method=alist)

treesummary = table(fit_ctpred$frame$var)

if(length(treesummary) == 1){
  logwarn("no split", logger=paste("simtest.", iteration, ".", "CT", sep=""))
}

for(i in 1:length(treesummary)){
  loginfo("split covariates: %s, number: %d", rownames(treesummary)[i],treesummary[i] ,logger=paste("simtest.", iteration, ".", "CT", sep=""))
}
