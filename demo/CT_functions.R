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
  for (i in 1:k){
    trainingset <- subset(crxvdata, id %in% list[-i])
    testset <- subset(crxvdata, id %in% c(i))
    fit1 = rpart (cbind(modelOutcome,treatment.status,m1.pscore,transOutcome)  ~ modelVar + coord1 + coord2 +trans_dist,
                  trainingset,
                  control = rpart.control(cp = alpha,minsplit = tree_split_lim),
                  method=alist)
    
    
    treesize = treesize + dim(fit1$frame[which(fit1$frame$var=="<leaf>"),])[1]
    pt = predict(fit1,testset,type = "matrix")
    y = data.frame(pt)
    val = data.matrix(y)
    idx = as.numeric(rownames(testset))
    dbidx = as.numeric(rownames(crxvdata))
    for(pid in 1:(dim(testset)[1])){
      id = match(idx[pid],dbidx)
      error = error + (crxvdata$transOutcome[id] - val[pid])^2
    }
  }
  
  tsize = c(tsize,treesize/k)
  if(error == 0){
    errset[l] = 1000000
  }
  else{
    errset[l] = error/k
  }
  msg = paste(l,": ",errset[l]*k,sep="")
}

tsize = tsize[-1]
alpha_res = alphacandidate[which.min(errset)]


fit_ctpred <- rpart(cbind(modelOutcome,treatment.status,m1.pscore,transOutcome) ~ modelVar + coord1 + coord2 + trans_dist,
                    crxvdata, control=rpart.control(minsplit=tree_split_lim,cp=alpha_res),
                    method=alist)