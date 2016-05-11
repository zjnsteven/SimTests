rm(list = ls())
library(sp)
library(ncf)
library(gstat)
library(devtools)
library(rpart)
library(rpart.plot)
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp("/home/aiddata/Desktop/Github/MatchIt/demo/splitc.cpp")

#detach("package:MatchIt", unload=TRUE)
#load_all("~/Desktop/Github/MatchIt/R")
library(devtools)
install_github("itpir/matchit")
library(MatchIt)




iterations <- 25

results <- data.frame(
  id=c(1:iterations)
)

results_out <- data.frame(
  id=c(1:iterations)
)


results_nospill <- data.frame(
  id=c(1:iterations)
)


ptm <- proc.time()
  for(p in 1:iterations)
{

  results["id"] <- p
  results_out["id"] <- p
  # -----------------------------------------------------------------------------
  # Simulation Settings
  # -----------------------------------------------------------------------------
  
  
  #General Options
  # set dataframe size (number of points)
  nrandom <- 2000
  
  #General psill
  #Note: setting this to 0 will result in no data randomization.
  #Larger values indicate more autocorrelation.
  psill <- 1 + runif(1, -.95, 20)
  
  # define bounding box
  minx <- -45
  maxx <- 45
  miny <- -22.5
  maxy <- 22.5
  
  #Covariate Spatial Correlation
  var1.vrange <- 1 + runif(1, -.95, 20)
  
  
  #Degree to which the covariate
  #explains the propensity to receive
  #treatment (1.0 = perfect correlation, or no error)
  prop_acc = runif(1, 0.05, .95)
  
  #Spatial pattern in the PSM error, if any 
  #(vrange = 1 approximates random noise)
  var1_error.vrange <- runif(1, 0.1, 20)
  
  #Define the spatial pattern of any model error.
  #Magnitue of error is 0-1 (0 = no error)
  mod_error.vrange <-  1.0 
  mod_error.magnitude <- runif(1, 0.1, 2)
  
  #Percent of locations which are to be defined as
  #treated
  trt_prc = runif(1, 0.45, 0.5)
  
  
  #Theta coefficient for treatment effect
  #used for defining the outcome
  theta <- 1#runif(1, -.95, 10)
  
  #Beta coefficient for ancillary data
  #in defining the treatment
  beta <- 1#runif(1,1,5) * theta
  
  #Spillover Range
  spill.vrange <- 1.0 + runif(1, -.95, 2)
  
  #Spillover Magnitude (relative to theta)
  spill.magnitude <- 1 * runif(1, 0, 3)
  
  #Caliper for Matching
  cal = .5

  sample_size = runif(1,0.05,1.0)
  
  tree_split_lim = sample(2:50,1) 
  
  # -----------------------------------------------------------------------------
  # Data Simulation
  # -----------------------------------------------------------------------------
  
  #Given your settings, outputs:
  #spdf@data$modelVar - an ancillary variable with error (defined by var1_error)
  #spdf@data$treatment.status - Binary 1/0 treated
  #spdf@data$trueOutcome - the measured outcome
  #spdf@data$modelOutcome - measured outcome including measurement error (model.error)
  source("/home/aiddata/Desktop/Github/MatchIt/demo/simulation_spatial_data.R")
  
  
  
  # -----------------------------------------------------------------------------
  # Data Generation Visualizations (Per-iteration; currently not saved)
  # Disable for a large speed boost.
  # Note visualizations use MatchIt for PSM calculations (nearest/logit/caliper=0.25)
  # -----------------------------------------------------------------------------
  #Creates a figure describing your iteration (PSM vs. observed and parameterized spillovers
  #and maps).
 # source("/home/aiddata/Desktop/Github/MatchIt/demo/test.spatial.figures.R")
  
  
  # -----------------------------------------------------------------------------
  # Model Tests
  # Save all Model Predictions into a SPDF for comparison.
  # -----------------------------------------------------------------------------
  outcome.predictions <- spdf
  outcome.predictions@data <- outcome.predictions@data[c(1,8)]
  
  treatment.predictions <- spdf
  treatment.predictions@data <- treatment.predictions@data[c(1,6,7)]
  treatment.predictions@data$trueTreatment <- (treatment.predictions@data$treatment.status)*
                                               theta + treatment.predictions$trueSpill
  
  nospill.t.pred <- spdf
  nospill.t.pred@data <- nospill.t.pred@data[1]
  nospill.t.pred@data$trueTreatment <- theta
  
  model_dta <- spdf[sample(nrow(spdf), (nrandom * sample_size)), ]
  
  #No Matching
  baseline <- lm(modelOutcome ~ treatment.status +  modelVar, data=model_dta@data)
  outcome.predictions@data$baseline <- predict(baseline, newdata=spdf@data)
  treatment.predictions@data$baseline <- summary(baseline)$coefficients[2]
  nospill.t.pred@data$baseline <- summary(baseline)$coefficients[2]
    
  #Baseline for Comparison
  baseline.matchit <- matchit(treatment.status ~ modelVar, data= model_dta@data,
                    method="nearest", distance="logit", 
                    caliper=cal, calclosest=FALSE, calrandom=FALSE)
  
  baseline.model <- lm(modelOutcome ~ treatment.status +  modelVar, 
                       data=match.data(baseline.matchit))
  
  outcome.predictions@data$baseline.matchit <- predict(baseline.model, newdata=spdf@data)
  treatment.predictions@data$baseline.matchit <- summary(baseline.model)$coefficients[2]
  nospill.t.pred@data$baseline.matchit <- summary(baseline.model)$coefficients[2]
  
  
  #Cheating Spatial PSM - we give the accurate vrange, and use it as a threshold.
  spatial.opts <- list(decay.model = "threshold",
                       threshold = spill.vrange)

  spatial.trueThreshold <- matchit(treatment.status ~ modelVar, data= model_dta,
                    method = "nearest", distance = "logit", 
                    caliper=cal, calclosest=FALSE, calrandom=FALSE,
                    spatial.options=spatial.opts)
  
  spatial.trueThreshold.model <- lm(modelOutcome ~ treatment.status +  modelVar, 
                              data=match.data(spatial.trueThreshold))
  outcome.predictions@data$spatial.trueThreshold <- predict(spatial.trueThreshold.model,
                                                            newdata=spdf@data)
  treatment.predictions@data$spatial.trueThreshold <- 
    summary(spatial.trueThreshold.model)$coefficients[2]
  
  nospill.t.pred@data$spatial.trueThreshold <- summary(spatial.trueThreshold.model)$coefficients[2]
  
  
  #PSM-approximating Traditional and Spatial PSMs

  
  #Propensity Correlogram
  p_cor_spdf <- model_dta
  p_cor_spdf$m1.pscore <- baseline.matchit$distance
  correlog.pscore.spillover <- correlog(x=p_cor_spdf@coords[, 1],
                                        y=p_cor_spdf@coords[, 2],
                                        z=p_cor_spdf$treatment.status,
                                        increment=500,
                                        latlon=TRUE, na.rm=TRUE, resamp=5,
                                        quiet=FALSE)
  
  pscore.spillover.model.dta <- data.frame(
    mean.of.class = c(correlog.pscore.spillover$mean.of.class),
    correlation = c(correlog.pscore.spillover$correlation)
  )
  
  correlog.m1.polynomial <- lm(correlation ~ poly(mean.of.class, degree=5, raw=TRUE),
                               data=pscore.spillover.model.dta)
  
  estimated.spillover.weights <- c()
  for (k in 1:length(p_cor_spdf)) {
    
    correlog.dist <- spDists(p_cor_spdf[k, ]@coords, p_cor_spdf@coords, longlat=TRUE)
    
    
    correlog.neighbors <- correlog.dist > 0
    correlog.treated <- p_cor_spdf$treatment.status
    
    
    
    correlog.newdata <- data.frame(
      mean.of.class = c(correlog.dist)
    )
    
    #Limit to the min (after the final predicted distance all values are 0)
    correlog.newdata$weight_dist <- predict(correlog.m1.polynomial, newdata=correlog.newdata)
    correlog.newdata$weight_dist[correlog.newdata$weight_dist < 0] <- 0
    
    estimated.spillover.weights[k] <- sum(correlog.newdata$weight_dist * correlog.neighbors * 
                                           correlog.treated)
    
  }
  
  p_cor_spdf$spillover.est <- estimated.spillover.weights
  
  #Spatial Matchit
  spatial.opts <- list(decay.model = "threshold",
                       threshold = correlog.pscore.spillover$x.intercept+1)#(correlog.pscore.spillover$x.intercept+1))
  #print(correlog.pscore.spillover$x.intercept)
  spatial.matchit.spill <- matchit(treatment.status ~ modelVar, data=p_cor_spdf@data,
                                    method="nearest", distance="logit", 
                                    caliper=cal, calclosest=FALSE, calrandom=FALSE)
  
  spatial.matchit.spill.model <- lm(modelOutcome ~ treatment.status +  modelVar, 
                                     data=match.data(spatial.matchit.spill))
  
  outcome.predictions@data$spatial.matchit.spill <- predict(spatial.matchit.spill.model, 
                                                          newdata=spdf@data)
  treatment.predictions@data$spatial.matchit.spill <- 
    ((summary(spatial.matchit.spill.model)$coefficients[2] * nrandom)) / 
    nrandom
  
  nospill.t.pred@data$spatial.matchit.spill <- summary(spatial.matchit.spill.model)$coefficients[2]
  
  p_cor_spdf@data["coord1"] <- coordinates(p_cor_spdf)[,1]
  p_cor_spdf@data["coord2"] <- coordinates(p_cor_spdf)[,2]
  
  #TOT - Non Random Forest
  trans_dta <- p_cor_spdf
  trans_dta <- trans_dta[(trans_dta@data$m1.pscore != 0 & 
                            trans_dta@data$m1.pscore != 1),]
  
  transOutcome <- list(rep(0,nrow(trans_dta)))
  
  for(i in 1:nrow(p_cor_spdf))
  {
    if(trans_dta$treatment.status[i] == 1)
    {
      transOutcome[i] = trans_dta@data$modelOutcome[i] / trans_dta@data$m1.pscore[i]
    }
    else
    {
      transOutcome[i] = trans_dta@data$modelOutcome[i] / (1-trans_dta@data$m1.pscore[i])
    }
  }
  trans_dta@data$transOutcome <- unlist(transOutcome)
  
  tot.fit.spill <- rpart(transOutcome ~ modelVar + coord1 + coord2,
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
  #res[res_leaf2prune] = NA
  
  ## end
  
  
  #Total Outcome
  outcome.predictions@data$tot.spill <- NA

  treatment.predictions@data$tot.spill <-  res * treatment.predictions@data$treatment.status 

  nospill.t.pred@data$tot.spill <- NA
  
  
  
  #CT
  source("/home/aiddata/Desktop/Github/MatchIt/demo/CT_functions.R")
  alist <- list(eval=ctev, split=ctsplit, init=ctinit)
  
  #trans_dta
  dbb = trans_dta@data
  k = 10 
  n = dim(dbb)[1]
  #sample_size = floor(n)
  #ridx = sample(1:n,sample_size,replace=FALSE)
  #crxvdata= dbb[ridx,]
  crxvdata = dbb
  crxvdata$id <- sample(1:k, nrow(crxvdata), replace = TRUE)
  list = 1:k
  fit1 = rpart(cbind(modelOutcome,treatment.status,m1.pscore,transOutcome) ~ modelVar + coord1 + coord2,
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
  for(j in 2:(length(alphalist)-1)){
    res[j] = sqrt(alphalist[j]*alphalist[j+1])
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
      fit1 = rpart (cbind(modelOutcome,treatment.status,m1.pscore,transOutcome)  ~ modelVar + coord1 + coord2,
                    trainingset,
                    control = rpart.control(cp = alpha,minsplit = tree_split_lim),
                    method=alist)
      
      if(dim(fit1$frame)[1] == 1){
        error = 0
        break
      }
      
      else{
        treesize = treesize + dim(fit1$frame[which(fit1$frame$var=="<leaf>"),])[1]
        pt = predict(fit1,testset,type = "matrix")
        y = data.frame(pt)
        val = data.matrix(y)
        idx = as.numeric(rownames(y))
        dbidx = as.numeric(rownames(dbb))
        
        for(pid in 1:(dim(y)[1])){
          id = match(idx[pid],dbidx)
          error = error + (crxvdata$transOutcome[id] - val[pid])^2
          #print(error)
        }
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
    #print(msg)
  }
  
  tsize = tsize[-1]
  alpha_res = alphacandidate[which.min(errset)]
  fit_ctpred <- rpart(cbind(modelOutcome,treatment.status,m1.pscore,transOutcome) ~ modelVar + coord1 + coord2,
                crxvdata, control=rpart.control(minsplit=tree_split_lim,cp=alpha_res),
                method=alist)
  #prp(fit_ctpred)
  #res = rep(0,length(alphalist)-1)
  #for(j in 2:(length(alphalist)-1)){
  #  res[j] = sqrt(alphalist[j]*alphalist[j+1])
  #}
  
  #Total Outcome - CT
  outcome.predictions@data$ct.spill <- NA
  
  treatment.predictions@data$ct.spill <-  predict(fit_ctpred,newdata=spdf@data) * treatment.predictions@data$treatment.status 
  print("CT Nodes:")
  print(length(unique(fit_ctpred$where)))
  
  nospill.t.pred@data$ct.spill <- NA
  
  #Save Summary Results
  
  for(i in 3:length(nospill.t.pred@data))
  {

    if(p == 1)
    {
      results_nospill[names(nospill.t.pred@data)[i]] <- NA
      results_nospill["trueTreatment"] <- NA
    }
    results_nospill[names(nospill.t.pred@data)[i]][p,] <- sum(nospill.t.pred@data[,i]) / (nrandom*trt_prc)
    results_nospill["trueTreatment"][p,] <- theta
  }
  
  for(i in 4:length(treatment.predictions@data))
  {
    if(p == 1)
    {
      results[names(treatment.predictions@data)[i]] <- NA
      results["trueTreatment"] <- NA
    }
    results[names(treatment.predictions@data)[i]][p,] <- sum(treatment.predictions@data[,i],na.rm=T) / (nrandom*trt_prc)
    results["trueTreatment"][p,] <- sum(treatment.predictions@data$trueTreatment) / (nrandom*trt_prc)
  }

  for(i in 3:length(outcome.predictions@data))
  {
    if(p == 1)
    {
      results_out[names(outcome.predictions@data)[i]] <- NA
      results_out["trueTreatment"] <- NA
    }
    results_out[names(outcome.predictions@data)[i]][p,] <- sum(outcome.predictions@data[,i]) / (nrandom*trt_prc)
    results_out["trueTreatment"][p,] <- sum(treatment.predictions@data$trueTreatment) / (nrandom*trt_prc)
    }

  #Save relevant parameters
  if(p == 1)
  {
    results["spill.magnitude"] <- NA
    results["psill"] <- NA
    results["var1.vrange"] <- NA
    results["prop_acc"] <- NA
    results["mod_error.vrange"] <- NA
    results["mod_error.magnitude"] <- NA
    results["spill.vrange"] <- NA
    results["beta"] <- NA
    results["var1_error.vrange"] <- NA
    results["theta"] <- NA
    results["caliper"] <- NA
    results["sample_size"] <- NA
    results["tree_split_lim"] <- NA
    
    
    results_out["spill.magnitude"] <- NA
    results_out["psill"] <- NA
    results_out["var1.vrange"] <- NA
    results_out["prop_acc"] <- NA
    results_out["mod_error.vrange"] <- NA
    results_out["mod_error.magnitude"] <- NA
    results_out["spill.vrange"] <- NA
    results_out["beta"] <- NA
    results_out["var1_error.vrange"] <- NA
    results_out["theta"] <- NA
    results_out["caliper"] <- NA
    results_out["sample_size"] <- NA
    results_out["tree_split_lim"] <- NA
    
    results_nospill["spill.magnitude"] <- NA
    results_nospill["psill"] <- NA
    results_nospill["var1.vrange"] <- NA
    results_nospill["prop_acc"] <- NA
    results_nospill["mod_error.vrange"] <- NA
    results_nospill["mod_error.magnitude"] <- NA
    results_nospill["spill.vrange"] <- NA
    results_nospill["beta"] <- NA
    results_nospill["var1_error.vrange"] <- NA
    results_nospill["theta"] <- NA
    results_nospill["caliper"] <- NA
    results_nospill["sample_size"] <- NA
    results_nospill["tree_split_lim"] <- NA
    
    
  }
  results["spill.magnitude"][p,] <- spill.magnitude
  results["psill"][p,] <- psill
  results["var1.vrange"][p,] <- var1.vrange
  results["prop_acc"][p,] <- prop_acc
  results["mod_error.vrange"][p,] <- mod_error.vrange
  results["mod_error.magnitude"][p,] <- mod_error.magnitude
  results["spill.vrange"][p,] <- spill.vrange
  results["beta"][p,] <- beta
  results["theta"][p,] <- theta
  results["var1_error.vrange"][p,] <- var1_error.vrange
  results["caliper"][p,] <- cal
  results["sample_size"][p,] <- sample_size
  results["tree_split_lim"][p,] <- tree_split_lim
  
  results_out["spill.magnitude"][p,] <- spill.magnitude
  results_out["psill"][p,] <- psill
  results_out["var1.vrange"][p,] <- var1.vrange
  results_out["prop_acc"][p,] <- prop_acc
  results_out["mod_error.vrange"][p,] <- mod_error.vrange
  results_out["mod_error.magnitude"][p,] <- mod_error.magnitude
  results_out["spill.vrange"][p,] <- spill.vrange
  results_out["beta"][p,] <- beta
  results_out["theta"][p,] <- theta
  results_out["var1_error.vrange"][p,] <- var1_error.vrange
  results_out["caliper"][p,] <- cal
  results_out["sample_size"][p,] <- sample_size
  results_out["tree_split_lim"][p,] <- tree_split_lim
  
  results_nospill["spill.magnitude"][p,] <- spill.magnitude
  results_nospill["psill"][p,] <- psill
  results_nospill["var1.vrange"][p,] <- var1.vrange
  results_nospill["prop_acc"][p,] <- prop_acc
  results_nospill["mod_error.vrange"][p,] <- mod_error.vrange
  results_nospill["mod_error.magnitude"][p,] <- mod_error.magnitude
  results_nospill["spill.vrange"][p,] <- spill.vrange
  results_nospill["beta"][p,] <- beta
  results_nospill["theta"][p,] <- theta
  results_nospill["var1_error.vrange"][p,] <- var1_error.vrange
  results_nospill["caliper"][p,] <- cal
  results_nospill["sample_size"][p,] <- sample_size
  results_nospill["tree_split_lim"][p,] <- tree_split_lim

  
#Compare Maps
#spplot(outcome.predictions, zcol=names(outcome.predictions)
#       [names(outcome.predictions) != "id"])
#spplot(treatment.predictions, zcol=names(treatment.predictions)
#       [names(treatment.predictions) != "id"])
  print("Iteration Complete")
  print(p)
  print("----------------")
}

print(proc.time() - ptm)



