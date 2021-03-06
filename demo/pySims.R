rm(list=ls())
library(sp)
library(ncf)
library(gstat)
library(devtools)
library(rpart)
library(rpart.plot)
library(Rcpp)
library(classInt)
library(RColorBrewer)
library(methods)
library(partykit)
library(logging)
basicConfig()


path_base = "/sciclone/home00/geogdan/SimTests/"
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")



sourceCpp(paste(path_base, "demo/splitc.cpp", sep=""))
CT_src <- paste(path_base, "demo/CT_functions.R", sep="")
sim_src <- paste(path_base, "demo/simulation_spatial_data.R", sep="")
map_out <- "/sciclone/home00/geogdan/jM/"


#detach("package:MatchIt", unload=TRUE)
#load_all("/sciclone/home00/geogdan/SimTests/R")
load_all(paste(path_base,"R",sep=""))

#1 3800.41856568 0.90376081592 -45.0 45.0 -22.5 22.5 3.21749654825 0.250852506018 0.448021052911 4.27592030555 0.0684864449219 0.29100048171 1 0.330411927736 3.83573033709 1.88067542642 0.698254286741 0.437623061042 10 2.58494466138 /sciclone/home00/geogdan/AlphaSims/test_0.csv 0.954552979835 0.539550663469 0.164665770447
Args <- commandArgs(trailingOnly = TRUE)

if(length(Args) == 0)
{
  #Exception for running the script directly without a call.
Args <- c("1", "3800.41856568", "0.90376081592", "-45.0", "45.0", "-22.5", "22.5", "50000", "0.250852506018", "0.448021052911", "4.27592030555", "0.0684864449219", "0.29100048171", "1", "0.330411927736", "50000", "1.88067542642", "0.698254286741", "0.437623061042", "10", "2.58494466138", "/sciclone/home00/geogdan/may_a/test_0.csv", "0.954552979835", "0.539550663469", "50000")
}

print(Args)
out_path=Args[22]
out_itDta_path = paste(substr(out_path, 1, nchar(out_path)-4),"_dta.csv",sep="")
nums = as.numeric(Args)

addHandler(writeToFile, logger=paste("simtest.",version,sep=""), file="/sciclone/home00/geogdan/logs/simtest.log")


# -----------------------------------------------------------------------------
# Simulation Settings
# -----------------------------------------------------------------------------

version = Args[1]
nrandom = round(as.numeric(Args[2]), 0)
xvar_error_psill = as.numeric(Args[3])
minx = as.numeric(Args[4])
maxx = as.numeric(Args[5])
miny = as.numeric(Args[6])
maxy = as.numeric(Args[7])

var1.vrange = as.numeric(Args[8])
var1.error = as.numeric(Args[9])
prop_acc = as.numeric(Args[10])
var1_error.vrange = as.numeric(Args[11])
mod_error.magnitude = as.numeric(Args[12])
trt_prc = as.numeric(Args[13])
theta = as.numeric(Args[14])
beta =as.numeric( Args[15])
spill.vrange = as.numeric(Args[16])
spill.magnitude= as.numeric(Args[17])
cal= as.numeric(Args[18])
sample_size = as.numeric(Args[19])
tree_split_lim = 5 
mod_error.vrange= as.numeric(Args[21])
xvar_psill=as.numeric(Args[23])
mod_error_psill=as.numeric(Args[24])
trt_spill_sill=as.numeric(Args[25])
p <- 1

loginfo("maximum amount of spatial covariation for covariates %f",xvar_error_psill, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("latitude and longitude bound %f, %f, %f, %f",minx, maxx, miny, maxy , logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )

loginfo("distance threshold for covariates %f",var1_vrange, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("total error in the covarite  %f",var1_error, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("prop_acc  %f",prop_acc, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("the range of error for the covariate  %f",var1_error_vrange, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("error coefficient  %f",mod_error_magnitude, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("the pecentage of treated points  %f",trt_prc, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("treatment coefficient  %f",theta, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )

loginfo("covariate coefficient  %f",beta, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("distance of treatment effect %f",spill_vrange, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("coefficient of spill_vrange %f",spill_magnitude, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("distance of cov error %f",mod_error_vrange, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("maximum amount of spatial covariation for covariates for errors %f",xvar_error_psill, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("maximum amount of saptial cov error  %f",mod_error_psill, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("magnitude of treatment spillover  %f",trt_spill_sill, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )

#print("Total Size / Sample Size / Tree Split Limit:")
#print(nrandom)
#print(sample_size)
#print(tree_split_lim)
#loginfo("Total Size %d / Sample Size %d / Tree Split Limit %d :", nrandom, sample_size, tree_split_lim, logger="")


per_split_lim <- tree_split_lim / sample_size


iterations <- 1
iteration = iterations
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


loginfo("Total Size %d",nrandom, logger=paste("simtest.", iteration, ".", "Model.parameter", sep="") )
loginfo("Sample Size %f",sample_size, logger=paste("simtest.", iteration, ".", "Model.parameter", sep=""))
loginfo("Tree Split Limit %d",tree_split_lim,paste("simtest.", iteration, ".", "Model.parameter", sep="") )


# -----------------------------------------------------------------------------
# Data Simulation
# -----------------------------------------------------------------------------

#Given your settings, outputs:
#spdf@data$modelVar - an ancillary variable with error (defined by var1_error)
#spdf@data$treatment.status - Binary 1/0 treated
#spdf@data$trueOutcome - the measured outcome
#spdf@data$modelOutcome - measured outcome including measurement error (model.error)
source(sim_src)



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
outcome.predictions@data <- outcome.predictions@data[c(1,9)]

treatment.predictions <- spdf
treatment.predictions@data <- treatment.predictions@data[c(1,6,8)]
treatment.predictions@data$trueTreatment <- (treatment.predictions$treatment.status * theta) + treatment.predictions$trueSpill


nospill.t.pred <- spdf
nospill.t.pred@data <- nospill.t.pred@data[1]
nospill.t.pred@data$trueTreatment <- theta

model_dta <- spdf[sample(nrow(spdf), (nrandom * sample_size)), ]

#No Matching
baseline <- lm(modelOutcome ~ treatment.status +  modelVar, data=model_dta@data)
outcome.predictions@data$baseline <- predict(baseline, newdata=spdf@data)
treatment.predictions@data$baseline <- summary(baseline)$coefficients[2] * treatment.predictions$treatment.status
nospill.t.pred@data$baseline <- summary(baseline)$coefficients[2]



#Baseline for Comparison
baseline.matchit <- matchit(treatment.status ~ modelVar, data= model_dta@data,
                            method="nearest", distance="logit",
                            caliper=cal, calclosest=FALSE, calrandom=FALSE)

baseline.model <- lm(modelOutcome ~ treatment.status +  modelVar,
                     data=match.data(baseline.matchit))

outcome.predictions@data$baseline.matchit <- predict(baseline.model, newdata=spdf@data)
treatment.predictions@data$baseline.matchit <- summary(baseline.model)$coefficients[2] * treatment.predictions$treatment.status
nospill.t.pred@data$baseline.matchit <- summary(baseline.model)$coefficients[2]


#Cheating Spatial PSM - we give the accurate vrange, and use it as a threshold.
spatial.opts <- list(decay.model = "threshold",
                     threshold = (spill.vrange/111.32))

spatial.trueThreshold <- matchit(treatment.status ~ modelVar, data= model_dta,
                                 method = "nearest", distance = "logit",
                                 caliper=cal, calclosest=FALSE, calrandom=FALSE,
                                 spatial.options=spatial.opts)

spatial.trueThreshold.model <- lm(modelOutcome ~ treatment.status +  modelVar,
                                  data=match.data(spatial.trueThreshold))
outcome.predictions@data$spatial.trueThreshold <- predict(spatial.trueThreshold.model,
                                                          newdata=spdf@data)
treatment.predictions@data$spatial.trueThreshold <- summary(spatial.trueThreshold.model)$coefficients[2] * treatment.predictions$treatment.status

nospill.t.pred@data$spatial.trueThreshold <- summary(spatial.trueThreshold.model)$coefficients[2]


#PSM-approximating Traditional and Spatial PSMs


#Propensity Correlogram
p_cor_spdf <- model_dta
p_cor_spdf$m1.pscore <- baseline.matchit$distance
correlog.pscore.spillover <- correlog(x=p_cor_spdf@coords[, 1],
                                      y=p_cor_spdf@coords[, 2],
                                      z=p_cor_spdf$treatment.status,
                                      increment=500,
                                      latlon=TRUE, na.rm=TRUE, resamp=2,
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
  summary(spatial.matchit.spill.model)$coefficients[2] * treatment.predictions$treatment.status

nospill.t.pred@data$spatial.matchit.spill <- summary(spatial.matchit.spill.model)$coefficients[2]

p_cor_spdf@data["coord1"] <- coordinates(p_cor_spdf)[,1]
p_cor_spdf@data["coord2"] <- coordinates(p_cor_spdf)[,2]

#TOT - Non Random Forest

#Simple P-Score for Testing
#p.score <- lm(treatment.status ~ 0 + modelVar, data=model_dta@data)
#p_cor_spdf@data$m1.pscore <- predict(p.score, newdata=p_cor_spdf@data)



trans_dta <- p_cor_spdf

#Create a pscore using the rpart
pscore.Calc <- matchit(treatment.status ~ modelVar + coord1 + coord2 + trans_dist, data= trans_dta@data,
                            method="nearest", distance="logit")
                            

trans_dta@data$m1.pscore = pscore.Calc$distance 
#Simple caliper for testing
upper_lim <- (sd(p_cor_spdf@data$m1.pscore, na.rm=TRUE) * cal) + mean(p_cor_spdf@data$m1.pscore)
lower_lim <-  mean(p_cor_spdf@data$m1.pscore) - (sd(p_cor_spdf@data$m1.pscore, na.rm=TRUE) * cal) 


#trans_dta <- trans_dta[trans_dta@data$m1.pscore < upper_lim,]
#trans_dta <- trans_dta[trans_dta@data$m1.pscore > lower_lim,]

trans_dta <- trans_dta[(trans_dta@data$m1.pscore != 0 &
                          trans_dta@data$m1.pscore != 1),]



transOutcome <- list(rep(0,nrow(trans_dta)))

for(i in 1:nrow(trans_dta))
{
  if(trans_dta$treatment.status[i] == 1)
  {
    #Treated
    transOutcome[i] = trans_dta@data$modelOutcome[i] * 
      (1 / trans_dta@data$m1.pscore[i])
  }
  else
  {
    #Untreated
    transOutcome[i] = -1 * (trans_dta@data$modelOutcome[i] * 
      ((1-0) / (1 - trans_dta@data$m1.pscore[i])))
  }
}
trans_dta@data$transOutcome <- unlist(transOutcome)

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
#res[res_leaf2prune] = NA

## end


#Total Outcome


treatment.predictions@data$tot.spill <-  res * treatment.predictions$treatment.status





#CT
source(CT_src)
alist <- list(eval=ctev, split=ctsplit, init=ctinit)

#trans_dta
dbb = trans_dta@data
#write.csv(trans_dta@data, out_itDta_path)
k = 10
n = dim(dbb)[1]
#sample_size = floor(n)
#ridx = sample(1:n,sample_size,replace=FALSE)
#crxvdata= dbb[ridx,]
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
#if(alphalist[1]==0 & alphalist[2]==0){
#  alphalist = alphalist[-1]
#}
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
print("AlphaCandidate:")
print(alphacandidate)
for(l in 1:length(alphacandidate)){
  alpha = alphacandidate[l]
  error = 0
  treesize = 0
  #print("--")
  #print("Alpha (alphacandidate[l]:)")
  #print(l)
  #print(alpha)
  #print("--")
  loginfo("alpha index %d value %f", l, alpha, logger=paste("simtest.", iteration, ".", "CT", sep="") )
  for (i in 1:k){
    trainingset <- subset(crxvdata, id %in% list[-i])
    testset <- subset(crxvdata, id %in% c(i))
    fit1 = rpart (cbind(modelOutcome,treatment.status,m1.pscore,transOutcome)  ~ modelVar + coord1 + coord2 +trans_dist,
                  trainingset,
                  control = rpart.control(cp = alpha,minsplit = tree_split_lim),
                  method=alist)
    
    #if(dim(fit1$frame)[1] == 1){
    #    error = 0
    #    break
    #  }
    
    #  else{
    temperr = 0
    treesize = treesize + dim(fit1$frame[which(fit1$frame$var=="<leaf>"),])[1]
    pt = predict(fit1,testset,type = "matrix")
    y = data.frame(pt)
    val = data.matrix(y)
    idx = as.numeric(rownames(testset))
    dbidx = as.numeric(rownames(crxvdata))
    #print(dim(y)[1])
    #print(paste("error:", error))
    for(pid in 1:(dim(testset)[1])){
      id = match(idx[pid],dbidx)
      #print(paste("error before:", error,"val",val[pid],"tranoutcome",crxvdata$transOutcome[id]))
      
      temperr = temperr + (crxvdata$transOutcome[id] - val[pid])^2
      
    }
    loginfo("alpha id %d fold %d error %f", l, i ,temperr, logger=paste("simtest.", iteration, ".", "CT", sep=""))
    error = error + temperr
    #}
    # print(paste("error after:", error))
  }
  
  tsize = c(tsize,treesize/k)
  if(error == 0){
    errset[l] = 1000000
  }
  else{
    errset[l] = error/k
  }
  #msg = paste(l,": ",errset[l]*k,sep="")
  #print(msg)
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



#Total Outcome - CT

treatment.predictions@data$ct.spill <-  predict(fit_ctpred,newdata=spdf@data) * treatment.predictions$treatment.status
#print("CT Nodes:")
#print(length(unique(fit_ctpred$where)))




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






#Save Summary Results
for(i in 4:length(treatment.predictions@data))
{
  if(p == 1)
  {
    results[names(treatment.predictions@data)[i]] <- NA
    results["trueTreatment"] <- NA
  }
  results[names(treatment.predictions@data)[i]][p,] <- sum(abs(treatment.predictions@data[,i] * 
                                                                 treatment.predictions$treatment.status) - 
                                                             (treatment.predictions@data$trueTreatment * 
                                                                treatment.predictions$treatment.status))
  results["trueTreatment"][p,] <- sum(treatment.predictions@data$trueTreatment * treatment.predictions$treatment.status) / sum(treatment.predictions$treatment.status)
}






results["spill.magnitude"] <- NA
results["xvar_error_psill"] <- NA
results["mod_error_psill"] <- NA
results["trt_spill_sill"] <- NA
results["xvar_psill"] <- NA
results["var1.vrange"] <- NA
results["prop_acc"] <- NA
results["mod_error.vrange"] <- NA
results["mod_error.magnitude"] <- NA
results["spill.vrange"]<- NA
results["beta"] <- NA
results["theta"] <- NA
results["var1_error.vrange"] <- NA
results["caliper"] <- NA
results["sample_size"]<- NA
results["tree_split_lim"] <- NA
results["nrandom"] <- NA
results["ct_split_count"] <- NA

results["spill.magnitude"][p,] <- spill.magnitude
results["xvar_error_psill"][p,]  <- xvar_error_psill
results["mod_error_psill"][p,]  <- mod_error_psill
results["trt_spill_sill"][p,]  <- trt_spill_sill
results["xvar_psill"][p,]  <- xvar_psill
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
results["tree_split_lim"][p,] <- per_split_lim
results["nrandom"][p,] <- nrandom
results["ct_split_count"][p,] <- length(unique(fit_ctpred$where))




#Save Trees
#png()
#prp(fit_ctpred)
#prp(tot.fit.spillB)


#Compare Maps
 map_trt <- treatment.predictions[names(treatment.predictions) != "tot.spill"]
 map_trt <- map_trt[names(map_trt) != "spatial.trueThreshold"]
 map_trt <- map_trt[names(map_trt) != "spatial.matchit.spill"]
 map_trt <- map_trt[names(map_trt) != "baseline"]
 map_trt <- map_trt[names(map_trt) != "baseline.matchit"]
 #map_trt <- map_trt[names(map_trt) != "treatment.status"]
 #map_trt <- map_trt[names(map_trt) != "ct.spill"]
 #map_trt <- map_trt[names(map_trt) != "trueSpill"]
 #map_trt <- map_trt[names(map_trt) != "id"]
# 
# 
 map_trt@data$trueTreatment <- map_trt@data$trueTreatment * map_trt@data$treatment.status
 map_trt@data$ct <- map_trt@data$ct.spill * map_trt@data$treatment.status
names(map_trt)
 
map_trt <- map_trt[names(map_trt) != "ct.spill"]

 pal = brewer.pal(9,"Greens")
 brks = c(0.0,1.0,1.25,1.5,1.75,2.0,2.5,100)
map_out_path <- paste(map_out, "vsm",version, "_", spill.magnitude, ".png", sep="")
title.v <- paste("Mag:",spill.magnitude," Treat Range:",spill.vrange," RunID:", version, sep="")
png(filename=map_out_path)
print(spplot(map_trt, zcol=names(map_trt)[names(map_trt) != "id"],cuts=brks,col.regions=pal, col="transparent",
       main = list(label=title.v), cex=0.5))

 dev.off()

 
#print("Iteration Complete")


#print(proc.time() - ptm)

write.csv(results,file=out_path)
write.csv(trans_dta@data,file=out_itDta_path)


