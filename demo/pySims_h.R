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


Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp("/sciclone/home00/geogdan/SimTests/demo/splitc.cpp")
CT_src <- "/sciclone/home00/geogdan/SimTests/demo/CT_functions.R"
TOT_src <- "/sciclone/home00/geogdan/SimTests/demo/TOT_functions.R"
sim_src <- "/sciclone/home00/geogdan/SimTests/demo/simulation_spatial_data.R"
map_out <- "/sciclone/home00/geogdan/M4/"

#detach("package:MatchIt", unload=TRUE)
load_all("/sciclone/home00/geogdan/SimTests/R")



#1 3800.41856568 0.90376081592 -45.0 45.0 -22.5 22.5 3.21749654825 0.250852506018 0.448021052911 4.27592030555 0.0684864449219 0.29100048171 1 0.330411927736 3.83573033709 1.88067542642 0.698254286741 0.437623061042 10 2.58494466138 /sciclone/home00/geogdan/AlphaSims/test_0.csv 0.954552979835 0.539550663469 0.164665770447
Args <- commandArgs(trailingOnly = TRUE)

if(length(Args) == 0)
{
  #Exception for running the script directly without a call.
Args <- c("1", "3800.41856568", "0.90376081592", "-45.0", "45.0", "-22.5", "22.5", "750", "0.250852506018", "0.448021052911", "4.27592030555", "0.0684864449219", "0.29100048171", "1", "0.330411927736", "750", "1.88067542642", "0.698254286741", "0.437623061042", "10", "2.58494466138", "/sciclone/home00/geogdan/may_a/test_0.csv", "0.954552979835", "0.539550663469", "50000")
}


out_path=Args[22]
out_itDta_path = paste(substr(out_path, 1, nchar(out_path)-4),"_dta.csv",sep="")
nums = as.numeric(Args)



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


per_split_lim <- tree_split_lim / (sample_size * nrandom)

results <- data.frame(
  id=c(1:1)
)



ptm <- proc.time()

# -----------------------------------------------------------------------------
# Data Simulation
# -----------------------------------------------------------------------------
source(sim_src)

# -----------------------------------------------------------------------------
# Calculate true values, subsample data for this iteration.
# -----------------------------------------------------------------------------

treatment.predictions <- spdf
treatment.predictions@data <- treatment.predictions@data[c(1,6,8)]
treatment.predictions@data$trueTreatment <- 
  treatment.predictions$treatment.status *
  ((treatment.predictions$treatment.status * theta) + 
     treatment.predictions$trueSpill)
model_dta <- spdf[sample(nrow(spdf), (nrandom * sample_size)), ]

# -----------------------------------------------------------------------------
# Basic Linear Model Estimate
# -----------------------------------------------------------------------------
baseline <- lm(modelOutcome ~ treatment.status +  modelVar, 
               data=model_dta@data)
treatment.predictions@data$baseline.lm <- summary(baseline)$coefficients[2] * 
  treatment.predictions$treatment.status

# -----------------------------------------------------------------------------
# Basic PSM Matching Estimate
# -----------------------------------------------------------------------------
baseline.matchit <- matchit(treatment.status ~ modelVar, data= model_dta@data,
                            method="nearest", distance="logit",
                            caliper=cal, calclosest=FALSE, calrandom=FALSE)

baseline.model <- lm(modelOutcome ~ treatment.status +  modelVar,
                     data=match.data(baseline.matchit))

treatment.predictions@data$baseline.matchit <- 
  summary(baseline.model)$coefficients[2] * treatment.predictions$treatment.status



# -----------------------------------------------------------------------------
# Spatial PSM Matching Estimate
# Note this is a best case for Spatial PSM, as we give the true vrange.
# -----------------------------------------------------------------------------
spatial.opts <- list(decay.model = "threshold",
                     threshold = (spill.vrange/111.32/1000))

spatial.trueThreshold <- matchit(treatment.status ~ modelVar, data= model_dta,
                                 method = "nearest", distance = "logit",
                                 caliper=cal, calclosest=FALSE, calrandom=FALSE,
                                 spatial.options=spatial.opts)

spatial.trueThreshold.model <- lm(modelOutcome ~ treatment.status +  modelVar,
                                  data=match.data(spatial.trueThreshold))

treatment.predictions@data$spatialPSM <- 
  summary(spatial.trueThreshold.model)$coefficients[2] * 
  treatment.predictions$treatment.status

# -----------------------------------------------------------------------------
# GWR with MatchIt
#Note this is a best case for GWR, as we give it the true range.
# -----------------------------------------------------------------------------
treatment.predictions@data$gwr <- NA
sum_trt <- 0
cnt_trt <- 0
for (i in 1:length(treatment.predictions))
{
  it_dfa <- model_dta
  if(treatment.predictions@data$treatment.status[i] == 1)
  {
  it_dfa@data$dists <- as.vector(spDists(model_dta@coords, treatment.predictions[i, ]@coords, longlat=TRUE))
  it_dfa <- it_dfa[it_dfa@data$dists <= spill.vrange,]
  
  mod_it_dfa <- it_dfa[,names(it_dfa@data) %in% c("treatment.status", "modelVar", "modelOutcome")]
  

  result = tryCatch(
    gwr.matchit <- matchit(treatment.status ~ modelVar, data= mod_it_dfa@data,
                              method="nearest", distance="logit",
                              caliper=cal, calclosest=FALSE, calrandom=FALSE), 
    error = function(e) {
      print(e) 
    return("Error")})


  if(result == "Error")
  {
    treatment.predictions@data$gwr[[i]] <- NA
    result = "Reset"
  }
  
  else
  {
  gwr.model <- lm(modelOutcome ~ treatment.status +  modelVar,
                       data=match.data(gwr.matchit))
  
  treatment.predictions@data$gwr[[i]]<- 
    summary(gwr.model)$coefficients[2]
  
  sum_trt <- sum_trt + summary(gwr.model)$coefficients[2]
  cnt_trt <- cnt_trt + 1
  }
  }
  else
  {
    treatment.predictions@data$gwr[[i]] <- 0
  }
  
}


for (i in 1:length(treatment.predictions))
{
  if(is.na(treatment.predictions@data$gwr[[i]]))
  {
  treatment.predictions@data$gwr[[i]] <- sum_trt / cnt_trt
  }
}

# -----------------------------------------------------------------------------
# Create coordinates for the trees to split along.
# -----------------------------------------------------------------------------
trans_dta <- model_dta
trans_dta@data["coord1"] <- coordinates(model_dta)[,1]
trans_dta@data["coord2"] <- coordinates(model_dta)[,2]

# -----------------------------------------------------------------------------
# Calculate the Pscore for the Trees, remove 0 and 1 cases.
# -----------------------------------------------------------------------------
pscore.Calc <- matchit(treatment.status ~ modelVar + coord1 + coord2, data= trans_dta@data,
                            method="nearest", distance="logit")
                            
trans_dta@data$m1.pscore = pscore.Calc$distance 

trans_dta <- trans_dta[(trans_dta@data$m1.pscore != 0 &
                          trans_dta@data$m1.pscore != 1),]


# -----------------------------------------------------------------------------
# Calculate the transformed outcome for the trees.
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Transformed Outcome Tree
# -----------------------------------------------------------------------------
source(TOT_src)
treatment.predictions@data$tot.pred <-  res * treatment.predictions$treatment.status

# -----------------------------------------------------------------------------
# Causal Tree
# -----------------------------------------------------------------------------
source(CT_src)
treatment.predictions@data$ct.pred <-  predict(fit_ctpred,newdata=spdf@data) * 
  treatment.predictions$treatment.status







# -----------------------------------------------------------------------------
# Save Summary Results for this Iteration
# Quantity Error: Difference in total treated units estimate
# Allocation error: The spatial mis-match in estimated locations of treated and
# untreated
# Total Error: Sum of Quantity and Allocation error
# -----------------------------------------------------------------------------
for(i in 4:length(treatment.predictions@data))
{
  q_name <- paste(names(treatment.predictions@data)[i],"_quant", sep="")
  a_name <- paste(names(treatment.predictions@data)[i],"_tot", sep="")
  t_name <- paste(names(treatment.predictions@data)[i],"_alloc", sep="")

  results[q_name] <- NA
  results[a_name] <- NA
  results[t_name] <- NA
  results[names(treatment.predictions@data)[i]] <- NA
  results["trueTreatment"] <- NA
    
  
  results[q_name][1,] <- sum((treatment.predictions@data[,i] - 
                               treatment.predictions@data$trueTreatment) * 
                               treatment.predictions$treatment.status)
  
  results[t_name][1,] <- sum(abs(treatment.predictions@data[,i] - 
                                   treatment.predictions@data$trueTreatment) * 
                               treatment.predictions$treatment.status)
  
  results[a_name][1,] <- results[t_name][1,] - abs(results[q_name][1,])
  
 
    
  results[names(treatment.predictions@data)[i]][1,] <- sum(treatment.predictions@data[,i] * treatment.predictions$treatment.status, na.rm=TRUE) / sum(treatment.predictions$treatment.status)
  results["trueTreatment"][1,] <- sum(treatment.predictions@data$trueTreatment * treatment.predictions$treatment.status) / sum(treatment.predictions$treatment.status)
}



# -----------------------------------------------------------------------------
# Save Summary Parameters
# -----------------------------------------------------------------------------

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
results["trt_prc"] <- NA

results["spill.magnitude"][1,] <- spill.magnitude
results["xvar_error_psill"][1,]  <- xvar_error_psill
results["mod_error_psill"][1,]  <- mod_error_psill
results["trt_spill_sill"][1,]  <- trt_spill_sill
results["xvar_psill"][1,]  <- xvar_psill
results["var1.vrange"][1,] <- var1.vrange
results["prop_acc"][1,] <- prop_acc
results["mod_error.vrange"][1,] <- mod_error.vrange
results["mod_error.magnitude"][1,] <- mod_error.magnitude
results["spill.vrange"][1,] <- spill.vrange
results["beta"][1,] <- beta
results["theta"][1,] <- theta
results["var1_error.vrange"][1,] <- var1_error.vrange
results["caliper"][1,] <- cal
results["sample_size"][1,] <- sample_size
results["tree_split_lim"][1,] <- per_split_lim
results["nrandom"][1,] <- nrandom
results["ct_split_count"][1,] <- length(unique(fit_ctpred$where))
results["trt_prc"][1,] <- trt_prc




# -----------------------------------------------------------------------------
# Output Visualizations
# -----------------------------------------------------------------------------
 map_trt <- treatment.predictions[names(treatment.predictions) != "tot.pred"]
 map_trt <- map_trt[names(map_trt) != "baseline"]
 
 
 #Make the treated areas more evident
 map_trt@data$treatment.status <- map_trt@data$treatment.status * 99
 
 pal = brewer.pal(9,"Greens")
 brks = c(0.0,0.75,1.25,1.5,1.75,2.0,2.5,100)
map_out_path <- paste(map_out, "vsm",version, "_", spill.magnitude, ".png", sep="")
title.v <- paste("Mag:",spill.magnitude," Treat Range:",spill.vrange," RunID:", version, sep="")
png(filename=map_out_path)
print(spplot(map_trt, zcol=names(map_trt)[names(map_trt) != "id"],cuts=brks,col.regions=pal, col="transparent",
       main = list(label=title.v), cex=0.5))

 dev.off()

 
 
# -----------------------------------------------------------------------------
# Record Simulation Results
# -----------------------------------------------------------------------------

write.csv(results,file=out_path)
write.csv(trans_dta@data,file=out_itDta_path)


