
# -----------------------------------------------------------------------------
#Generating Dataframe
# -----------------------------------------------------------------------------

# generate random coordinates
random.longitude <- runif(nrandom, minx, maxx)
random.latitude <- runif(nrandom, miny, maxy)

# create dataframe
spdf <- data.frame(id = 1:nrandom,
                   longitude = random.longitude,
                   latitude = random.latitude)

# convert to spatial points dataframe
coordinates(spdf) <- c("longitude", "latitude")
proj4string(spdf) <- CRS("+proj=longlat +datum=WGS84")
prj <- proj4string(spdf)
spdf <- spTransform(spdf, CRS(prj))

# -----------------------------------------------------------------------------
#Generating Ancillary Data
# -----------------------------------------------------------------------------


# using gstat to generate fields with spatial autocorrelation
# source:
#   http://santiago.begueria.es/2010/10/
#     generating-spatially-correlated-random-fields-with-r/

# Define the spatial correlation of the x-variable
var.g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                     model=vgm(psill=xvar_psill, model="Sph", var1.vrange),
                     nmax=100)

# make simulations based on the gstat object
var.sim <- predict(var.g.dummy, newdata=spdf, nsim=1)

var.sim$sim1 <- var.sim$sim1 / max(var.sim$sim1)

spdf@data$trueVar <- var.sim$sim1

# -----------------------------------------------------------------------------
#Generating Ancillary Error Data
# -----------------------------------------------------------------------------

# Define the spatial correlation of the x-variable
var.error.g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                           model=vgm(psill=xvar_error_psill, model="Sph", var1_error.vrange),nmax=100)

# make simulations based on the gstat object
var.error.sim <- predict(var.error.g.dummy, newdata=spdf, nsim=1)

var.error.sim$sim1 <- var.error.sim$sim1  / max(var.error.sim$sim1)

spdf@data$modelVarError <- var.error.sim$sim1* (1-prop_acc)
spdf@data$modelVar <- (var.error.sim$sim1 * (1-prop_acc)) + var.sim$sim1


# -----------------------------------------------------------------------------
#Generating General Model Error Data
# -----------------------------------------------------------------------------

# Define the spatial correlation of the x-variable
mod.error.g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                           model=vgm(psill=mod_error_psill, model="Sph", mod_error.vrange), nmax=100)

# make simulations based on the gstat object
mod.error.sim <- predict(mod.error.g.dummy, newdata=spdf, nsim=1)
if(max(mod.error.sim$sim1) > 0)
{
mod.error.sim$sim1 <- mod.error.sim$sim1 / max(mod.error.sim$sim1)
}

spdf@data$modelError <- mod.error.sim$sim1 * mod_error.magnitude

# -----------------------------------------------------------------------------
#Generating the Treatment Variable
# -----------------------------------------------------------------------------

temp_rand <- spdf@data$trueVar + (runif(length(spdf@data$trueVar),0.0,trtcon_overlap) * spdf@data$trueVar)

treatment.binary <- ifelse(temp_rand> quantile(temp_rand, 
                                                        (1-trt_prc)), 1, 0)

spdf$treatment.status = treatment.binary




# -----------------------------------------------------------------------------
#Outcome Model
# -----------------------------------------------------------------------------

# Define the spatial spillover variogram
#trt.spillover <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
#                       model=vgm(psill=trt_spill_sill, model="Sph", spill.vrange), nmax=100)
#trt.spillover.sim <- predict(trt.spillover, newdata=spdf, nsim=1)
#vgm.spillover <- variogram(sim1~1, trt.spillover.sim, covariogram=TRUE)


#vgm.polynomial <- lm(gamma ~ poly(dist, degree=3, raw=TRUE),
#                     data=vgm.spillover)

#vgm.spillover$fit <- predict(vgm.polynomial, vgm.spillover)

#print(summary(vgm.polynomial))
#print(summary(vgm.spillover$fit))


#Calculate the distance to treated units for each unit (excluding itself)
spdf@data$tmp.spillover.weights <- NA
spdf@data$trans_dist <- NA
for (i in 1:nrandom) {
  
  dist <- spDists(spdf[i, ]@coords, spdf@coords, longlat=TRUE)[1,]
  
  
  neighbors <- dist > 0

  #Limit to the min (after the final predicted distance all values are 0)
  gamma <- trt_spill_sill - (trt_spill_sill * (((3/2)*((dist)/spill.vrange)) - 
                                                 ((1/2) * (dist)/spill.vrange)^3))
  #spdf@data$gamma <- predict(vgm.polynomial, newdata=t)
  gamma[gamma < 0] <- 0
  neighbors[dist > spill.vrange] <- 0
  #tmp.newdata$gamma <- 1/tmp.newdata$gamma
  #tmp.newdata$gamma[is.infinite(tmp.newdata$gamma)] <- 0#max(tmp.newdata$gamma[!is.infinite(tmp.newdata$gamma)])
  spdf@data$tmp.spillover.weights[i] <- sum(gamma * neighbors * spdf@data$treatment.status)
  noZdist <- (dist * neighbors)
  noZdist <- noZdist[noZdist != 0]
  spdf@data$trans_dist[i] <- 1 / (((spdf@data$treatment.status[i]) * min(noZdist))^3 - ((1-spdf@data$treatment.status[i]) * min(noZdist))^3)
  spdf@data$dist_avg[i] <- mean(dist)
  spdf@data$neighbors <- sum(neighbors)
  }
print(head(spdf@data))

spdf@data$tmp.spillover.weights[is.na(spdf@data$tmp.spillover.weights)] <- 0
#vgm.spillover$gamma[vgm.spillover$gamma < 0] <- 0
print("tmp.spillover")
print(summary(spdf@data$tmp.spillover.weights))



#Total treated effect across the study area
tmp.spillover.t1 <- sum(theta * spdf$treatment.status)

#Ratio of spillover to assign to each unit
tmp.spillover.t2 <- c()
sum.spillover.weights <- sum(spdf@data$tmp.spillover.weights, na.rm=TRUE)
for (i in 1:nrandom) {
  spdf@data$tmp.spillover.t2[i] <- (spdf@data$tmp.spillover.weights[i] / sum.spillover.weights)
  if(is.na(spdf@data$tmp.spillover.t2[i])){
    tmp.spillover.t2[i] = 0
  }
}

#Calculate total spillover
tot_spill <- tmp.spillover.t1 * spill.magnitude
print("Tot spill")
print(tot_spill)
print(summary(spdf@data$tmp.spillover.t2))
#Calculate per-unit spillover
spdf$trueSpill <-spdf@data$tmp.spillover.t2 * tot_spill

print(summary(spdf$trueSpill))

#Calculate per-unit true outcome
spdf$trueOutcome <- spdf$trueSpill + (beta * spdf$trueVar) + (spdf$treatment.status * theta)

print(summary(spdf$trueOutcome))
print(summary(spdf$treatment.status))

spdf$modelOutcome <- spdf$trueOutcome + spdf$modelError

back_df.2 <- spdf@data[,names(spdf@data) %in% c("tmp.spillover.weights", "tmp.spillover.t2", "dist_avg", "neighbors")]
spdf@data <- spdf@data[,!names(spdf@data) %in% c("tmp.spillover.weights", "tmp.spillover.t2", "dist_avg", "neighbors")]






