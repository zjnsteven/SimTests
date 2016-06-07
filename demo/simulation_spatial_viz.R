#Folder that contains CSVs to be visualized
rm()
library(Cairo)
res_folder <- "/mnt/sc/H"
files <- list.files(path=res_folder, full.names=T, recursive=FALSE, pattern="\\.csv$")

files <- files[-grep("dta",files)]

print(length(files))

for(path in 1:length(files))
{
  if(path == 1)
  {
    full_df <- read.csv(files[path])
  }
  else
  {
    full_df <-  rbind(full_df,read.csv(files[path]))
  }
}

results <- full_df

#Visualization of Simulation results
library(plotly)
viz.sims <- function(results, varH, mtitle, pre="")
{
  results.plot <- results
  alpha.set = 0.1
  cex.set = 0.1
  
  
  eval(parse(text=paste("results.plot$v1 <- results.plot$",varH,sep="")))
  results.plot <- results.plot[order(results.plot$v1),]
  ylower <- 1.0#min(results.plot[paste(pre,"ct.spill",sep="")]) - (2*abs(min(results.plot[paste(pre,"trueTreatment",sep="")])))
  yupper <- 1.5#max(results.plot[paste(pre,"trueTreatment",sep="")]) 
  #if(mtitle == "ATE by Model")
  #{
  #  yupper <- max(results.plot[paste(pre,"propensity.tree",sep="")]) * 1.25
  #}
  plot(ylim=c(ylower,yupper), 
       results.plot$v1, 
       results.plot[paste(pre,"baseline.lm",sep="")][[1]], 
       col=rgb(1,0,0,alpha=alpha.set), pch=3, cex=cex.set,
       main=mtitle,
       ylab="Estimate",
       xlab=varH)
  lines(lowess(results.plot$v1,results.plot[paste(pre,"baseline.lm",sep="")][[1]]), col=rgb(1,0,0), pch=3)
  
  lines(lowess(results.plot$v1, 
               results.plot[paste(pre,"baseline.matchit",sep="")][[1]]), col=rgb(0,0,1), pch=4)
  points(results.plot$v1, 
         results.plot[paste(pre,"baseline.matchit",sep="")][[1]], col=rgb(0,0,1,alpha=alpha.set), pch=4, cex=cex.set)
  
  lines(lowess(results.plot$v1, 
               results.plot[paste(pre,"trueTreatment",sep="")][[1]]), col=rgb(0,1,0), pch=4)
  points(results.plot$v1, 
         results.plot[paste(pre,"trueTreatment",sep="")][[1]], col=rgb(0,1,0,alpha=alpha.set), pch=4, cex=cex.set)
  
  lines(lowess(results.plot$v1, 
               results.plot[paste(pre,"ct.pred",sep="")][[1]]), col=rgb(1,0.5,0), pch=2)
  points(results.plot$v1, 
         results.plot[paste(pre,"ct.pred",sep="")][[1]], col=rgb(1,0.5,0,alpha=alpha.set), pch=2, cex=cex.set)
  
  lines(lowess(results.plot$v1, 
               results.plot[paste(pre,"gwr",sep="")][[1]]), col=rgb(0,1,1), pch=4)
  points(results.plot$v1, 
                results.plot[paste(pre,"gwr",sep="")][[1]], col=rgb(0,1,1,alpha=alpha.set), pch=4, cex=cex.set)
  
  lines(lowess(results.plot$v1, 
               results.plot[paste(pre,"spatialPSM",sep="")][[1]]), col=rgb(0.5,0.5,0.5), pch=4)
  points(results.plot$v1, 
         results.plot[paste(pre,"spatialPSM",sep="")][[1]], col=rgb(0.5,0.5,0.5,alpha=alpha.set), pch=4, cex=cex.set)
  
  #lines(lowess(results.plot$v1, 
  #               results.plot[paste(pre,"tot.pred",sep="")][[1]]), col=rgb(0,0,0), pch=3)
  #points(results.plot$v1, 
  #         results.plot[paste(pre,"tot.pred",sep="")][[1]], col=rgb(0,0,0,alpha=alpha.set), pch=3, cex=cex.set)

  
  legend("topleft",
         cex = 0.65,
         legend=c("Baseline LM","True ATE", "Baseline MatchIt", 
                  "CT", "NA", 
                  "GWR", "Spatial Thresh"), 
         pch=c(pch = 3, pch=1, pch=4, pch=2, pch=3, pch=4, pch=2),
         col=c(col="red", col="green", col="blue", col="orange", 
               col="black", col=109, col=144), title = "Legend")
}

#Viz Creation

viz.list=c("nrandom","beta", "tree_split_lim", "caliper", "spill.magnitude", "spill.vrange", "var1.vrange", "trt_spill_sill", "tree_thresh", "trt_prc", "thresh_est")

viz.list=c("trt_prc", "beta", "spill.magnitude", "caliper", "tree_thresh", "thresh_est", "trtcon_overlap")

for(type in viz.list)
{
fname = paste("/home/aiddata/Desktop/SimViz/",type,"_",length(files),".png",sep="")
CairoPNG(1600,900,file=fname, bg="white")
title <- paste("ATE by Model", length(files), sep="")
print(type)
print("__")
viz.sims(results, type, title)
dev.off()
}
#Break out of Cairo dev
dev.off()




quantity_alloc <- function(results)
{
hist(results$ct_split_count)
allocs <- cbind(results$baseline.lm_alloc, results$spatialPSM_alloc, results$ct.pred_alloc, results$baseline.matchit_alloc, results$gwr_alloc)
quants <- cbind(results$baseline.lm_quant, results$spatialPSM_quant, results$ct.pred_quant, results$baseline.matchit_quant, results$gwr_quant)
tots <- cbind(results$baseline.lm_tot, results$spatialPSM_tot, results$ct.pred_tot, results$baseline.matchit_tot, results$gwr_tot)

allocs <- allocs / (results$trt_prc * results$nrandom)
quants <- quants / (results$trt_prc * results$nrandom)
tots <- tots / (results$trt_prc * results$nrandom)

par(xpd=TRUE)

boxplot(allocs, use.cols=TRUE, col=c("red", "blue", "green", "yellow", "orange"), xaxt="n", main="Allocation Error")
legend("bottom",c("OLS", "Spatial PSM", "CT", "PSM", "GWR"), fill=c("red", "blue", "green", "yellow", "orange"), inset=-0.1, horiz=TRUE, cex=0.5)

boxplot(quants, use.cols=TRUE, col=c("red", "blue", "green", "yellow", "orange"), xaxt="n", main="Quantity Error")
legend("bottom",c("OLS", "Spatial PSM", "CT","PSM", "GWR"), fill=c("red", "blue", "green", "yellow", "orange"), inset=-0.1, horiz=TRUE, cex=0.5)
#dev.off()
#boxplot(tots, use.cols=TRUE, col=c("red", "blue", "green", "yellow"), xaxt="n", main="Total Error")
#legend("bottom",c("OLS", "Spatial PSM", "CT","PSM", "GWR"), fill=c("red", "blue", "green", "yellow"), inset=-0.1, horiz=TRUE, cex=0.5)
}

quantity_alloc(results)


