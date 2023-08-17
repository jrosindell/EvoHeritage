#Code for drawing graphs showing results of primary productivity example using EvoHeritage, PD and species richness

# clear workspace
rm(list=ls()) 

# source the functions we're wanting to test
source("../EvoHeritage Tools.R")
source("./Primary Productivity Configuration.R") # keep paths and various other details in a separate config file

# open PDF to make the graph in 
pdf(file = paste(graphs.file.path,"Figure8.pdf"),  width = 8, height = 8)
# prepare for 4 panels
par(mfrow = c(2, 2))

# make first plot and draw first line (all species case)
load(file=paste(results.file.path,"EH.PProd.results.rda",sep=""))
plot(EvoH.results$rho,EvoH.results$EvoH_est,log="x",type="l",ylim=c(0.44,0.67),xlim=c(10^-4,1),xlab=expression(paste("Information erosion (",rho,")")),ylab=expression(paste("Correlation ",varphi[rho]," with biomass")),main="Unfiltered data")
points(1,EvoH.results$sr_est[2],pch = 1,cex=1.3)
points(10^-4,EvoH.results$pd_est[2],pch = 16,cex=1.3)
# draw a vertical line to show where the max is
abline(v=10^-4,lty=1,lwd = 0.5)

# save the data to a CSV for sharing 
write.csv(EvoH.results, paste(csv.file.path,"Fig.8.Unfiltered.all.csv"), row.names=FALSE)

# add a new data series to the graph for the case excluding grasses (only Pentapetalae)
load(file=paste(results.file.path,"EH.PProd.results.nograss.rda",sep=""))
lines(EvoH.results$rho,EvoH.results$EvoH_est,type="l",col="red",lty=3,lwd=2)
points(1,EvoH.results$sr_est[2],pch = 1,cex=1.3,col="red")
points(10^-4,EvoH.results$pd_est[2],pch = 16,cex=1.3,col="red")

# get the rho that gives best correlation
max.corr <- EvoH.results$EvoH_est[1]
max.rho <- EvoH.results$rho[1]
for(i in 2:length(EvoH.results$rho)) {
  if ( EvoH.results$EvoH_est[i] > max.corr) {
    max.corr <- EvoH.results$EvoH_est[i]
    max.rho <- EvoH.results$rho[i]
  }
}
abline(v=max.rho,lty=3,col="red",lwd = 0.5)

vopt.unfiltered <- max.rho

max.corr.unfiltered <- max.corr

write.csv(EvoH.results, paste(csv.file.path,"Fig.8.Unfiltered.Pentapetalae.csv"), row.names=FALSE)


# add a new data series to the graph for the case excluding Pentapetalae (only grasses)
load(file=paste(results.file.path,"EH.PProd.results.onlygrass.rda",sep=""))
lines(EvoH.results$rho,EvoH.results$EvoH_est,type="l",col="blue",lty=2,lwd=2)
points(1,EvoH.results$sr_est[2],pch = 1,cex=1.3,col="blue")
points(10^-4,EvoH.results$pd_est[2],pch = 16,cex=1.3,col="blue")

max.corr <- EvoH.results$EvoH_est[1]
max.rho <- EvoH.results$rho[1]
for(i in 2:length(EvoH.results$rho)) {
  if ( EvoH.results$EvoH_est[i] > max.corr) {
    max.corr <- EvoH.results$EvoH_est[i]
    max.rho <- EvoH.results$rho[i]
  }
}
abline(v=max.rho,lty=2,col="blue",lwd = 0.5)

write.csv(EvoH.results, paste(csv.file.path,"Fig.8.Unfiltered.Poaceae.csv"), row.names=FALSE)


# make the second plot
load(file=paste(results.file.path,"EH.PProd.results.drop0.rda",sep=""))
plot(EvoH.results$rho,EvoH.results$EvoH_est,log="x",type="l",ylim=c(0.39,0.47),xlim=c(10^-4,1),xlab=expression(paste("Information erosion (",rho,")")),ylab=expression(paste("Correlation ",varphi[rho]," with biomass")),main="Filtered data")
points(1,EvoH.results$sr_est[2],pch = 1,cex=1.3)
points(10^-4,EvoH.results$pd_est[2],pch = 16,cex=1.3)

legend(10^-4, 0.47, legend=c("all species", expression(paste("Poaceae")),expression(paste("Pentapetalae"))),
       col=c("black","blue", "red"), lty=c(1,2,3), cex=0.8, box.lty=0)

max.corr <- EvoH.results$EvoH_est[1]
max.rho <- EvoH.results$rho[1]
for(i in 2:length(EvoH.results$rho)) {
  if ( EvoH.results$EvoH_est[i] > max.corr) {
    max.corr <- EvoH.results$EvoH_est[i]
    max.rho <- EvoH.results$rho[i]
  }
}
abline(v=max.rho,lty=1,lwd = 0.5)


load(file=paste(results.file.path,"EH.PProd.results.nograss.drop0.rda",sep=""))
lines(EvoH.results$rho,EvoH.results$EvoH_est,type="l",col="red",lty=3,lwd=2)
points(1,EvoH.results$sr_est[2],pch = 1,cex=1.3,col="red")
points(10^-4,EvoH.results$pd_est[2],pch = 16,cex=1.3,col="red")

max.corr <- EvoH.results$EvoH_est[1]
max.rho <- EvoH.results$rho[1]
for(i in 2:length(EvoH.results$rho)) {
  if ( EvoH.results$EvoH_est[i] > max.corr) {
    max.corr <- EvoH.results$EvoH_est[i]
    max.rho <- EvoH.results$rho[i]
  }
}
abline(v=max.rho,lty=3,col="red",lwd = 0.5)

vopt.filtered <- max.rho

max.corr.filtered <- max.corr

write.csv(EvoH.results, paste(csv.file.path,"Fig.8.Filtered.all.csv"), row.names=FALSE)

write.csv(EvoH.results, paste(csv.file.path,"Fig.8.Filtered.Pentapetalae.csv"), row.names=FALSE)

load(file=paste(results.file.path,"EH.PProd.results.onlygrass.drop0.rda",sep=""))
lines(EvoH.results$rho,EvoH.results$EvoH_est,type="l",col="blue",lty=2,lwd=2)
points(1,EvoH.results$sr_est[2],pch = 1,cex=1.3,col="blue")
points(10^-4,EvoH.results$pd_est[2],pch = 16,cex=1.3,col="blue")

max.corr <- EvoH.results$EvoH_est[1]
max.rho <- EvoH.results$rho[1]
for(i in 2:length(EvoH.results$rho)) {
  if ( EvoH.results$EvoH_est[i] > max.corr) {
    max.corr <- EvoH.results$EvoH_est[i]
    max.rho <- EvoH.results$rho[i]
  }
}
abline(v=max.rho,lty=2,col="blue",lwd=0.5)


write.csv(EvoH.results, paste(csv.file.path,"Fig.8.Filtered.Poaceae.csv"), row.names=FALSE)

par(mar = c(9,4, 0, 2))
load(file=paste(results.file.path,"EH.PProd.results.nograss.rda",sep=""))
plot(best.clean.data$EvoH,best.clean.data$biomass,pch=4,cex=0.75,xlab=expression(paste(varphi[rho],"(","Pentapetalae",") for optimal ",rho ,"=0.0016")),ylab="Total biomass")
load(file=paste(results.file.path,"EH.PProd.results.nograss.drop0.rda",sep=""))
plot(best.clean.data$EvoH,best.clean.data$biomass,pch=4,cex=0.75,xlab=expression(paste(varphi[rho],"(","Pentapetalae",") for optimal ",rho ,"=0.16")),ylab="Total biomass")

dev.off()


