#Code for drawing graphs showing results of detecting living fossils on a phylogentic tree using partitioned EvoHeritage and comparing with ED

# clear workspace
rm(list=ls()) 

# install useful packages
require(stringr)

# source the functions we're wanting to test
source("../EvoHeritage Tools.R")
source("./Living Fossil Configuration.R") # keep paths and various other details in a separate config file

# start by defining useful functions for data processing ready to make the plots

# this function prepares a list of file names given the way we know file names were constructed for a given scenario
prepare.filenames <- function(scenario.name,num.files) {
  filenames <- c()
  for (i in 1:num.files) {
    # calculate the next file name
    filename.temp <- paste(results.file.path,"LivingFossil",str_pad(i, 4, pad = "0"),scenario.name,sep="")
    filenames <- c(filenames,filename.temp) # append next file name to final result
  }
  return(filenames)
}

# this function processes files with names given by the input filename.list
# it returns summary data from all files - each file is considered to be an output from a different tree in a distribution of trees
process.LF.files <- function (filename.list) {
  # load the first file just to get the list of tip labels
  result.data <- read.csv(file = filename.list[1], header = TRUE)
  result.data <- result.data[order(result.data$tip.label),]
  # make a dataframe to hold all the data from all files
  partitioned.EvoHeritage.temp <- data.frame(result.data$tip.label) 
  
  for (i in filename.list) {
    # loop over list of files
    result.data <- read.csv(file = i, header = TRUE) # read in data
    result.data <- result.data[order(result.data$tip.label),] #order by tip.label (for later matching)
    partitioned.EvoHeritage.temp <- data.frame(partitioned.EvoHeritage.temp,result.data$partitioned.EvoHeritage) # add data as an extra column
  }
  # now we have a column for each tree (that means for each input file)
  # we need a new data frame to hold summary data - initialise with the correct columns but filled with rubbish initially
  partitioned.EvoHeritage.summary <- data.frame(result.data$tip.label) 
  partitioned.EvoHeritage.summary$LF.min <- rep(-1,length=nrow(partitioned.EvoHeritage.summary))
  partitioned.EvoHeritage.summary$LF.lq <- rep(-1,length=nrow(partitioned.EvoHeritage.summary))
  partitioned.EvoHeritage.summary$LF.median <- rep(-1,length=nrow(partitioned.EvoHeritage.summary))
  partitioned.EvoHeritage.summary$LF.uq <- rep(-1,length=nrow(partitioned.EvoHeritage.summary))
  partitioned.EvoHeritage.summary$LF.max <- rep(-1,length=nrow(partitioned.EvoHeritage.summary))
  for (j in 1:nrow(partitioned.EvoHeritage.summary)) {
    # loop over all tips in the files
    # extract a vector of the values for this top across all files then transform
    transformed.data <- log(as.numeric(partitioned.EvoHeritage.temp[j,2:(length(filename.list)+1)])*1000000)/log(10)
    #get the quantiles
    quantiles.temp <- quantile(transformed.data) 
    # store the quantiles in the final data frame 
    partitioned.EvoHeritage.summary$LF.min[j] <- as.numeric(quantiles.temp[1])
    partitioned.EvoHeritage.summary$LF.lq[j] <- as.numeric(quantiles.temp[2])
    partitioned.EvoHeritage.summary$LF.median[j] <- as.numeric(quantiles.temp[3])
    partitioned.EvoHeritage.summary$LF.uq[j] <- as.numeric(quantiles.temp[4])
    partitioned.EvoHeritage.summary$LF.max[j] <- as.numeric(quantiles.temp[5])
  }
  # add rank data
  partitioned.EvoHeritage.summary$rank <- 1+length(partitioned.EvoHeritage.summary$LF.median)-rank(partitioned.EvoHeritage.summary$LF.median)
 
  # replace tip labels text "_" with a " " to make later graphs easier to make look good
  partitioned.EvoHeritage.summary$result.data.tip.label <- str_replace_all(partitioned.EvoHeritage.summary$result.data.tip.label, "[_]" ," ")
  return(partitioned.EvoHeritage.summary)
}

# this is another way to combine the input files - easier for box plots across small numbers of tips.
# the idea is to flatten all the original data into one data frame
# but only storing data for a subset of tips given by tip.list so it doesn't get out of hand
# which has a range of results for each tip - one result for each tree / input file
flat.process.LF.files <- function (filename.list,tip.list) {
  # create an empty data frame
  partitioned.EvoHeritage.temp <- data.frame()
  for (i in 1:length(filename.list)) {
    # loop over the filename list
    temp.data <- read.csv(file = filename.list[i], header = TRUE) # open current file
    temp.data$tip.label <- str_replace_all(temp.data$tip.label, "[_]" ," ") # replace "_" with " " in tip labels
    subsetted.data <- data.frame() # create another data frame to hold the sub setted data
    for (j in 1:length(tip.list)) {
      # loop over the tip list
      if (j == 1) {
        # we're on the first one so just save the subset
        subsetted.data <- subset(temp.data, temp.data$tip.label == tip.list[j])
      } else {
        # we're on some other one so glue the rows of the main result dataframe to the rows of the newly prepared subset
        subsetted.data <- rbind(subsetted.data, subset(temp.data, temp.data$tip.label == tip.list[j]))
      }
    }
    # re-build the dataframe with new column names
    subsetted.data <- data.frame(subsetted.data$tip.label,subsetted.data$partitioned.EvoHeritage)
    names(subsetted.data) <- c('tip.label', 'partitioned.EvoHeritage')
    if (i == 1) {
      #we're processing the first file just put what we have in the results
      partitioned.EvoHeritage.temp <- subsetted.data
    } else {
      # need to glue the summary across all tips from the current file to the main dataframe of existing results
      partitioned.EvoHeritage.temp <- rbind(partitioned.EvoHeritage.temp, subsetted.data)
    }
  }
  # re order the tip labels and transform the final scores
  partitioned.EvoHeritage.temp$tip.label <- factor(partitioned.EvoHeritage.temp$tip.label , levels=tip.list)
  partitioned.EvoHeritage.temp$partitioned.EvoHeritage <- log(partitioned.EvoHeritage.temp$partitioned.EvoHeritage*1000000)/log(10)
  
  return(partitioned.EvoHeritage.temp)
}

# prepare all date for a given scenario
full.process <- function(scenario.name,num.files=100,top.limit=16) {
  result <- list()
  print(scenario.name)
  files <- prepare.filenames(scenario.name,num.files)
  result[[1]] <- process.LF.files(files) # this is the basic data - used for paired scatter plot
  result[[2]] <- result[[1]][order(-result[[1]]$LF.median),] # reorded data by median used to pick top few
  tip.list <-  result[[2]][1:top.limit,]$result.data.tip.label # list top few species for more data
  tip.list <- rev(tip.list) # so final box plot goes top to bottom 
  result[[3]] <- flat.process.LF.files(files,tip.list) # for box plots to show uncertainty around top few
  write.csv(result[[1]], paste(csv.file.path,"Livingfossilness_results_",scenario.name,".csv",sep=""), row.names=FALSE)
  return(result) # return the full set of pre processed data for this scenario
}

# now use the files for data processing. The scneario names are hard coded and need to match the data generation code

CretaceousAll2.res <- full.process("CretaceousAll2") # rho = 10^-2 using EvoHeritage from Cretaceous and before
JurassicAll2.res <- full.process("JurassicAll2") # rho = 10^-2 using EvoHeritage from Jurassic and before
JurassicAll3.res <- full.process("JurassicAll3") # rho = 10^-3 using EvoHeritage from Jurassic and before
CretaceousAll3.res <- full.process("CretaceousAll3") # rho = 10^-3 using EvoHeritage from Cretaceous and before
ED.res <- full.process("ED") # classic ED
EDcomp.res <- full.process("EDcomp") # rho = 0 using EvoHeritage from Quaternary and before (this is mathematically identical to ED but included as a demonstration of this fact)
EDcomprho.res <- full.process("EDcomprho") # rho = 10^-2 using EvoHeritage from Quaternary and before (this is all EvoHeritage as Quaternary stretches to present day)
CretaceousAllED.res <- full.process("CretaceousAllED") # rho = 0 using EvoHeritage from Cretaceous and before (this is an augmented verison of ED)
JurassicAllED.res <- full.process("JurassicAllED") # rho = 0 using EvoHeritage from Jurassic and before (this is an augmented verison of ED)

# make paired scatter plot

# start by compiling the data for it based on the ranks for all cases
Jurassic0 <- JurassicAllED.res[[1]]$rank
Jurassic0.01 <- JurassicAll2.res[[1]]$rank
Cretaceous0.01 <- CretaceousAll2.res[[1]]$rank
Quaternary0.01 <- EDcomprho.res[[1]]$rank
Cretaceous0 <- CretaceousAllED.res[[1]]$rank
Quaternary0 <- EDcomp.res[[1]]$rank
ED <- ED.res[[1]]$rank
full.data <- data.frame(Jurassic0,Jurassic0.01,Cretaceous0,Cretaceous0.01,Quaternary0,Quaternary0.01,ED)

# now open a file and make the plot
pdf(file = paste(graphs.file.path,"Supplementary Figure Pairs.pdf",sep=""),  width = 8, height = 8)
pairs(~Jurassic0.01 + Cretaceous0.01 + Quaternary0.01 + Jurassic0 + Cretaceous0 + Quaternary0 + ED, data = full.data,cex=0.1,labels=c(expression(paste("Jurassic , ",rho,"=",0.01)),expression(paste("Cretaceous , ",rho,"=",0.01)),expression(paste("Quaternary , ",rho,"=",0.01)),expression(paste("Jurassic , ",rho,"=",0)),expression(paste("Cretaceous , ",rho,"=",0)),expression(paste("Quaternary , ",rho,"=",0)), "ED"))
dev.off()

# make rank and box plots

# open a file for a 4 panel plot
pdf(file = paste(graphs.file.path,"Figure9.pdf",sep=""),  width = 8, height = 8)
par(mfrow = c(2, 2))

# make rank plot for Jurassic case
plot(1:length(JurassicAll3.res[[2]]$LF.median),JurassicAll3.res[[2]]$LF.median,type = "l",log="x",col=rgb(0.65,0.75,1),lty=1,lwd=2, xlab="Species Rank", ylab="Median Living-fossil-ness", ylim = c(3.5,7.2))
lines(1:length(JurassicAll2.res[[2]]$LF.median),JurassicAll2.res[[2]]$LF.median,type = "l",col="black",lty=1)
legend(200, 7, legend=c(expression(paste(rho,"=",10^-2)), expression(paste(rho,"=",10^-3))), col=c("black", rgb(0.65,0.75,1)), lty=1, lwd=c(1,2), cex=1,
       box.lty=0)
rect(1, 0, 16, 8, col = rgb(0, 0, 0, max = 255, alpha = 20),border=FALSE)

# make boxplot for Jurassic case
text(4,4.3,"Top 16", cex = 1)
text(20000000,8,"Jurassic Mammalian Living Fossils", cex = 1.5, pos = 2,xpd=NA)
op <- par(mar = c(5,11,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
boxplot(partitioned.EvoHeritage~tip.label,data=JurassicAll2.res[[3]],  horizontal=TRUE, 
        xlab=expression(paste("Living-fossil-ness ",rho,"=",10^-2)), ylab="",las=1,cex=0.5,yaxt="n")
axis(2,font.axis=3,labels = JurassicAll2.res[[3]]$tip.label[1:16], at = 1:16,las=2)
op <- par(mar = c(5,4,4,2) + 0.1) # set back to default values

# make the rank plot for Cretaceous case
plot(1:length(CretaceousAll2.res[[2]]$LF.median),CretaceousAll3.res[[2]]$LF.median,type = "l",log="x",col=rgb(0.65,0.75,1),lty=1, lwd=2,xlab="Species Rank", ylab="Median Living-fossil-ness", ylim = c(4,7.5))
lines(1:length(CretaceousAll3.res[[2]]$LF.median),CretaceousAll2.res[[2]]$LF.median,type = "l",col="black", lty=1)
legend(200, 7, legend=c(expression(paste(rho,"=",10^-2)), expression(paste(rho,"=",10^-3))), col=c("black", rgb(0.65,0.75,1)), lty=c(1,1), lwd=c(1,2), cex=1,
       box.lty=0)
rect(1, 0, 16, 8, col = rgb(0, 0, 0, max = 255, alpha = 20),border=FALSE)

# make boxplot for Cretaceous case
text(4,4.3,"Top 16", cex = 1)
text(20000000,8,"Cretaceous Mammalian Living Fossils", cex = 1.5, pos = 2,xpd=NA)
op <- par(mar = c(5,11,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
boxplot(partitioned.EvoHeritage~tip.label,data=CretaceousAll2.res[[3]],  horizontal=TRUE, 
        xlab=expression(paste("Living-fossil-ness ",rho,"=",10^-2)), ylab="",las=1,cex=0.5,yaxt="n")
axis(2,font.axis=3,labels = CretaceousAll2.res[[3]]$tip.label[1:16], at = 1:16,las=2)

# close file
dev.off()

# supplementary data giving ranks for ED

# open the file and prepare for a 2 panel plot
pdf(file = paste(graphs.file.path,"Supplementary Figure ED.pdf",sep=""),  width = 10, height = 5)
par(mfrow = c(1, 2))

# make the rank plot
plot(1:length(ED.res[[2]]$LF.median),ED.res[[2]]$LF.median,type = "l",log="x",col="black",xlab="Species Rank", ylab="Median Log(ED)", ylim = c(5.75,8))
rect(1, 0, 16, 8.2, col = rgb(0, 0, 0, max = 255, alpha = 20),border=FALSE)

# make the box plot
text(4,4.3,"Top 16", cex = 1)
text(20000000,8.3,"ED-based Mammalian Living Fossils", cex = 1.5, pos = 2,xpd=NA)
op <- par(mar = c(5,11,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
boxplot(partitioned.EvoHeritage~tip.label,data=ED.res[[3]],  horizontal=TRUE, 
        xlab=expression(paste("Log(ED)")), ylab="",las=1,cex=0.5,yaxt="n")
axis(2,font.axis=3,labels = ED.res[[3]]$tip.label[1:16], at = 1:16,las=2)

# close the file
dev.off()


