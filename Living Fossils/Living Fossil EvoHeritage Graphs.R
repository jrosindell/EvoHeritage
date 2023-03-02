#Code for drawing graphs showing results of detecting living fossils on a phylogentic tree using partitioned EvoHeritage

# clear workspace
rm(list=ls()) 

# install useful packages
require(stringr)

# source the functions we're wanting to test
source("../EvoHeritage Tools.R")
source("./Living Fossil Configuration.R")

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
  # re order by median value
  partitioned.EvoHeritage.summary <- partitioned.EvoHeritage.summary[order(-partitioned.EvoHeritage.summary$LF.median),]
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

full.process <- function(scenario.name,num.files=100,top.limit=16) {
  result <- list()
  print(scenario.name)
  files <- prepare.filenames(scenario.name,num.files)
  result[[1]] <- process.LF.files(files)
  tip.list <-  result[[1]][1:top.limit,]$result.data.tip.label
  tip.list <- rev(tip.list) # so final box plot goes top to bottom 
  result[[2]] <- flat.process.LF.files(files,tip.list)
  return(result)
}

JurassicNew1.res <- full.process("JurassicNew1")
CretaceousNew1.res <- full.process("CretaceousNew1")
JurassicAll1.res <- full.process("JurassicAll1")
CretaceousAll1.res <- full.process("CretaceousAll1")

JurassicNew2.res <- full.process("JurassicNew2")
CretaceousNew2.res <- full.process("CretaceousNew2")
JurassicAll2.res <- full.process("JurassicAll2")
CretaceousAll2.res <- full.process("CretaceousAll2")

JurassicNew3.res <- full.process("JurassicNew3")
CretaceousNew3.res <- full.process("CretaceousNew3",49)
JurassicAll3.res <- full.process("JurassicAll3")
CretaceousAll3.res <- full.process("CretaceousAll3",39)

PaleogeneAll2.res <- full.process("PaleogeneAll2",11)

JurassicAll2.res.more <- full.process("JurassicAll2",top.limit=30)
CretaceousAll2.res.more <- full.process("CretaceousAll2",top.limit=30)
JurassicAll3.res.more <- full.process("JurassicAll3",top.limit=30)
CretaceousAll3.res.more <- full.process("CretaceousAll3",39,top.limit=30)

JurassicNew2.res.more <- full.process("JurassicNew2",top.limit=30)
CretaceousNew2.res.more <- full.process("CretaceousNew2",top.limit=30)
JurassicNew3.res.more <- full.process("JurassicNew3",top.limit=30)
CretaceousNew3.res.more <- full.process("CretaceousNew3",39,top.limit=30)

pdf(file = paste(graphs.file.path,"Main.pdf",sep=""),  width = 8, height = 8)


par(mfrow = c(2, 2))

plot(1:length(JurassicNew3.res[[1]]$LF.median),JurassicNew3.res[[1]]$LF.median,type = "l",log="x",col="red",xlab="Species Rank", ylab="Median Living-fossil-ness", ylim = c(4,7.5))
lines(1:length(JurassicNew2.res[[1]]$LF.median),JurassicNew2.res[[1]]$LF.median,type = "l",col="blue")
legend(200, 7, legend=c(expression(paste(rho,"=",10^-2)), expression(paste(rho,"=",10^-3))), col=c("blue", "red"), lty=1, cex=1,
       box.lty=0)


rect(1, 0, 16, 8, col = rgb(0, 0, 0, max = 255, alpha = 25),border=FALSE)

text(4,4.3,"Top 16", cex = 1)

text(20000000,8,"Jurassic Mammalian Living Fossils", cex = 1.5, pos = 2,xpd=NA)

op <- par(mar = c(5,11,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
boxplot(partitioned.EvoHeritage~tip.label,data=JurassicNew2.res[[2]],  horizontal=TRUE, 
        xlab=expression(paste("Living-fossil-ness ",rho,"=",10^-2)), ylab="",las=1,cex=0.5,yaxt="n")
axis(2,font.axis=3,labels = JurassicNew2.res[[2]]$tip.label[1:16], at = 1:16,las=2)

op <- par(mar = c(5,4,4,2) + 0.1) # set back to default values


plot(1:length(CretaceousNew2.res[[1]]$LF.median),CretaceousNew2.res[[1]]$LF.median,type = "l",log="x",col="blue",xlab="Species Rank", ylab="Median Living-fossil-ness", ylim = c(4,7.5))
lines(1:length(CretaceousNew3.res[[1]]$LF.median),CretaceousNew3.res[[1]]$LF.median,type = "l",col="red")
legend(200, 7, legend=c(expression(paste(rho,"=",10^-2)), expression(paste(rho,"=",10^-3))), col=c("blue", "red"), lty=1, cex=1,
       box.lty=0)

rect(1, 0, 16, 8, col = rgb(0, 0, 0, max = 255, alpha = 25),border=FALSE)

text(4,4.3,"Top 16", cex = 1)

text(20000000,8,"Cretaceous Mammalian Living Fossils", cex = 1.5, pos = 2,xpd=NA)

op <- par(mar = c(5,11,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
boxplot(partitioned.EvoHeritage~tip.label,data=CretaceousNew2.res[[2]],  horizontal=TRUE, 
        xlab=expression(paste("Living-fossil-ness ",rho,"=",10^-2)), ylab="",las=1,cex=0.5,yaxt="n")
axis(2,font.axis=3,labels = CretaceousNew2.res[[2]]$tip.label[1:16], at = 1:16,las=2)



dev.off()




pdf(file = paste(graphs.file.path,"SI_1.pdf",sep=""),  width = 8, height = 8)


par(mfrow = c(2, 2))

plot(1:length(JurassicAll3.res[[1]]$LF.median),JurassicAll3.res[[1]]$LF.median,type = "l",log="x",col="red",xlab="Species Rank", ylab="Median Living-fossil-ness", ylim = c(4,7.5))
lines(1:length(JurassicAll2.res[[1]]$LF.median),JurassicAll2.res[[1]]$LF.median,type = "l",col="blue")
legend(200, 7, legend=c(expression(paste(rho,"=",10^-2)), expression(paste(rho,"=",10^-3))), col=c("blue", "red"), lty=1, cex=1,
       box.lty=0)


rect(1, 0, 16, 8, col = rgb(0, 0, 0, max = 255, alpha = 25),border=FALSE)

text(4,4.3,"Top 16", cex = 1)

text(20000000,8,"Pre-Cretaceous Mammalian Living Fossils", cex = 1.5, pos = 2,xpd=NA)

op <- par(mar = c(5,11,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
boxplot(partitioned.EvoHeritage~tip.label,data=JurassicAll2.res[[2]],  horizontal=TRUE, 
        xlab=expression(paste("Living-fossil-ness",rho,"=",10^-2)), ylab="",las=1,cex=0.5,yaxt="n")
axis(2,font.axis=3,labels = JurassicAll2.res[[2]]$tip.label[1:16], at = 1:16,las=2)

op <- par(mar = c(5,4,4,2) + 0.1) # set back to default values


plot(1:length(CretaceousAll2.res[[1]]$LF.median),CretaceousAll2.res[[1]]$LF.median,type = "l",log="x",col="blue",xlab="Species Rank", ylab="Median Living-fossil-ness", ylim = c(4,7.5))
lines(1:length(CretaceousAll3.res[[1]]$LF.median),CretaceousAll3.res[[1]]$LF.median,type = "l",col="red")
legend(200, 7, legend=c(expression(paste(rho,"=",10^-2)), expression(paste(rho,"=",10^-3))), col=c("blue", "red"), lty=1, cex=1,
       box.lty=0)

rect(1, 0, 16, 8, col = rgb(0, 0, 0, max = 255, alpha = 25),border=FALSE)

text(4,4.3,"Top 16", cex = 1)

text(20000000,8,"Pre-Paleogene Mammalian Living Fossils", cex = 1.5, pos = 2,xpd=NA)

op <- par(mar = c(5,11,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
boxplot(partitioned.EvoHeritage~tip.label,data=CretaceousAll2.res[[2]],  horizontal=TRUE, 
        xlab=expression(paste("Living-fossil-ness ",rho,"=",10^-2)), ylab="",las=1,cex=0.5,yaxt="n")
axis(2,font.axis=3,labels = CretaceousAll2.res[[2]]$tip.label[1:16], at = 1:16,las=2)



dev.off()




pdf(file = paste(graphs.file.path,"SI_2.pdf",sep=""),  width = 12, height = 8)


par(mfrow = c(1, 2))

op <- par(mar = c(5,11.5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

boxplot(partitioned.EvoHeritage~tip.label,data=JurassicNew2.res.more[[2]],  horizontal=TRUE, 
        xlab=expression(paste("Living-fossil-ness")), ylab="",las=1,cex=0.5,yaxt="n")
axis(2,font.axis=3,labels = JurassicNew2.res.more[[2]]$tip.label[1:30], at = 1:30,las=2)

text(7.1,34,expression(paste("Jurassic Mammalian Living Fossils ",rho,"=",10^-2)), cex = 1.5, pos = 2,xpd=NA)


boxplot(partitioned.EvoHeritage~tip.label,data=CretaceousNew2.res.more[[2]],  horizontal=TRUE, 
        xlab=expression(paste("Living-fossil-ness")), ylab="",las=1,cex=0.5,yaxt="n")
axis(2,font.axis=3,labels = CretaceousNew2.res.more[[2]]$tip.label[1:30], at = 1:30,las=2)

text(7.6,34,expression(paste("Cretaceous Mammalian Living Fossils ",rho,"=",10^-2)), cex = 1.5, pos = 2,xpd=NA)

dev.off()


pdf(file = paste(graphs.file.path,"SI_3.pdf",sep=""),  width = 12, height = 8)


par(mfrow = c(1, 2))

op <- par(mar = c(5,11.5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

boxplot(partitioned.EvoHeritage~tip.label,data=JurassicNew3.res.more[[2]],  horizontal=TRUE, 
        xlab=expression(paste("Living-fossil-ness")), ylab="",las=1,cex=0.5,yaxt="n")
axis(2,font.axis=3,labels = JurassicNew3.res.more[[2]]$tip.label[1:30], at = 1:30,las=2)

text(7.5,34,expression(paste("Jurassic Mammalian Living Fossils ",rho,"=",10^-3)), cex = 1.5, pos = 2,xpd=NA)


boxplot(partitioned.EvoHeritage~tip.label,data=CretaceousNew3.res.more[[2]],  horizontal=TRUE, 
        xlab=expression(paste("Living-fossil-ness")), ylab="",las=1,cex=0.5,yaxt="n")
axis(2,font.axis=3,labels = CretaceousNew3.res.more[[2]]$tip.label[1:30], at = 1:30,las=2)

text(7.8,34,expression(paste("Cretaceous Mammalian Living Fossils ",rho,"=",10^-3)), cex = 1.5, pos = 2,xpd=NA)

dev.off()



