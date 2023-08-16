#Code for detecting living fossils on a phylogentic tree using ED (for comparison to EvoHeritage)

# clear workspace
rm(list=ls()) 

# source the functions we're wanting to use
source("../EvoHeritage Tools.R")

# a separate file giving the paths and filenames
source("./Living Fossil Configuration.R")

# record the time when the job started
start.time <- as.numeric(proc.time()[3])

for (j in 1:length(tree.file.paths)) {
  print(paste("simulating tree",j))
  # read in current tree
  current.tree <- ape::read.tree(file=tree.file.paths[j]) 
  result.data <- (Calculate_ED(input.tree = current.tree))
  # save the results
  write.csv(result.data , file = paste(results.file.path,"LivingFossil",str_pad(j, 4, pad = "0"),"ED",sep="") , row.names = TRUE) 
  # print out the elapsed time is helpful 
  current.time <-  as.numeric(proc.time()[3])
  print(paste("time elapsed (hours)",(current.time-start.time)/3600,sep=" "))
}


