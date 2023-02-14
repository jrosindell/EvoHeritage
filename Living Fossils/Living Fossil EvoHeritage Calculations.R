#Code for detecting living fossils on a phylogentic tree using partitioned EvoHeritage

# clear workspace
rm(list=ls()) 

# source the functions we're wanting to test
source("../EvoHeritage Tools.R")
source("./Living Fossil Configuration.R")

# record the time when the job started
start.time <- as.numeric(proc.time()[3])

for(i in 1: nrow(scenarios)) {
  print("simulating new scenario")
  print(scenarios[i,])
  for (j in 1:length(tree.file.paths)) {
    print(paste("simulating tree",j))
    # read in current tree
    current.tree <- ape::read.tree(file=tree.file.paths[j]) 
    # set seed for repeat ability
    set.seed(scenarios$seed[i]) 
    # make an ancestral EvoHeritage Tree
    current.EH.tree <- make.ancestral.evoheritage.tree(current.tree,rho=scenarios$rho[i],lambda = scenarios$lambda[i],min.age = scenarios$min.age[i],max.age = scenarios$max.age[i])
    # do the calculation 
    result.data <- (Partitioned.EvoHeritage(input.tree = current.EH.tree, num.repeats = number.repeats))
    # save the results
    write.csv(result.data , file = paste(results.file.path,"LivingFossil",str_pad(j, 4, pad = "0"),scenarios$name[i],sep="") , row.names = TRUE) 
    # print out the elapsed time is helpful 
    current.time <-  as.numeric(proc.time()[3])
    print(paste("time elapsed (hours)",(current.time-start.time)/3600,sep=" "))
  }
}

# I need to use standard units TODO
