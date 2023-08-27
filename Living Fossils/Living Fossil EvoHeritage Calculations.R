#Code for detecting living fossils on a phylogentic tree using partitioned EvoHeritage

# clear workspace
rm(list=ls()) 

# source the functions we're wanting to use
source("../EvoHeritage Tools.R")
source("./Living Fossil Configuration.R") # a separate file giving the paths and filenames

# make a data frame of scenarios this is te
min.age <- c(145,66,145,66,145,66,0,0)
max.age <- rep(4025,length=8)
lambda <- rep(1,length=8)
rho <- c(0.01,0.01,0.001,0.001,0,0,0,0.01)
number.repeats <- c(rep(10000,length=6),rep(100,length=2))
seed <- seq(from=0 ,to=7000,by=1000) # increasing by 1000 because will add an index for the tree as well
name <- c("JurassicAll2","CretaceousAll2","JurassicAll3","CretaceousAll3","JurassicAllED","CretaceousAllED","EDcomp","EDcomprho")
scenarios <- data.frame(name,seed,rho,lambda,min.age,max.age,number.repeats)

# record the time when the job started
start.time <- as.numeric(proc.time()[3])

for(i in 1:(length(scenarios[,1]))) {
  print("simulating new scenario")
  print(scenarios[i,])
  for (j in 1:length(tree.file.paths)) {
    print(paste("simulating tree",j))
    # read in current tree
    current.tree <- ape::read.tree(file=tree.file.paths[j]) 
    # set seed for repeat ability
    set.seed(scenarios$seed[i]+j) 
    # make an ancestral EvoHeritage Tree
    current.EH.tree <- make.ancestral.evoheritage.tree(current.tree,rho=scenarios$rho[i],lambda = scenarios$lambda[i],min.age = scenarios$min.age[i],max.age = scenarios$max.age[i])
    # do the calculation 
    result.data <- (Partitioned.EvoHeritage(input.tree = current.EH.tree, num.repeats = scenarios$number.repeats[i]))
    # save the results
    write.csv(result.data , file = paste(results.file.path,"LivingFossil",str_pad(j, 4, pad = "0"),scenarios$name[i],sep="") , row.names = TRUE) 
    # print out the elapsed time is helpful 
    current.time <-  as.numeric(proc.time()[3])
    print(paste("seed",scenarios$seed[i]+j,"time elapsed (hours)",(current.time-start.time)/3600,sep=" "))
  }
}
