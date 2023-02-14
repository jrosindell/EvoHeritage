# This code will configure the Living Fossil calculations example with study specific parameters and local paths
# Do not change the file name - it is sourced from the calculation scripts and GitIgnored

# make a vector of paths to the tree files that will be used
tree.file.paths <- c()
for (i in 0:99) {
  tree.file.paths <- c(tree.file.paths,paste("Data/Mammal Phylogenies/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree",str_pad(i, 4, pad = "0"),".tre",sep=""))
}

# make a path for where to store all the results files
results.file.path <- "Data/Results/"

# make a data frame of scenarios
rho <- c(rep(0.01,length=6),rep(0.1,length=6)) # these are the two I did before
lambda <- rep(1,length=12)
min.age <- c(145,66,23,145,66,23)
max.age <- c(201.4,145,66,4025,4025,4025)
min.age <- c(min.age,min.age)
max.age <- c(max.age,max.age)
seed <- 1:12
name <- c("JurassicNew2","CretaceousNew2","PaleogeneNew2","JurassicAll2","CretaceousAll2","PaleogeneAll2","JurassicNew1","CretaceousNew1","PaleogeneNew1","JurassicAll1","CretaceousAll1","PaleogeneAll1")
scenarios <- data.frame(name,seed,rho,lambda,min.age,max.age)

number.repeats <- 10000# we'll make this 10,000 for the main run this is a test.