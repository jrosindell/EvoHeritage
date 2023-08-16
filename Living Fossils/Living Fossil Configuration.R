# This code will configure the Living Fossil calculations example with study specific parameters and local paths
# Do not change the file name - it is sourced from the calculation scripts

# make a vector of paths to the tree files that will be used
tree.file.paths <- c()
for (i in 0:99) {
  tree.file.paths <- c(tree.file.paths,paste("Data/Mammal Phylogenies/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree",str_pad(i, 4, pad = "0"),".tre",sep=""))
}

# make a path for where to store all the results files
results.file.path <- "Data/Results/"
graphs.file.path <- "Data/Graphs/"
csv.file.path <- "Data/CSV/"
