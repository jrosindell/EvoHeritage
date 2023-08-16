#GPD calculations V1.2
require(caper)

# store the date of the origin of life in millions of years ago for future use
origin.life <- (3.77+4.28)*10^3/2

# this function returns the crown age of an ultrametric phylogeny
# phylogeny.in - a phylogeny object that should be ultrametric with edge lengths in millions of years
crown.age <- function(phylogeny.in, start_leaf = 1)  {
  # first check the data is in a sensible format
  # first check the data is in a sensible format (I'm suppressing this temporarily because the pez sub object doesn't meet the condition)
  #if (class(phylogeny.in) != "phylo") {
  #  stop("a phylogeny is needed to calculate Phi_rho.")
  #} else {
  #  if(!is.ultrametric(phylogeny.in)){
  #    warning("phylogeny was not ultrametric.")
  #  }
  #}
  
  cm.in <- clade.matrix(phylogeny.in)
  # convert the phylogeny to a clade matrix
  
  # order the clade matrix so that any edges starting with a node come after the edges that end with the same node
  # this enables us to work through the edge list methodically after the re ordering.
  looping <- TRUE
  while(looping) {
    looping <- FALSE # we will now test everything is in order
    for (i in 1: (dim(cm.in$edge)[1]-1)) {
      for (j in (i+1): dim(cm.in$edge)[1]) {
        from.edge <- cm.in$edge[i,2] # edge is connecting from this vertex
        to.edge <- cm.in$edge[j,1] # edge to connecting to this vertex
        if(from.edge == to.edge) {
          looping <- TRUE # we found something not right, so fix it and set looping = TRUE to test again
          temp.edge <- cm.in$edge[i,]
          cm.in$edge[i,] <- cm.in$edge[j,]
          cm.in$edge[j,] <- temp.edge
        }
      }
    }
  }
  
  # find the crown age - start by getting a list of parent nodes
  nodes.to.root <- c(start_leaf)
  current.node <- start_leaf
  for (i in 1:(dim(cm.in$edge)[1])) {
    from.edge <- cm.in$edge[i,2]
    to.edge <- cm.in$edge[i,1]
    if(from.edge == current.node){
      current.node <- to.edge
      nodes.to.root <- c(nodes.to.root,to.edge)
    }
  }
  # now add up all the edge lengths except for the stem
  crown.age <- 0
  for (i in nodes.to.root)
  {
    if (i != (dim(cm.in$clade.matrix)[2]+1)) {
      # this if statement excludes the stem
      crown.age <- crown.age + cm.in$edge.length[i]
    }
  }
  
  return(crown.age)
}
  
# this function return Phi_rho for a given phylogeny and value of n
# phylogeny.in - a phylogeny object that should be ultrametric with edge lengths in millions of years
# tips.in - a vector of the tip node numbers corresponding to the set we want to calculate Phi_rho for 
# rho - the value of the parameter rho of information erosion.... large numbers converge to richness, small numbers to PD with origin of life included
# origin.life - the date in millions of years ago of the origin of life
# std.units - a flag to tell if we want the answer in standard units or not
Phi_rho <- function(phylogeny.in , tips.in , rho , origin.life = (3.77+4.28)*10^3/2 , std.units = TRUE)  {
  
  # first check the data is in a sensible format
  # (I'm suppressing this temporarily because the pez sub object doesn't meet the condition)
  #if (class(phylogeny.in) != "phylo") {
  #  stop("a phylogeny is needed to calculate Phi_rho.")
  #} else {
  #  if(!is.ultrametric(phylogeny.in)){
  #    warning("phylogeny was not ultrametric.")
  #  }
  #}

  # calculate the value of rho from the value of n
  cm.in <- clade.matrix(phylogeny.in)
  # convert the phylogeny to a clade matrix

  # order the clade matrix so that any edges starting with a node come after the edges that end with the same node
  # this enables us to work through the edge list methodically after the re ordering.
  # note that this step is probably entirely redundant if the object is ordered 'cladewise'
  looping <- TRUE
  while(looping) {
    looping <- FALSE # we will now test everything is in order
    for (i in 1: (dim(cm.in$edge)[1]-1)) {
      for (j in (i+1): dim(cm.in$edge)[1]) {
        from.edge <- cm.in$edge[i,2] # edge is connecting from this vertex
        to.edge <- cm.in$edge[j,1] # edge to connecting to this vertex
        if(from.edge == to.edge) {
          looping <- TRUE # we found something not right, so fix it and set looping = TRUE to test again
          temp.edge <- cm.in$edge[i,]
          cm.in$edge[i,] <- cm.in$edge[j,]
          cm.in$edge[j,] <- temp.edge
        }
      }
    }
  }
  
  # find the crown age - start by getting a list of parent nodes
  nodes.to.root <- c(1)
  current.node <- 1
  for (i in 1:(dim(cm.in$edge)[1])) {
    from.edge <- cm.in$edge[i,2]
    to.edge <- cm.in$edge[i,1]
    if(from.edge == current.node){
      current.node <- to.edge
      nodes.to.root <- c(nodes.to.root,to.edge)
    }
  }
  # now add up all the edge lengths except for the stem
  crown.age <- 0
  for (i in nodes.to.root)
  {
    if (i != (dim(cm.in$clade.matrix)[2]+1)) {
      # this if statement excludes the stem
      crown.age <- crown.age + cm.in$edge.length[i]
    }
  }
  
  #now we make the stem reach back to the origin of life.
  if (!is.null(origin.life)){
  cm.in$edge.length[dim(cm.in$clade.matrix)[2]+1] = (origin.life - crown.age)
  }
  
  alpha.calc <- function(edge.length) {
    # precalculate rho * edge.length for efficiency
    rho.edge.length <- rho * edge.length
    if (rho.edge.length > 10^(-5)) {
      # use exponential function as no danger of numerical errors
      return((1-exp(-rho.edge.length))/rho)
    } else {
      # we are close to evaluating exp(0) so use polynomial expansion to avoid numerical errors
      return(edge.length*(1 - rho.edge.length /2 + rho.edge.length^2 /6 - rho.edge.length^3 /24 + rho.edge.length^4 /120)) 
    }
  }
  
  alpha <- sapply(cm.in$edge.length,alpha.calc)
  # accumulation along a branch
  beta <- exp(-rho * cm.in$edge.length)
  # erosion along a branch
  
  prop.eroded <- seq(1,1,length.out=length(cm.in$edge.length)) 
  # will give the proportion of information at this vertex that will be eroded
  prop.survive <- seq(-1,-1,length.out=length(cm.in$edge.length))
  # will give the proportion of information at this vertex that will survive erosion
  
  # calculate prop.eroded and prop.survive by looping over all edges
  # note that the ordering of edges is what makes this possible
  for(i in tips.in) {
    prop.eroded[i] <- 0
  }
  for (i in 1:(dim(cm.in$edge)[1])) {
    from.edge <- cm.in$edge[i,2]
    to.edge <- cm.in$edge[i,1]
    prop.survive[from.edge] <- 1-prop.eroded[from.edge]
    prop.eroded[to.edge] <- prop.eroded[to.edge]*(1- prop.survive[from.edge]*beta[from.edge])
  }
  # we have to do the root also
  to.edge <- cm.in$edge[dim(cm.in$edge)[1],1]
  prop.survive[to.edge] <- 1-prop.eroded[to.edge]
  
  # total up the GPD
  GPD_tot <- 0
  for (i in 1:(length(cm.in$edge.length))) {
    GPD_tot <- GPD_tot + (alpha[[i]]* prop.survive[i])
  }
  
  if (std.units) {
    # return Phi_rho in standard units
    std.unit <- alpha.calc(1) # this is now for 1 million years which is a new standard unit
    return(GPD_tot / std.unit)
  } else {
    # return Phi_rho without standardisation
    return(GPD_tot) 
  }
}

# now let's try out the function to check that it behaves as expected

Test_GPD <- function() {
  
  # load in the test data from caper
  data(BritishBirds)

  # rho small WITHOUT unit standardisation gives us PD with added stem to the origin of life
  print(c("Phi_rho = " , Phi_rho(BritishBirds.tree,1:250,0,std.units = FALSE)))
  PD <- pd.calc(clade.matrix(BritishBirds.tree),tip.subset = 1:250,root.edge=FALSE)
  print(c("PD = " , PD))
  print(c("PD with stem = " , (PD + origin.life - crown.age(BritishBirds.tree))))      

  # Now standard units still give us PD correctly
  print(c("Phi_rho = " , Phi_rho(BritishBirds.tree,5:200,0,std.units = TRUE)))
  PD <- pd.calc(clade.matrix(BritishBirds.tree),tip.subset = 5:200,root.edge=FALSE)
  print(c("PD = " , PD))
  print(c("PD with stem = " , (PD + origin.life - crown.age(BritishBirds.tree))))   

  # rho large WITH unit standardisation gives us species richness
  print(c("Phi_rho = " , Phi_rho(BritishBirds.tree,1:250,10,std.units = TRUE)))
  print(c("species richness = " , length(1:250)))  

  print(c("Phi_rho = " , Phi_rho(BritishBirds.tree,5:200,10,std.units = TRUE)))
  print(c("species richness = " , length(5:200)))   

  # small rho WITH or WITHOUT unit standardisation for a single species gives us the origin of life
  print(c("Phi_rho = " , Phi_rho(BritishBirds.tree,c(42),0,std.units = FALSE)))
  print(c("origin of life (Millions of years) = " , origin.life))

  # LARGE POSITIVE n WITH unit standardisation for a single species gives us one
  print(c("Phi_rho = " , Phi_rho(BritishBirds.tree,c(42,43),0,std.units = TRUE)))

  print(c("crown age",crown.age(BritishBirds.tree)))
  print(c("crown age",crown.age(BritishBirds.tree,1)))
  print(c("crown age",crown.age(BritishBirds.tree,20)))

  #so, it all works, though note that it might not be as fast as would be ideal it's a start

  }

# now I will write a function gpd_pez which will perform like the .pd function of pez and consider a community matrix
# x should be a comparative.comm object as defined in pez
# I'm not trying to make the code efficient I'm just trying to get it to work at this stage
gpdpez <- function(x, rho , origin.life = NULL , std.units = TRUE) {
  number.communities <- nrow(x$comm)
  # pez returns some kind of dataframe, I'll just return a vector of values
  gpd.result <- c()
  for ( i in 1:number.communities) {
    current.tips <- as.numeric(which(x$comm[i,]>=1)) # get a list of the tips present in this community
    current.phy <- x$phy 
    class(current.phy) <- "phylo" # dangerous hack thing to do but annoyingly the x object phylo doesn't have a class
    current.gpd <- Phi_rho(current.phy, current.tips , rho , origin.life , std.units)
    gpd.result <- c(gpd.result,current.gpd)
  }
  return(gpd.result)
}


# now I will write a function gpd_pez which will perform like the .pd function of pez and consider a community matrix
# x should be a comparative.comm object as defined in pez
# I'm not trying to make the code efficient I'm just trying to get it to work at this stage
richpez <- function(x, rho , origin.life = (3.77+4.28)*10^3/2 , std.units = TRUE) {
  number.communities <- nrow(x$comm)
  # pez returns some kind of dataframe, I'll just return a vector of values
  gpd.result <- c()
  for ( i in 1:number.communities) {
    current.tips <- as.numeric(which(x$comm[i,]>=1)) # get a list of the tips present in this community
    gpd.result <- c(gpd.result,length(current.tips))
  }
  return(gpd.result)
}