#General Functions for EvoHeritage Calculations

# install useful packages
require(ape)
require(stringr)

# make a quick test funciton to see that two vectors are the same
vector.equal <- function(v1,v2) {
  if (length(v1) == length(v2)) { # first check that the vectors are of the same length
    for (i in 1:(length(v1))) {
      if (v1[i] != v2[i]) { # compare each element of the vectors
        return(FALSE) # we've found an element that does not match up
      }
    }
    return(TRUE) # we've checked every element and everything was identical 
  } else {
    return(FALSE) # they are not so return false
  }
}

# I want a function that will sum across the tree (assuming it's ultrametric) and give me dates for each node based on branch lengths
get.node.ages <- function(input.tree) {
  # start by calculating parameters for the tree - this is probably unnecessary but it's useful for readability
  num.interior.nodes <- input.tree$Nnode
  num.leaf.nodes <- length(input.tree$tip.label)
  num.nodes <- num.interior.nodes + num.leaf.nodes 
  num.edges <- length(input.tree$edge.length)
  # store a vector for the node ages starting with 0 everywhere and change later
  node.ages <- rep(0,length=num.nodes) 
  # loop over all edges starting from terminal edges
  for (i in num.edges:1) {
    if (node.ages[input.tree$edge[i,1]] == 0) {
      # this part is not already filled out from another route up the tree (we assume the tree is ultra metric)
      node.ages[input.tree$edge[i,1]] <- input.tree$edge.length[i] + node.ages[input.tree$edge[i,2]]
      # add branch length to total so far and save the result in the correct place of node.ages
    }
  }
  return(node.ages)
}

# add a stem to a given tree
add.stem <- function(input.tree,stem.length) {
  # store number of edges and leaf nodes for convenience
  num.edges <- length(input.tree$edge.length)
  num.leaf.nodes <- length(input.tree$tip.label)
  # increase the number of nodes
  input.tree$Nnode <- input.tree$Nnode+1
  # loop over edges and change node pointers
  for (i in 1:num.edges) {
    # the first node of an edge is never a leaf so needs to be incremented to make space for a new root 
    input.tree$edge[i,1] <- input.tree$edge[i,1]+1
    # the second node of an edge needs to be checked, if it's a leaf we leave it unchanged
    if ( input.tree$edge[i,2] > num.leaf.nodes)
      # otherwise we increment it like the rest
      input.tree$edge[i,2] <- input.tree$edge[i,2] +1
  }
  # add a new length for the stem
  input.tree$edge.length <- c(stem.length,input.tree$edge.length)
  # add a new edge to the matrix of edges
  input.tree$edge <- rbind(c(num.leaf.nodes+1,num.leaf.nodes+2),input.tree$edge)
  return(input.tree)
}

# return true if the tree has a stem otherwise return false
has.stem <- function(input.tree) {
  # get the root node number
  root.node <- input.tree$edge[1,1]
  # look to see how many edges start with that root
  if(length(which(input.tree$edge[, 1] == root.node))==1) {
    # only one so it's a stem
    return(TRUE)
  } else {
    # not one (and zero or negative is impossible) so the tree does not have a stem
    return(FALSE)
  }
}

# return a vector of the beta values for a tree given a value of rho
get.beta.vals <- function(input.tree,rho) {
  # I'll assume rho is already per unit length with lengths likely in millions of years
  beta <- exp(-rho * input.tree$edge.length)
  # beta will automatically be a vector because input.tree$edge.length is a vector and R is clever like that
  return(beta) 
}

# calculate net alpha for a single edge
# start.node.age is the node at the start of the edge
# end.node.age is the node at the end of the edge
# min.age gives the minimum number of years ago that EvoHeritage accumulation is turned on for 
# max.age gives the maximum number of years ago that EvoHeritage accumulation is turned on for
# so we can study the fate of information generated between these (min.age, max.age) dates and not more.
alpha.clac <- function(start.node.age,end.node.age,rho,lambda,min.age,max.age) {
  # if the edge doesn't overlap at all with the period of interest alpha will be zero
  if (end.node.age > max.age) {
    return(0)
  }
  if (start.node.age < min.age) {
    return(0)
  }
  # find out when accumulation starts and ends given the restrictions
  accumulation.start <- min(start.node.age,max.age)
  accumulation.end <- max(min.age,end.node.age)
  # calculate an alpha for the given period
 
  alpha <- (lambda/rho)*(1-exp(-rho*(accumulation.start-accumulation.end)))
  # find out if the edge continues after accumulation stops
  if (end.node.age < min.age) {
    # need to adjust for this with attrition on the remainder of the edge
    alpha <- alpha * exp(-rho * (min.age-end.node.age))
  } 
  return(alpha)
}

# I want to build a function that will return me the range of edges descended from every edge of a tree
# These will become pre-calculations to reduce the side of program loops later on
# It requires cladewise ordering which guarantees that any clade is represented by a fixed range of edges
get.descendent.edges <- function(input.tree) {
  if (attributes(input.tree)$order != "cladewise") {
    stop("your tree needs to be ordered cladewise for this function to work")
  } else {
    # it's convenient to pre calculate dimension of the tree of sizing loops and vectors in a moment
    num.edges <- length(input.tree$edge.length)
    num.leaf.nodes <- length(input.tree$tip.label)
    # set up the output parameters - each edge must include itself as a minimum
    descendant.edge.from <- 1:num.edges
    descendant.edge.to <- 1:num.edges
    # loop over the edges starting with the tips and working up
    for (i in num.edges:1) {
      ancestor.node <- input.tree$edge[i,1] # this is the ancestor node
      ancestor.edge <- which(input.tree$edge[, 2] == ancestor.node) # this is the edge that connects to the ancestor node
      if (length(ancestor.edge)> 0) {
        # there is an edge above us update its from and to values
        # the from value is the smaller one so if it's bigger than something meant to be in range then we update it
        if ( descendant.edge.from[ancestor.edge] > descendant.edge.from[i]) { 
          descendant.edge.from[ancestor.edge] = descendant.edge.from[i] 
        }
        # the to value is the bigger one so if it's small than something meant to be in range then we update it
        if ( descendant.edge.to[ancestor.edge] < descendant.edge.to[i]) { 
          descendant.edge.to[ancestor.edge] = descendant.edge.to[i] 
        }
      }
    }
    return(list(descendant.edge.from,descendant.edge.to))
  }
}

# SPECIFICATIONS OF AN ANCESTRAL EVOHERITAGE TREE
# a tree stored in the normal way with the following additional conditions
# cladewise ordering
# has branch lengths measured in millions of years
# has a stem going back to the origin of life 4025 Mya
# has $tip.label defined
# to this we have to add the following properties that must be consisten with each other
# $node.ages - gives the age of all nodes in millions of years - so leaves have a value of 0 and the root should be 4025
# $alpha - net accumulation for each edge
# $beta - attrition for each edge
# $descendant.edge.from and $descendant.edge.to - giving the range of edges descended from the given edge
# $rho > 0 , $lambda > 0 - the parameters for accumulation and attrition that $alpha and $beta were calculated for
# $min.age > $max.age > 0 - the period of history during which accumulation counts 

# make a function to check that a tree object is consistent and meets the definition of being an evoheritage tree
# this function can also be used for testing the function to make an evoheritage tree
check.ancestral.evoheritage.tree <- function(input.tree) {
  # first we'll do basic tests that need to stop further testing if they fail
  
  # check cladewise ordering
  if (attributes(input.tree)$order != "cladewise") { 
    print("not ordered cladewise")
    return(false)
  } 
  if(is.null(input.tree$rho)) {
    print("no rho was specified")
    return(false)
  }
  if(is.null(input.tree$lambda)) {
    print("no lambda was specified")
    return(false)
  }
  if(is.null(input.tree$min.age)) {
    print("no min.age was specified")
    return(false)
  }
  if(is.null(input.tree$max.age)) {
    print("no max.age was specified")
    return(false)
  }
  if(is.null(input.tree$tip.label)) {
    print("no tip.label was specified")
    return(false)
  }
  if(is.null(input.tree$node.ages)) {
    print("no node.ages was specified")
    return(false)
  }
  if(is.null(input.tree$alpha)) {
    print("no alpha was specified")
    return(false)
  }
  if(is.null(input.tree$beta)) {
    print("no beta was specified")
    return(false)
  }
  if(is.null(input.tree$descendant.edge.from)) {
    print("no descendant.edge.from was specified")
    return(false)
  }
  if(is.null(input.tree$descendant.edge.to)) {
    print("no descendant.edge.to was specified")
    return(false)
  }
  if(is.null(input.tree$standardised.units)) {
    print("no standardised.units were specified")
    return(false)
  }
 
  # now we'll do more advanced testing 
  is.EHT <- TRUE # we'll set this to false and print a message if anything doesn't check out.
  
  if(input.tree$rho < 0) {
    print("rho was negative")
    is.EHT <- FALSE
  }
  if(input.tree$lambda < 0) {
    print("lambda was negative")
    is.EHT <- FALSE
  }
  if(input.tree$min.age < 0) {
    print("min.age was negative")
    is.EHT <- FALSE
  }
  if(input.tree$max.age < 0) {
    print("max.age was negative")
    is.EHT <- FALSE
  }
  # check $node.ages
  test.node.ages <- get.node.ages(input.tree)
  if (!(vector.equal(test.node.ages,input.tree$node.ages))) {
    is.EHT <- FALSE
    print("node.ages was not set correctly")
  }
  # check has a stem going back to the origin of life 4025 Mya
  if (!(max(test.node.ages)==4025)) {
    is.EHT <- FALSE
    print("oldest node did not correspond to the origin of life")
  }
  # check standardised units
  if (!(input.tree$standardised.units==input.tree$lambda/input.tree$rho*(1-exp(-input.tree$rho)))) {
    is.EHT <- FALSE
    print("standardised units were wrong")
  }
    
  test.descendant.edge <- get.descendent.edges(input.tree)
  if (!(vector.equal(test.descendant.edge[[1]],input.tree$descendant.edge.from))) {
    is.EHT <- FALSE
    print("descendant.edge.from was not set correctly")
  }
  if (!(vector.equal(test.descendant.edge[[2]],input.tree$descendant.edge.to))) {
    is.EHT <- FALSE
    print("descendant.edge.to was not set correctly")
  }
  
  test.beta <- get.beta.vals(input.tree,input.tree$rho)
  if (!(vector.equal(test.beta,input.tree$beta))) {
    is.EHT <- FALSE
    print("beta was not set correctly")
  }
 test.alpha <- rep(0,length=length(test.beta))
  for (i in 1:length(test.alpha)) {
    test.alpha[i] <- alpha.clac(start.node.age = test.node.ages[input.tree$edge[i,1]],end.node.age = test.node.ages[input.tree$edge[i,2]],rho = input.tree$rho,lambda = input.tree$lambda,min.age = input.tree$min.age,max.age = input.tree$max.age)
  }
 if (!(vector.equal(test.alpha,input.tree$alpha))) {
   is.EHT <- FALSE
   print("alpha was not set correctly")
 } 

  return(is.EHT)
}

# make a function to process a tree by adding alpha and beta values to it.
make.ancestral.evoheritage.tree <- function(input.tree,rho,lambda = 1,min.age = 0,max.age = 4025) {
  # it's convenient to pre calculate dimension of the tree of sizing loops and vectors in a moment
  num.interior.nodes <- input.tree$Nnode
  num.leaf.nodes <- length(input.tree$tip.label)
  num.nodes <- num.interior.nodes + num.leaf.nodes 
  num.edges <- length(input.tree$edge.length)
  
  # we take the origin of life as the midpoint of current min and max estimates 
  origin.life <- 4025
  # (3.77*10^9,4.28*10^9) years and measure in millions of years 
  
  # calculate the node ages
  node.ages <- get.node.ages(input.tree)
  # figure out how much we need to add to the root of the tree to include the origin of life
  to.add <- origin.life-node.ages[num.leaf.nodes+1] # origin of life take the root age 
  if (!has.stem(input.tree)) {
    # we have to add a stem
    input.tree <- add.stem(input.tree,stem.length=to.add)
    # we also have to add to the node ages for the new node - do this and at the same time save the node.ages data to the tree object
    input.tree$node.ages <- c(node.ages[1:num.leaf.nodes],origin.life,node.ages[(num.leaf.nodes+1):num.nodes])
  } else {
    # add to the stem edge length
    input.tree$edge.length[1] <- input.tree$edge.length[1]+to.add
    # update the node ages accordingly
    node.ages[num.leaf.nodes+1] <- origin.life
    # save the node.ages data to the tree object
    input.tree$node.ages <- node.ages
  }
  # update the stored values where needed
  num.edges <- length(input.tree$edge.length)
  
  # calculate and save the alpha and beta vectors using the alpha.clac function
  # generate a vector of alpha values starting with zeros
  alpha <- rep(0,length=num.edges)
  # loop over all edges
  for (i in num.edges:1) {
    # calculate each one using the separate function
    alpha[i] <- alpha.clac(input.tree$node.ages[input.tree$edge[i,1]],input.tree$node.ages[input.tree$edge[i,2]],rho,lambda,min.age,max.age)
  }
  
  standardised.units <- alpha.clac(1,0,rho,lambda,0,1)
  
  input.tree$alpha <- alpha
  input.tree$beta <- get.beta.vals(input.tree,rho) 
  descendent.edges <- get.descendent.edges(input.tree)
  input.tree$descendant.edge.from <- descendent.edges[[1]]
  input.tree$descendant.edge.to <- descendent.edges[[2]]
  input.tree$rho <- rho
  input.tree$lambda <- lambda
  input.tree$min.age <- min.age
  input.tree$max.age <- max.age
  input.tree$standardised.units <- standardised.units
  return(input.tree)
}

# return a vector indicating which terminal vertices get copies of EvoHeritage
# from the edge indicated by edge.index and surviving attrition along all edges
# the input.tree does need to be an ancestral EvoHeritage tree made and verified by the functions
Random.EvoHeritage.copies <- function(input.tree,edge.index) {
  # it's convenient to pre calculate dimension of the tree of sizing loops and vectors in a moment
  num.interior.nodes <- input.tree$Nnode
  num.leaf.nodes <- length(input.tree$tip.label)
  num.nodes <- num.interior.nodes + num.leaf.nodes 
  num.edges <- length(input.tree$edge.length)
  
  current.result <- rep(0,length=num.nodes) # this will show if a node has a copy of the ancestral information or not
  current.result[input.tree$edge[edge.index,2]] <- 1 # the descendant node of the indicated edge gets the EvoHeritage to start off with
  
  i <- input.tree$descendant.edge.from[edge.index] # this index i will be used to check all descendant edges to see if they inherit 
  end.loop <- input.tree$descendant.edge.to[i] # we use descendant.edge.to to see the last edge with a chance of inherited this EvoHeritage
  i <- i+1 # we start from the next edge along because the indicated edge already has been handed in the initial conditions 
  while (i <= end.loop) { 
    # we use a while loop to make it easier to jump ahead if some edges are not needed to be evaluated
    if (current.result[input.tree$edge[i,1]] > 0) {
      if (runif(1) <= input.tree$beta[i]) {
        # information survives along this edge on this occasion
        current.result[input.tree$edge[i,2]] <- 1
      } else {
        # as soon as an edge loses EvoHeritage from attrition none of its descendants are worth evaluating 
        # because they won't have any chance of inheriting the EvoHeritage
        i <- input.tree$descendant.edge.to[i]
      }
    } else {
      # this edge doesn't have any EvoHeritage - its descendants are not worth evaluating 
      i <- input.tree$descendant.edge.to[i]
    }
    i <- i+1 # advance to the next edge noting we're using cladewise ordering
  }
  return(current.result)
}

# This function takes an EvoHeritage tree and produces a dataframe of partitioned EvoHeritage on each tip as a return
# loop over all edges with alpha > 0
# ripple down the tree looking to see if EvoHeritage survives based on stochastic draws
# this will then be done with the cladewise ordering in a single sweep for each edge with alpha > 0
# it uses a while loop rather than a for loop so that we can jump forward a big chunk at a time if nothing survives.
Partitioned.EvoHeritage <- function(input.tree,num.repeats) {
  if (attributes(input.tree)$order != "cladewise") {
    stop("your tree needs to be ordered cladewise for this function to work")
  } else {
    # because of the cladewise ordering I can work down the edge list in figuring out which tips inherit ancestral information
    # I'll pre-calculate the beta value for each edge as this can be used again with each repeat of the monty carlo method
    # so now I've calculated the beta values to use for all edges 
    # it's convenient to pre calculate dimension of the tree of sizing loops and vectors in a moment
    num.interior.nodes <- input.tree$Nnode
    num.leaf.nodes <- length(input.tree$tip.label)
    num.nodes <- num.interior.nodes + num.leaf.nodes 
    num.edges <- length(input.tree$edge.length)

    # initialise final result after all repeats
    partitioned.EvoHeritage <- rep(0,length=num.leaf.nodes) # this will store the final result for tips
    ancestral.EvoHeritage <- rep(0,length=num.leaf.nodes) # this will store the total ancestral information
    for (i in 1:num.repeats) {
      # here we're looping over the repeats
      # build a vector to store the interim result from the current repeat
      for (j in 1:num.edges) {
        if (input.tree$alpha[j] > 0) {
          current.result <- Random.EvoHeritage.copies(input.tree,edge.index=j)
          num.copies <- sum(current.result[1:num.leaf.nodes]) 
          if (num.copies > 0) {
            # avoid divide by zero in case no information survives anywhere
            partitioned.EvoHeritage <- partitioned.EvoHeritage + (current.result[1:num.leaf.nodes]/num.copies)*input.tree$alpha[j]
            ancestral.EvoHeritage <- ancestral.EvoHeritage + (current.result[1:num.leaf.nodes])*input.tree$alpha[j]
            # sum paritioned information on leaf nodes
          }
        }
      } # end of loop over number of edges
      
    } # close repeat calcs
    tip.label <- input.tree$tip.label # get tip labels
    partitioned.EvoHeritage <- partitioned.EvoHeritage/num.repeats/input.tree$standardised.units # normalise for number of repeats
    ancestral.EvoHeritage <- ancestral.EvoHeritage/num.repeats/input.tree$standardised.units  # normalise for number of repeats
    result.data <- data.frame(tip.label,partitioned.EvoHeritage,ancestral.EvoHeritage)
    result.data <- result.data[order(-result.data$partitioned.EvoHeritage),]
    result.data$lf.rank <- 1:num.leaf.nodes
    return(result.data) # return a data data frame with results
  } # close else
} # close function     
      

