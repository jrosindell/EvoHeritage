#Unit testing of general Functions for EvoHeritage

# clear workspace
rm(list=ls()) 

# source the functions we're wanting to test
source("./EvoHeritage Tools.R")

test.pass <- TRUE # we'll make this false if any test fails

# make a small test tree from scratch
make.test.tree <- function() {
  test.tree<- list()
  class(test.tree) <- "phylo"
  test.tree$edge <- matrix(c(
    6,7,
    7,1,
    7,8,
    8,2,
    8,3,
    6,9,
    9,4,
    9,5
  ),8,2,byrow=TRUE)
  test.tree$Nnode <- 4
  test.tree$tip.label <- c("A","B","C","D","E")
  test.tree$edge.length <- c(1,2,1,1,1,1,2,2)
  attr(test.tree,"order") <- "cladewise"
  return(test.tree)
}

# testing the vector.equal function which is used in the package and also used in the testing routine

if(vector.equal(c(0,0,0,0,0,3,2,1,2),c(0,0,0,0,0,3,2,1,2))) {
  print("PASSED - vector.equal function")
} else {
  print("WARNING vector.equal function FAILED unit testing")
  test.pass <- FALSE
}
if(vector.equal(c(300,2,1,2),c(300,2,1,2))) {
  print("PASSED - vector.equal function")
} else {
  print("WARNING vector.equal function FAILED unit testing")
  test.pass <- FALSE
}
if(!(vector.equal(c(300,2,1,2,1),c(300,2,1,2)))) {
  print("PASSED - vector.equal function")
} else {
  print("WARNING vector.equal function FAILED unit testing")
  test.pass <- FALSE
}
if(!(vector.equal(c(300,2,1,2),c(300,2,1,2,50)))) {
  print("PASSED - vector.equal function")
} else {
  print("WARNING vector.equal function FAILED unit testing")
  test.pass <- FALSE
}
if(!(vector.equal(c(301,2,1,2),c(300,2,1,2)))) {
  print("PASSED - vector.equal function")
} else {
  print("WARNING vector.equal function FAILED unit testing")
  test.pass <- FALSE
}
if(!(vector.equal(c(300,2,1,2),c(300,2,1,20)))) {
  print("PASSED - vector.equal function")
} else {
  print("WARNING vector.equal function FAILED unit testing")
  test.pass <- FALSE
}

# testing the get.node.ages function

test.result <- get.node.ages(make.test.tree())
if(vector.equal(test.result,c(0,0,0,0,0,3,2,1,2))) {
  print("PASSED - get.node.ages function")
} else {
  print("WARNING get.node.ages function FAILED unit testing")
  test.pass <- FALSE
}

# testing the add.stem function

test.result <- add.stem(make.test.tree(),5)
if(vector.equal(test.result$edge[,1],c(6,7,8,8,9,9,7,10,10))) {
  print("PASSED - add.stem function - edge test 1")
} else {
  print("WARNING add.stem function FAILED unit testing")
  test.pass <- FALSE
}
if(vector.equal(test.result$edge[,2],c(7,8,1,9,2,3,10,4,5))) {
  print("PASSED - add.stem function - edge test 2")
} else {
  print("WARNING add.stem function FAILED unit testing")
  test.pass <- FALSE
}
if(vector.equal(test.result$edge.length,c(5,1,2,1,1,1,1,2,2))) {
  print("PASSED - add.stem function - edge length test")
} else {
  print("WARNING add.stem function FAILED unit testing")
  test.pass <- FALSE
}
if(test.result$Nnode == 5) {
  print("PASSED - add.stem function - node count test")
} else {
  print("WARNING add.stem function FAILED unit testing")
  test.pass <- FALSE
}
if(vector.equal(get.node.ages(test.result),c(0,0,0,0,0,8,3,2,1,2))) {
  print("PASSED - add.stem function - node age test")
} else {
  print("WARNING add.stem function FAILED unit testing")
  test.pass <- FALSE
}

# testing the has.stem function

test.result <- has.stem(make.test.tree())
if(test.result == FALSE) {
  print("PASSED - has.stem function")
} else {
  print("WARNING has.stem function FAILED unit testing")
  test.pass <- FALSE
}

test.result <- has.stem(add.stem(make.test.tree(),1))
if(test.result == TRUE) {
  print("PASSED - has.stem function combined with add.stem function")
} else {
  print("WARNING has.stem function FAILED unit testing")
  test.pass <- FALSE
}

# testing the get.beta.vals function

test.rho <- 0.1
test.result <- get.beta.vals(make.test.tree(),test.rho)
beta1 <- exp(-test.rho*1)
beta2 <- exp(-test.rho*2)
correct.answer <- c(beta1,beta2,beta1,beta1,beta1,beta1,beta2,beta2)
if(vector.equal(correct.answer,test.result)) {
  print("PASSED - get.beta.vals - rho = 0.1")
} else {
  print("WARNING get.beta.vals FAILED unit testing")
  test.pass <- FALSE
}

test.rho <- 5
test.result <- get.beta.vals(make.test.tree(),test.rho)
beta1 <- exp(-test.rho*1)
beta2 <- exp(-test.rho*2)
correct.answer <- c(beta1,beta2,beta1,beta1,beta1,beta1,beta2,beta2)
if(vector.equal(correct.answer,test.result)) {
  print("PASSED - get.beta.vals - rho = 5")
} else {
  print("WARNING get.beta.vals FAILED unit testing")
  test.pass <- FALSE
}

# testing the alpha.clac function

# here I'm not going to use the tree as a test
# Instead, i'll challenge the function with a range of scenarios that I also calculate by hand here

test.result <- alpha.clac(start.node.age=6,end.node.age=4,rho=0.4,lambda=2.5,min.age=0,max.age=7)
correct.answer <- (2.5/0.4)*(1-exp(-0.4*2))
if(test.result == correct.answer) {
  print("PASSED - alpha.clac no cuts")
} else {
  print("WARNING alpha.clac FAILED unit testing")
  test.pass <- FALSE
}
  
test.result <- alpha.clac(start.node.age=6,end.node.age=4,rho=0.3,lambda=2.2,min.age=0,max.age=5)
correct.answer <- (2.2/0.3)*(1-exp(-0.3*1))
if(test.result == correct.answer) {
  print("PASSED - alpha.clac cut edge at older end")
} else {
  print("WARNING alpha.clac FAILED unit testing")
  test.pass <- FALSE
}

test.result <- alpha.clac(start.node.age=6,end.node.age=4,rho=0.3,lambda=2.2,min.age=0,max.age=3)
correct.answer <- 0
if(test.result == correct.answer) {
  print("PASSED - alpha.clac no overlap old edge")
} else {
  print("WARNING alpha.clac FAILED unit testing")
  test.pass <- FALSE
}

test.result <- alpha.clac(start.node.age=6,end.node.age=4,rho=0.3,lambda=2.2,min.age=7,max.age=8)
correct.answer <- 0
if(test.result == correct.answer) {
  print("PASSED - alpha.clac no overlap young edge")
} else {
  print("WARNING alpha.clac FAILED unit testing")
  test.pass <- FALSE
}

test.result <- alpha.clac(start.node.age=6,end.node.age=4,rho=0.3,lambda=2.2,min.age=4.1,max.age=8)
correct.answer <- ((2.2/0.3)*(1-exp(-0.3*1.9)))*exp(-0.3*0.1)
if((test.result - correct.answer)< correct.answer*(10^-10)) {
  print("PASSED - alpha.clac cut edge at younger end")
} else {
  print("WARNING alpha.clac FAILED unit testing")
  test.pass <- FALSE
}

test.result <- alpha.clac(start.node.age=6,end.node.age=4,rho=0.3,lambda=2.2,min.age=4.1,max.age=5)
correct.answer <- ((2.2/0.3)*(1-exp(-0.3*0.9)))*exp(-0.3*0.1)
if((test.result - correct.answer)^2 < correct.answer*(10^-20)) {
  print("PASSED - alpha.clac cut edge at older and younger ends")
} else {
  print("WARNING alpha.clac FAILED unit testing")
  test.pass <- FALSE
}

# testing the get.descendent.edges function

test.result <- get.descendent.edges(make.test.tree())
correct.answer <- c(1,2,3,4,5,6,7,8)
if(vector.equal(correct.answer,test.result[[1]])) {
  print("PASSED - get.descendent.edges - from")
} else {
  print("WARNING get.descendent.edges FAILED unit testing")
  test.pass <- FALSE
}
correct.answer <- c(5,2,5,4,5,8,7,8)
if(vector.equal(correct.answer,test.result[[2]])) {
  print("PASSED - get.descendent.edges - to")
} else {
  print("WARNING get.descendent.edges FAILED unit testing")
  test.pass <- FALSE
}

# testing the make.ancestral.evoheritage.tree function and the check.ancestral.evoheritage.tree function

if(check.ancestral.evoheritage.tree(make.ancestral.evoheritage.tree(make.test.tree(),rho=0.02,lambda = 1.2,min.age = 0,max.age = 4025))){
  print("PASSED - make.ancestral.evoheritage.tree and check.ancestral.evoheritage.tree 1")
}  else {
  print("WARNING make.ancestral.evoheritage.tree or check.ancestral.evoheritage.tree FAILED unit testing")
  test.pass <- FALSE
}
if(check.ancestral.evoheritage.tree(make.ancestral.evoheritage.tree(make.test.tree(),rho=0.02,lambda = 1.2,min.age = 0,max.age = 20))){
  print("PASSED - make.ancestral.evoheritage.tree and check.ancestral.evoheritage.tree 2")
}  else {
  print("WARNING make.ancestral.evoheritage.tree or check.ancestral.evoheritage.tree FAILED unit testing")
  test.pass <- FALSE
}
if(check.ancestral.evoheritage.tree(make.ancestral.evoheritage.tree(make.test.tree(),rho=0.02,lambda = 2.9,min.age = 5,max.age = 20))){
  print("PASSED - make.ancestral.evoheritage.tree and check.ancestral.evoheritage.tree 3")
}  else {
  print("WARNING make.ancestral.evoheritage.tree or check.ancestral.evoheritage.tree FAILED unit testing")
  test.pass <- FALSE
}
if(check.ancestral.evoheritage.tree(make.ancestral.evoheritage.tree(make.test.tree(),rho=0.05,lambda = 2.9,min.age = 1,max.age = 20))){
  print("PASSED - make.ancestral.evoheritage.tree and check.ancestral.evoheritage.tree 4")
}  else {
  print("WARNING make.ancestral.evoheritage.tree or check.ancestral.evoheritage.tree FAILED unit testing")
  test.pass <- FALSE
}

# test Random.EvoHeritage.copies function
# this is a tough test because it's a stochastic function
# My apporach is to code the function again from scratch, but using a for loop (instead of a more efficient while loop)
# and making sure the random numbers get used in the exact same way as in the main function
# I'll make this easier by writing a function to do much of the work 

test.random.copies <- function(test.EH.tree, edge.index, test.seed) {
  
  set.seed(test.seed)
  current.result <- rep(0,(length(test.EH.tree$tip.label)+test.EH.tree$Nnode))  
  current.result[test.EH.tree$edge[edge.index,2]] <- 1
  for (i in 2:length(test.EH.tree$edge[,1]))  {
    if (current.result[test.EH.tree$edge[i,1]]>0) {
      if (runif(1) <= test.EH.tree$beta[i]) {
        current.result[test.EH.tree$edge[i,2]] <- 1
      }
    }
  }
  current.random <- runif(1)
  
  set.seed(test.seed)
  current.result.b <- (Random.EvoHeritage.copies(test.EH.tree,edge.index))
  current.random.b <- runif(1)
  
  if(!vector.equal(current.result.b,current.result)) {
    print(c(current.result.b,"!=",current.result))
    return(FALSE)
  }
  if (current.random.b != current.random) {
    print(c(current.random.b,"!=",current.random))
    return(FALSE)
  }
  return(TRUE)
}

test.tree <- make.test.tree()
test.EH.tree <- make.ancestral.evoheritage.tree(test.tree,rho=0.01,lambda=1,min.age = 0,max.age = 4025)

if(test.random.copies(test.EH.tree,edge.index = 1,test.seed = 1)){
  print("PASSED - Random.EvoHeritage.copies 1")
}  else {
  print("WARNING Random.EvoHeritage.copies FAILED unit testing")
  test.pass <- FALSE
}
if(test.random.copies(test.EH.tree,edge.index = 3,test.seed = 1)){
  print("PASSED - Random.EvoHeritage.copies 2")
}  else {
  print("WARNING Random.EvoHeritage.copies FAILED unit testing")
  test.pass <- FALSE
}
if(test.random.copies(test.EH.tree,edge.index = 5,test.seed = 1)){
  print("PASSED - Random.EvoHeritage.copies 3")
}  else {
  print("WARNING Random.EvoHeritage.copies FAILED unit testing")
  test.pass <- FALSE
}

test.EH.tree <- make.ancestral.evoheritage.tree(test.tree,rho=1,lambda=2,min.age = 1,max.age = 7)
if(test.random.copies(test.EH.tree,edge.index = 1,test.seed = 1)){
  print("PASSED - Random.EvoHeritage.copies 4")
}  else {
  print("WARNING Random.EvoHeritage.copies FAILED unit testing")
  test.pass <- FALSE
}
if(test.random.copies(test.EH.tree,edge.index = 3,test.seed = 1)){
  print("PASSED - Random.EvoHeritage.copies 5")
}  else {
  print("WARNING Random.EvoHeritage.copies FAILED unit testing")
  test.pass <- FALSE
}
if(test.random.copies(test.EH.tree,edge.index = 5,test.seed = 1)){
  print("PASSED - Random.EvoHeritage.copies 6")
}  else {
  print("WARNING Random.EvoHeritage.copies FAILED unit testing")
  test.pass <- FALSE
}

# test Partitioned.EvoHeritage function 
# approach is to make another function to perform the same calculation and compare the results
test.Partitioned.EvoHeritage <- function(test.EH.tree,num.repeats,test.seed) {
  set.seed(test.seed)
  test.partitioned.result <- Partitioned.EvoHeritage(test.EH.tree,num.repeats)
  set.seed(test.seed)
  
  partitioned.result <- rep(0,(length(test.EH.tree$tip.label)+test.EH.tree$Nnode))
  total.result <- rep(0,(length(test.EH.tree$tip.label)+test.EH.tree$Nnode))
  for (i in 1: num.repeats) {
    for (j in 1:length(test.EH.tree$edge[,1])) {
      if (test.EH.tree$alpha[j] > 0) {
        temp.result <- (Random.EvoHeritage.copies(test.EH.tree,edge.index=j))
        num.copies <- sum(temp.result[1:(length(test.EH.tree$tip.label))])
        if(num.copies > 0 ) {
          partitioned.result <- partitioned.result + test.EH.tree$alpha[j]*(temp.result/num.copies)
          total.result <- total.result + test.EH.tree$alpha[j]*(temp.result)
        }
      } 
    }
  }
  partitioned.result <- partitioned.result/num.repeats
  total.result <- total.result/num.repeats
  
  # my strategy will be to check each tip methodically
  # this will not be fast but I want it to be a different approach to the function being tested
  for (i in 1:(length(test.EH.tree$tip.label))) {
    j <- which(test.partitioned.result$tip.label == test.EH.tree$tip.label[i])
    if(test.partitioned.result$partitioned.EvoHeritage[j] != partitioned.result[i]) {
      return(FALSE)
    }
    if(test.partitioned.result$ancestral.EvoHeritage[j] != total.result[i]) {
      return(FALSE)
    }
  }
  # A separate testing routine to test the ordering is correct 
  if((length(test.EH.tree$tip.label)) >= 2) {
    for(i in 2:(length(test.EH.tree$tip.label))) {
      if (test.partitioned.result$partitioned.EvoHeritage[i-1] < test.partitioned.result$partitioned.EvoHeritage[i]) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

test.EH.tree <- make.ancestral.evoheritage.tree(test.tree,rho=0.015,lambda=1.5,min.age = 0.1,max.age = 70)
if(test.Partitioned.EvoHeritage(test.EH.tree,1000,42)){
  print("PASSED - Final partitioning calculation 1")
}  else {
  print("WARNING Final partitioning calculation FAILED unit testing")
  test.pass <- FALSE
}
test.EH.tree <- make.ancestral.evoheritage.tree(test.tree,rho=0.03,lambda=10,min.age = 3,max.age = 4)
if(test.Partitioned.EvoHeritage(test.EH.tree,500,43)){
  print("PASSED - Final partitioning calculation 2")
}  else {
  print("WARNING Final partitioning calculation FAILED unit testing")
  test.pass <- FALSE
}

print("===================================")
print("Testing is completed")

if (test.pass) {
  print("PASSED - All tests")
} else {
  print("FAIL - There were some issues - check the output above")
}
print("===================================")
