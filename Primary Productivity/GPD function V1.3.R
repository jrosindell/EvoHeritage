

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
