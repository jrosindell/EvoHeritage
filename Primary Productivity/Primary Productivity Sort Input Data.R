#Code for gathering input data for this analysis

# clear workspace
rm(list=ls()) 

# keep paths and various other details in a separate config file
source("./Primary Productivity Configuration.R") 

# Headers
library(pez)

# Load raw data
data <- read.delim("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-cdr.273.10&entityid=8cbc85eae79f460bafd377ffc940f0f8")
tree <- read.tree(text="((Magnolia_grandiflora:0.04749,Amborella_trichopoda:0.141321)0.877000:0.0080285,((((Bouteloua_gracilis:0.055299,Schizachyrium_scoparium:0.179071)0.058000:0.013591,(Panicum_virgatum:0.043921,((Sorghastrum_nutans:0.001392,Andropogon_gerardii:0.004429)0.826000:0.009624,((Sporobolus_cryptandrus:0.043752,Buchloe_dactyloides:0.030027)0.863000:0.01107,Poa_pratensis:0.096489)0.951000:0.018588)0.947000:0.024782)0.953000:0.022698)0.925000:0.029753,((Elymus_canadensis:0.00552,Pascopyrum_smithii:0.010382)0.988000:0.030761,Koeleria_cristata:0.023313)0.921000:0.039637)1.000000:0.268846,(Anemone_cylindrica:0.180639,((((Salvia_officinalis:0.005257,Monarda_fistulosa:0.012795)0.974000:0.051126,(Asclepias_tuberosa:0.065927,Asclepias_incarnata:0.002697)1.000000:0.189149)0.942000:0.061359,((((Liatris_aspera:0.064393,Coreopsis_palmata:0.037974,Rudbeckia_hirta:0.05118)0.998000:0.025406,(Achillea_millefolium:0.043562,Lactuca_sativa:0.028824)0.956000:0.023085)0.807000:0.007382,(Symphyotrichum_oolentangiense:0.054735,(Solidago_nemoralis:0.0,Oligoneuron_rigidum:0.00162)0.927000:0.015239)0.968000:0.028916)0.777000:0.006146,Symphyotrichum_cordifolium:0.020208)1.000000:0.077771)0.843000:0.04181,(Euphorbia_pubentissima:0.125511,(((Quercus_serrata:0.002244,Quercus_ellipsoidalis:0.0068)0.765000:0.036008,Quercus_macrocarpa:0.0)0.906000:0.077946,((Lupinus_perennis:0.058336,Lespedeza_capitata:0.113526)0.306000:0.005766,((Vicia_villosa:0.074364,Astragalus_canadensis:0.045926)1.000000:0.04969,(Amorpha_canadesis:0.023175,(Dalea_purpurea:0.05049,Dalea_candida:0.055157):0.05678)1.000000:0.019707)0.993000:0.02651)0.946000:0.073925)1.000000:0.015029)0.740000:0.037682)0.999000:0.039026)0.988000:0.019496)0.922000:0.0080285);") # Taken from supplementary materials of paper

# NOTE: chronogram (and so NPRS) is no longer supported in ape, and instead of back-porting it in from v.2.5-1 (which I tried and is slightly tricky because it depends on linked C code) I'm just using chronos (and so the newer Sanderson method) instead
tree <- chronos(tree)

# we can test the crown age with
# class(tree) <- "phylo"
# crown.age(tree)
# and see that it is 1 (note this function requires sourcing of the EvoHeritage tools if the reader wants to test)
# when we want it to be 234.5 million years ago
# so as the tree is already ultrametric let's multiply all edges accordingly

tree$edge.length <- tree$edge.length*234.5

# Community biomass and composition data
data$site.year <- with(data, paste(Plot, Year, sep="_"))
biomass <- with(data, tapply(Biomass..g.m2., site.year, sum))
comm <- data[,c(18:35,38)]
comm <- comm[!duplicated(comm$site.year),]
rownames(comm) <- comm$site.year
comm$site.year <- NULL
comm <- as.matrix(comm)
colnames(comm) <- c("Achillea_millefolium","Agropyron_smithii","Amorpha_canescens","Andropogon_gerardii","Asclepias_tuberosa","Elymus_canadensis","Koeleria_cristata","Lespedeza_capitata","Liatris_aspera","Lupinus_perennis","Monarda_fistulosa","Panicum_virgatum","Petalostemum_purpureum","Poa_pratensis","Quercus_ellipsoidalis","Quercus_macrocarpa","Schizachyrium_scoparium","Sorghastrum_nutans")
biomass <- biomass[rownames(comm)]

# Phylogenetic diversity and fix names on phylogeny (taxonomy updates)
tree$tip.label[tree$tip.label=="Amorpha_canadesis"] <- "Amorpha_canescens" # Seems to be a typo
tree$tip.label[tree$tip.label=="Dalea_purpurea"] <- "Petalostemum_purpureum" # Old name @ community data (not a phylogeny problem)
tree$tip.label[tree$tip.label=="Pascopyrum_smithii"] <- "Agropyron_smithii" # Old name @ community data (not a phylogeny problem)

c.data <- comparative.comm(tree, comm)

save(c.data,biomass,file = paste(inputs.file.path,"PProd.input.data.rda",sep=""))


