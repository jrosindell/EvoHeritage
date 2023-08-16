# Repeating Cadotte et al. (2009) PLoS One for GPD
# Will Pearse - 2022-06-22

rm(list=ls())

# Headers
library(pez)

source("GPD function V1.3.R")
source("./Primary Productivity Configuration.R") # keep paths and various other details in a separate config file

# Load raw data
data <- read.delim("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-cdr.273.10&entityid=8cbc85eae79f460bafd377ffc940f0f8")
tree <- read.tree(text="((Magnolia_grandiflora:0.04749,Amborella_trichopoda:0.141321)0.877000:0.0080285,((((Bouteloua_gracilis:0.055299,Schizachyrium_scoparium:0.179071)0.058000:0.013591,(Panicum_virgatum:0.043921,((Sorghastrum_nutans:0.001392,Andropogon_gerardii:0.004429)0.826000:0.009624,((Sporobolus_cryptandrus:0.043752,Buchloe_dactyloides:0.030027)0.863000:0.01107,Poa_pratensis:0.096489)0.951000:0.018588)0.947000:0.024782)0.953000:0.022698)0.925000:0.029753,((Elymus_canadensis:0.00552,Pascopyrum_smithii:0.010382)0.988000:0.030761,Koeleria_cristata:0.023313)0.921000:0.039637)1.000000:0.268846,(Anemone_cylindrica:0.180639,((((Salvia_officinalis:0.005257,Monarda_fistulosa:0.012795)0.974000:0.051126,(Asclepias_tuberosa:0.065927,Asclepias_incarnata:0.002697)1.000000:0.189149)0.942000:0.061359,((((Liatris_aspera:0.064393,Coreopsis_palmata:0.037974,Rudbeckia_hirta:0.05118)0.998000:0.025406,(Achillea_millefolium:0.043562,Lactuca_sativa:0.028824)0.956000:0.023085)0.807000:0.007382,(Symphyotrichum_oolentangiense:0.054735,(Solidago_nemoralis:0.0,Oligoneuron_rigidum:0.00162)0.927000:0.015239)0.968000:0.028916)0.777000:0.006146,Symphyotrichum_cordifolium:0.020208)1.000000:0.077771)0.843000:0.04181,(Euphorbia_pubentissima:0.125511,(((Quercus_serrata:0.002244,Quercus_ellipsoidalis:0.0068)0.765000:0.036008,Quercus_macrocarpa:0.0)0.906000:0.077946,((Lupinus_perennis:0.058336,Lespedeza_capitata:0.113526)0.306000:0.005766,((Vicia_villosa:0.074364,Astragalus_canadensis:0.045926)1.000000:0.04969,(Amorpha_canadesis:0.023175,(Dalea_purpurea:0.05049,Dalea_candida:0.055157):0.05678)1.000000:0.019707)0.993000:0.02651)0.946000:0.073925)1.000000:0.015029)0.740000:0.037682)0.999000:0.039026)0.988000:0.019496)0.922000:0.0080285);") # Taken from supplementary materials of paper
traits <- read.delim("https://www.cedarcreek.umn.edu/database/datadownload.php?input1=aafe120", skip=51, nrows=38)

# NOTE: chronogram (and so NPRS) is no longer supported in ape, and instead of back-porting it in from v.2.5-1 (which I tried and is slightly tricky because it depends on linked C code) I'm just using chronos (and so the newer Sanderson method) instead
tree <- chronos(tree)

# we can test the crown age with
# class(tree) <- "phylo"
# crown.age(tree)
# and see that it is 1
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

do.all.analyses <- function(mode,drop.zeros=FALSE) {

c.data <- comparative.comm(tree, comm)

if (drop.zeros) {
  rows.to.go <- c()
  for (i in 1:length(c.data$comm[,1])) {
    if (sum(c.data$comm[i,1:8]) ==0) {
      rows.to.go <- c(rows.to.go,i)
    } else {
      if (sum(c.data$comm[i,9:18]) ==0) {
        rows.to.go <- c(rows.to.go,i)
      } 
    }
  }
  
  if (length(rows.to.go) >0) {
    c.data$comm <- c.data$comm[-rows.to.go, ]
  }
}

if (mode == "onlygrass") {
  c.data$comm[,9:18] <- 0 # hack away all the non grasses
}
if (mode == "nograss") {
  c.data$comm[,1:8] <- 0 # hack away all the grasses
}

pd <- setNames(.pd(c.data)[,1], sites(c.data))
pd <- pd[names(biomass)]

gpd.results <- data.frame(matrix("NA", nrow = 0, ncol = 13), stringsAsFactors = FALSE)

df <- data.frame(matrix("NA", nrow = 0, ncol = 1000), stringsAsFactors = FALSE)

names(gpd.results) <- c("rho","pd_stat","gpd_stat","sr_stat","pd_param","gpd_param","sr_param","pd_pval","gpd_pval","sr_pval","pd_est","gpd_est","sr_est")

rho.list <- 10^(seq(from=-4,to=4,by=0.2))

best.rho <- rho.list[0]
best.corr <- 0
best.clean.data <- NULL

for (rho in rho.list) {
    
gpd <- setNames(gpdpez(c.data,rho,origin.life=NULL,std.units=TRUE), sites(c.data))
gpd <- gpd[names(biomass)]

# Match metrics and biomass data, then aggregate by years in Cadotte et al.
clean.data <- data.frame(ids=names(pd))
clean.data$site <- sapply(strsplit(clean.data$ids, "_"), function(x) x[1])
clean.data$year <- sapply(strsplit(clean.data$ids, "_"), function(x) x[2])
clean.data$pd <- pd[clean.data$ids]
clean.data$gpd <- gpd[clean.data$ids]
clean.data$biomass <- biomass[clean.data$ids]
clean.data$n.spp <- rowSums(c.data$comm>0)[clean.data$ids]
clean.data <- clean.data[clean.data$year %in% c(1996:2007),]
clean.data <- aggregate(clean.data, list(clean.data$site), mean) # Ignore warnings; that's the sites/ids being meaned (and they can't be)

clean.data$site <- NULL; clean.data$ids <- NULL; clean.data$year <- NULL; names(clean.data)[1] <- "site"

# Remove sites with no species present
# - this appears to be real - a site where nothnig was sorted and nothing is present as being sown (I assume 'negative controls' to see what is dispersing in)
clean.data <- clean.data[!clean.data$site %in% c(126,166,244,248),]


print(rho)
pd.res <- (cor.test(clean.data$biomass, clean.data$pd))
gpd.res <- (cor.test(clean.data$biomass, clean.data$gpd))
sr.res <- (cor.test(clean.data$biomass, clean.data$n.spp))
gpd.res.row <- c(rho,pd.res$statistic,gpd.res$statistic,sr.res$statistic,pd.res$parameter,gpd.res$parameter,sr.res$parameter,pd.res$p.value,gpd.res$p.value,sr.res$p.value,pd.res$estimate,gpd.res$estimate,sr.res$estimate)
gpd.results <- rbind(gpd.results,gpd.res.row) 

if (gpd.res$estimate > best.corr) {
  best.corr <- gpd.res$estimate
  best.rho <- rho
  best.clean.data <- clean.data
}

}
names(gpd.results) <- c("rho","pd_stat","gpd_stat","sr_stat","pd_param","gpd_param","sr_param","pd_pval","gpd_pval","sr_pval","pd_est","gpd_est","sr_est")

filename <- paste(results.file.path,"EH.PProd.results",sep="")

if (mode == "onlygrass") {
  filename <- paste(results.file.path,"EH.PProd.results.onlygrass",sep="")
}
if (mode == "nograss") {
  filename <- paste(results.file.path,"EH.PProd.results.nograss",sep="")
}

if (drop.zeros) {
  filename <- paste(filename,"drop0",sep=".")
}

filename <- paste(filename,"rda",sep=".")

save(gpd.results,best.rho,best.corr,best.clean.data,file = filename)

}

do.all.analyses("all",FALSE)
do.all.analyses("nograss",FALSE)
do.all.analyses("onlygrass",FALSE)

do.all.analyses("all",TRUE)
do.all.analyses("onlygrass",TRUE)
do.all.analyses("nograss",TRUE)

