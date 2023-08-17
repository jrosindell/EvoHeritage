# Repeating Cadotte et al. (2009) PLoS One for GPD
# Will Pearse - 2022-06-22

rm(list=ls())

# Headers
library(pez)

source("../EvoHeritage Tools.R")
source("./Primary Productivity Configuration.R") # keep paths and various other details in a separate config file

load(file=paste(inputs.file.path,"PProd.input.data.rda",sep=""))

do.all.analyses <- function(mode,drop.zeros=FALSE) {

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

