
###### Simulations ######

# Simulate data with different diversities and see what effect it has on network complexity

#this came from the RandomizationsFaust.R script which was taken from Faust et al 2015  http://psbweb05.psb.ugent.be/conet/documents/networksFromSimCounts.R

library(dirmult)
library(HMP)
library(vegan)

#Bacteria rarefied to 7987
#lo has 8442 taxa (8426 with no plants), 5854 bacteria, 230 make it into networks
#me has 10971 taxa (10937 with no plants), 7910 bacteria, 241 into networks
#hi has 11942 taxa (11893 with no plants), 7882 bacteria, 209 into networks

#First generate probabilities of each taxon, then simulate counts, then randomize the counts within each species
#the soil ssamples had sheldon's evenness of about .6, and a k=0.08 gives this. theta=0.0009 is from grassland
sim<-t(simMat(taxa=140, samples=25, counts=1000, distrib="dm", maxcount=100, mode=6, k=0.04, theta=0.0009, norm=F, shuffle.samples=F))
sim<-t(simMat(taxa=500, samples=25, counts=1000, distrib="dm", maxcount=100, mode=7, k=10, theta=0.0009, norm=F, shuffle.samples=F))
sort(colSums(sim>0),decreasing = T) #how many plots each taxon was present in
length(which(colSums(sim>0)>11)) # how many taxa would pass my 11 frequency cutoff
sim

sim[1:10,1:10]
plot(sim[,1],sim[,2])

#richness in the different plots (use t())
k=0.08: richness of 50, 54 pass filter 
k=0.04: richness of 97, 93 pass filter
rowSums(sim>0)

#just modifying k with "unlimited" taxa does not work, more richness mean more species pass filters
sim1<-t(simMat(taxa=140, samples=25, counts=1000, distrib="dm", maxcount=100, mode=6, k=0.08, theta=0.02, norm=F, shuffle.samples=F))
mean(rowSums(sim1>0))
length(which(colSums(sim1>0)>11)) # how many taxa would pass my 11 frequency 
mean(vegan::diversity(sim1)/log(specnumber(sim1)))
sort(colSums(sim1>0),decreasing = T) #how many plots each taxon was present in

sim2<-t(simMat(taxa=140, samples=25, counts=1000, distrib="dm", maxcount=100, mode=6, k=0.04, theta=0.02, norm=F, shuffle.samples=F))
mean(rowSums(sim2>0))
length(which(colSums(sim2>0)>11)) # how many taxa would pass my 11 frequency 
mean(vegan::diversity(sim2)/log(specnumber(sim2)))
sort(colSums(sim2>0),decreasing = T) #how many plots each taxon was present in

#just modifying k with "unlimited" taxa does not work, more richness mean more species pass filters
#less rich
sim1<-t(simMat(taxa=100, samples=25, counts=1000, distrib="dm", maxcount=100, mode=6, k=0.08, theta=0.02, norm=F, shuffle.samples=F))
mean(rowSums(sim1>0))
length(which(colSums(sim1>0)>11)) # how many taxa would pass my 11 frequency 
mean(vegan::diversity(sim1)/log(specnumber(sim1)))
sort(colSums(sim1>0),decreasing = T) #how many plots each taxon was present in
#more rich
sim2<-t(simMat(taxa=100, samples=25, counts=1000, distrib="dm", maxcount=100, mode=6, k=0.04, theta=0.02, norm=F, shuffle.samples=F))
mean(rowSums(sim2>0))
length(which(colSums(sim2>0)>11)) # how many taxa would pass my 11 frequency 
mean(vegan::diversity(sim2)/log(specnumber(sim2)))
sort(colSums(sim2>0),decreasing = T) #how many plots each taxon was present in

#the difference between these datasets structure and my actual data is that my rich plots have many more singletons, whereas here the rich plots have fewer singletons. I think I need to hack the distribution and just add some singletons on the end of probabs (and restandardize to 1) #this works a little bit
#33% of my taxa should be singletons

probabs<-getProbabs(taxa=100,mode=6,k=.08)
probabs2<-c(probabs[1:50],rep(probabs[50],50))
probabs3<-probabs2/sum(probabs2)
sim<-t(simCounts(samples=25,counts=1000,pi=probabs,theta=.02,distrib="dm",maxcount = 100))#theta=.0009 was from grassland
mean(rowSums(sim>0))
sort(colSums(sim>0),decreasing = T) #
length(which(colSums(sim>0)>11)) 

sim2<-t(simCounts(samples=25,counts=1000,pi=probabs3,theta=.02,distrib="dm",maxcount = 100))#theta=.0009 was from grassland
mean(rowSums(sim2>0))
sort(colSums(sim2>0),decreasing = T) #
length(which(colSums(sim2>0)>11)) 


#taking distribution from my actual data
#this works...preserves the differences in richness and in how many species pass my cutoff filter
probabslo<-templo
simlo<-t(simCounts(samples=25,counts=1000,pi=probabslo,theta=.02,distrib="dm",maxcount = 100))#theta=.0009 was from grassland
mean(rowSums(simlo>0)) #mean richness
sort(colSums(simlo>0),decreasing = T) #distribution of species frequencies
length(which(colSums(simlo>0)>11))  #passing filter

probabshi<-temphi
simhi<-t(simCounts(samples=25,counts=1000,pi=probabshi,theta=.02,distrib="dm",maxcount = 100))#theta=.0009 was from grassland
mean(rowSums(simhi>0)) #mean richness
sort(colSums(simhi>0),decreasing = T) #distribution of species frequencies
length(which(colSums(simhi>0)>11))  #passing filter







######## Functions ########

#functions

##############################################################################
# Generate a taxon probability vector
# 
# Arguments:
# taxa = the length of the taxon probability vector
# mode = 1-6
#        1 = perfect evenness (pi = rep(1/N,N))
#        2 = probabilities sampled from uniform distribution and normalized to one
#        3 = dominant species has a probability of 0.95 and all other species are equally distributed to add probabilities up to 1
#        4 = probabilities are sampled from a poisson distribution with lambda set to 0.5
#        5 = probabilities are generated using broken-stick function bstick() from vegan
#        6 = probabilities are generated with the geometric series of given evenness k
#             k    = evenness parameter for mode 6
#
# Value:
# A taxon probability vector
#
##############################################################################

probabs<-
sort(getProbabs(taxa=50,mode=7,k=10) )
getProbabs(taxa=50,mode=5) 
probabs<-getProbabs(taxa=50,mode=6,k=.08) #the soil samples had sheldon's evenness of about .6, and a k=0.08 gives this
sort(probabs,decreasing=T)
hist(probabs)
#sheldon(probabs)

probabs<-getProbabs(taxa=50,mode=4) 
hist(probabs)


getProbabs<-function(taxa = 10, mode = 1, k = 0.5){
  if(taxa < 1){
    stop("taxa should be at least 1.")
  }
  #  if(!is.wholenumber(mode)){
  #    stop("The mode should be an integer between 1 and 6")
  #  }
  if(mode == 1){
    pi=rep(1/taxa,taxa)
  }else if(mode == 2){
    pi = runif(taxa)
    pi =pi/sum(pi)
  }else if(mode == 3){
    dominant=0.95
    commoner=0.05/(taxa-1)
    pi=c(dominant,rep(commoner,(taxa-1)))
  }else if(mode == 4){
    pi=rpois(taxa,lambda=0.5)
    pi=pi/sum(pi)
  }else if(mode == 5){
    pi = bstick(taxa,1)
  }else if(mode == 6){
    pi=geom.series(l=taxa, counts=1, k = k)
  }else if(mode == 7){
    pi=rnbinom(n=taxa, mu=.5, size = k)
    pi=pi/sum(pi)
  }else if(mode < 1 || mode > 7){
    stop("Choose a mode between 1 and 6.")
  }
  pi
}

geom.series<-function(l=10, counts=1000, k=0.5){
  if(l < 1){
    stop("l should be at least 1.")
  }
  if(k > 1 || k < 0){
    stop("k should be a number between 0 and 1.")
  }
  pi=c()
  C = 1/(1 - ((1-k)^l))
  for(i in 1:l){
    pi=c(pi,counts*C*k*(1-k)^(i-1))
  }
  pi=pi/sum(pi)
  pi
}

##############################################################################
# Generate counts with Dirichlet-multinomial (dm) or uniform distribution (unif) 
# 
# Arguments:
# samples = number of samples to generate
# counts = total number of counts in each sample
# pi = taxon probability vector, its length indicates the number of taxa
# theta = overdispersion parameter for dm
# maxcount = maximal count number for any taxon (only for uniform distribution)
#
# Value:
# a count vector or matrix
#
# Note: taken partly from package HMP, function Dirichlet.multinomial and 
# package dirmult, function simPop  
#
##############################################################################

sim<-simCounts(samples=5,counts=500,pi=probabs,theta=.02,distrib="dm",maxcount = 100) ;sim#theta=.0009 was from grassland

#sprintf("%.20f",sum(pi)) #some rounding error somewhere so I commented out the three lines below

simCounts<-function(samples=1,counts=1000,pi=rep(1/10,10),theta=0.002, distrib="dm", maxcount=100){
  #  if(sum(pi) != 1){
  #    stop("Taxon probabilities should sum to one!")
  #  }
  if(theta > 1 || theta < 0){
    stop("Theta should be between 0 and 1.")
  }
  if(distrib == "dm"){
    # from a line in function ?simPop in package dirmult
    gamma = pi*(1 - theta)/theta
    # from a line in Dirichlet.multinomial from package HMP
    res=matrix(NA,nrow=length(pi),ncol=samples)
    for(i in 1:samples){#Emily put this in a for loop
      res[,i]=rmultinom(n=1, size=counts, prob=rdirichlet(1,gamma)) #Emily changed this from prob=rdirichlet(samples,gamma) to prob=rdirichlet(1,gamma) b/c the rmultinom function was not using the different rows as separate samples but rather concatenating them (i.e. if you have 50 taxa and 3 samples, the old way would end up with 150 taxa)
    }
  }else if(distrib == "unif"){
    if(maxcount < 1){
      stop("The maximal count should be at least 1.")
    }
    res=matrix(runif(length(pi)*samples,0,maxcount),nrow=length(pi),ncol=samples)
    res = round(res)
  }else{
    stop("Choose dm or unif as generating distribution.")
  }
  res
}

##############################################################################
# Simulate a count matrix
# 
# Arguments:
# taxa = the number of rows in the matrix
# samples = the number of columns in the matrix
# counts = either the total number of counts in each sample or a vector specifying the count number in each sample
# distrib = "dm" for Dirichlet-Multinomial distribution or "unif" for uniform distribution
# maxcount = maximal count number for any taxon (only for uniform distribution)
# mode   = the mode according to which the taxon probability vector for the DM distribution is to be generated, please find the description of the modes for function "getProbabs"
# k      = evenness parameter of mode 6
# theta = overdispersion parameter of the DM distribution
# norm = normalize matrix column-wise, such that the entries in each column add to one
# shuffle.samples = shuffle each sample
#
# Value:
# A count matrix or relative abundance matrix
#
##############################################################################
#sim<-simMat(taxa=50, samples=2, counts=1000, distrib="dm", maxcount=100, mode=6, k=0.08, theta=0.0009, norm=F, shuffle.samples=T)

simMat<-function(taxa=5, samples=10, counts=1000, distrib="dm", maxcount=100, mode=1, k=0.5, theta=0.002, norm=F, shuffle.samples=F){
  if(length(counts)==1 && counts < 1){
    stop("counts should be at least 1.")
  }
  if(taxa < 1){
    stop("There should be at least one taxon.")
  }
  if(samples < 1){
    stop("There should be at least one sample.")
  }
  pi = getProbabs(taxa=taxa, mode=mode, k=k)
  mat=matrix(nrow=taxa,ncol=samples)
  # get taxon probabilities for each column
  sample.count = counts
  for(i in 1:samples){ 
    if(length(counts)==samples){
      sample.count = counts[i]
    }
    mat[,i]=simCounts(pi=pi, theta=theta, counts = sample.count, distrib=distrib, maxcount=maxcount) 
    if(shuffle.samples){
      mat[,i]=sample(mat[,i])
    }
  }
  # normalize matrix column-wise
  if(norm){
    mat = normalize(mat)
  }
  mat
}


