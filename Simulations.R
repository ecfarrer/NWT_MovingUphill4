
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


###### Taking distribution from my actual data ######
#taking the distribution from one sample in lo and hi. this works...preserves the differences in richness and in how many species pass my cutoff filter
probabslo<-templo
simlo<-t(simCounts(samples=25,counts=2000,pi=probabslo,theta=.01,distrib="dm",maxcount = 100))#theta=.0009 was from grassland, originally I was using 0.02
mean(rowSums(simlo>0)) #mean richness
#sort(colSums(simlo>0),decreasing = T) #distribution of species frequencies
length(which(colSums(simlo>0)>11))  #passing filter

probabshi<-temphi
simhi<-t(simCounts(samples=25,counts=2000,pi=probabshi,theta=.01,distrib="dm",maxcount = 100))#theta=.0009 was from grassland
mean(rowSums(simhi>0)) #mean richness
#sort(colSums(simhi>0),decreasing = T) #distribution of species frequencies
length(which(colSums(simhi>0)>11))  #passing filter

#taking an average of all samples of lo and hi (this doesn't create the same probablity distribution as just taking one sample but it will allow me to use higher read counts)
probabslo<-lovec
simlo<-t(simCounts(samples=25,counts=8000,pi=probabslo,theta=.0037,distrib="dm",maxcount = 100))#theta=.0009 was from grassland
mean(rowSums(simlo>0)) #mean richness
#sort(colSums(simlo>0),decreasing = T) #distribution of species frequencies
length(which(colSums(simlo>0)>11))  #passing filter

probabshi<-hivec
simhi<-t(simCounts(samples=25,counts=8000,pi=probabshi,theta=.0031,distrib="dm",maxcount = 100))#theta=.0009 was from grassland
mean(rowSums(simhi>0)) #mean richness
#sort(colSums(simhi>0),decreasing = T) #distribution of species frequencies
length(which(colSums(simhi>0)>11))  #passing filter





#### Run loop for short vectors ####
#set probabilities if you want just one
#probabslo<-templo #from plot 1
#probabshi<-temphi

#generate x if you want to
#simx<-cbind(rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1))

#generate output dataframes
outputlo<-data.frame(rep=rep(NA,25),richness=rep(NA,25),passcutoff=rep(NA,25),pos=rep(NA,25),neg=rep(NA,25),tot=rep(NA,25),vertices=rep(NA,25),complexity=rep(NA,25))
outputlomod<-list(NA)
outputlorescor<-list(NA)
outputhi<-data.frame(rep=rep(NA,25),richness=rep(NA,25),passcutoff=rep(NA,25),pos=rep(NA,25),neg=rep(NA,25),tot=rep(NA,25),vertices=rep(NA,25),complexity=rep(NA,25))
outputhimod<-list(NA)
outputhirescor<-list(NA)

#lo
for(i in 1:25){
  temp<-sort(t(hmscYlo2Bac[i,]),decreasing = T)[1:1000]#800 I used 800 for my simulation runs
  probabslo<-temp/sum(temp)  #take this and use it as probab
  
  set.seed(i+4)
  simlo<-t(simCounts(samples=25,counts=2000,pi=probabslo,theta=.01,distrib="dm")) #counts=2000, theta=.01
  outputlo[i,"richness"]<-mean(rowSums(simlo>0)) #mean richness
  #sort(colSums(simlo>0),decreasing = T) #distribution of species frequencies
  outputlo[i,"passcutoff"]<-length(which(colSums(simlo>0)>11))  #passing filter
  
  ind<-which(colSums(simlo)>0);length(ind)
  simlo2<-simlo[,ind]
  #impute zeros
  simlo3 <- cmultRepl(simlo2,label=0, method="CZM")
  simlo3[simlo3<=0]<-min(simlo3[simlo3>0])
  #Calculate clr
  simlo4 <- t(apply(simlo3, 1, function(x){log(x) - mean(log(x))}))
  ind<-which(colSums(simlo2>0)>11)
  simlo5<-simlo4[,ind]
  
  #generate random variables for the environment
  #simx<-cbind(rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1))
  mod.lo<- boral(y = simlo5, X = NULL, lv.control = list(num.lv = 2), family = "normal", save.model = TRUE, calc.ics = F, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
  rescor.lo <- get.residual.cor(mod.lo) 
  
  outputlomod[[i]]<-mod.lo #this saves the data too $X and $y
  outputlorescor[[i]]<-rescor.lo
  
  colMatlo<-rescor.lo$sig.correlaton
  colMatlo[which(rescor.lo$sig.correlaton>0)]<-1
  colMatlo[which(rescor.lo$sig.correlaton<0)]<- -1
  
  graphlo1<-graph_from_adjacency_matrix(colMatlo, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
  myedgelistlo<-data.frame(as_edgelist(graphlo1),weight=E(graphlo1)$weight) #just the edges
  
  graphlo2<-graph.edgelist(as.matrix(myedgelistlo[,1:2]),directed=FALSE)
  length(E(graphlo2))/length(V(graphlo2))

  outputlo[i,"rep"]<-i
  outputlo[i,"pos"]<-length(which(myedgelistlo$weight==1))
  outputlo[i,"neg"]<-length(which(myedgelistlo$weight==-1))
  outputlo[i,"tot"]<-length(which(myedgelistlo$weight==1))+length(which(myedgelistlo$weight==-1))
  outputlo[i,"vertices"]<-length(V(graphlo2))
  outputlo[i,"complexity"]<-length(E(graphlo2))/length(V(graphlo2))
  print(outputlo)
}


#hi
for (j in 1:25){
  temp<-sort(t(hmscYhi2Bac[j,]),decreasing = T)[1:1500]
  probabshi<-temp/sum(temp)  #take this and use it as probab
  
  set.seed(j+4)
  simhi<-t(simCounts(samples=25,counts=2000,pi=probabshi,theta=.01,distrib="dm"))
  outputhi[j,"richness"]<-mean(rowSums(simhi>0)) #mean richness
  #sort(colSums(simhi>0),decreasing = T) #distribution of species frequencies
  outputhi[j,"passcutoff"]<-length(which(colSums(simhi>0)>11))  #passing filter
  
  ind<-which(colSums(simhi)>0);length(ind)
  simhi2<-simhi[,ind]
  #impute zeros
  simhi3 <- cmultRepl(simhi2,label=0, method="CZM")
  simlo3[simlo3<=0]<-min(simlo3[simlo3>0])
  #Calculate clr
  simhi4 <- t(apply(simhi3, 1, function(x){log(x) - mean(log(x))}))
  ind<-which(colSums(simhi2>0)>11)
  simhi5<-simhi4[,ind]
  
  #generate random variables for the environment
  #simx<-cbind(rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1))
  mod.hi<- boral(y = simhi5, X = NULL, lv.control = list(num.lv = 2), family = "normal", save.model = TRUE, calc.ics = F, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
  rescor.hi <- get.residual.cor(mod.hi) 

  outputhimod[[j]]<-mod.hi #this saves the data too $X and $y
  outputhirescor[[j]]<-rescor.hi
  
  colMathi<-rescor.hi$sig.correlaton
  colMathi[which(rescor.hi$sig.correlaton>0)]<-1
  colMathi[which(rescor.hi$sig.correlaton<0)]<- -1
  
  graphhi1<-graph_from_adjacency_matrix(colMathi, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
  myedgelisthi<-data.frame(as_edgelist(graphhi1),weight=E(graphhi1)$weight) #just the edges
  
  graphhi2<-graph.edgelist(as.matrix(myedgelisthi[,1:2]),directed=FALSE)
  length(E(graphhi2))/length(V(graphhi2))
  
  outputhi[j,"rep"]<-j
  outputhi[j,"pos"]<-length(which(myedgelisthi$weight==1))
  outputhi[j,"neg"]<-length(which(myedgelisthi$weight==-1))
  outputhi[j,"tot"]<-length(which(myedgelisthi$weight==1))+length(which(myedgelisthi$weight==-1))
  outputhi[j,"vertices"]<-length(V(graphhi2))
  outputhi[j,"complexity"]<-length(E(graphhi2))/length(V(graphhi2))
  print(outputhi)
}

#FINAL RESULTS, USE THIS FOR MANUSCRIPT
outputloSF<-outputlo
outputhiSF<-outputhi
outputSF<-rbind(outputloSF,outputhiSF)

plot(rbind(outputloSF,outputhiSF)$richness,rbind(outputloSF,outputhiSF)$tot)
abline(lm(rbind(outputloSF,outputhiSF)$tot~rbind(outputloSF,outputhiSF)$richness))
summary(lm(rbind(outputloSF,outputhiSF)$tot~rbind(outputloSF,outputhiSF)$richness))

plot(rbind(outputloSF,outputhiSF)$passcutoff,rbind(outputloSF,outputhiSF)$tot)
abline(lm(rbind(outputloSF,outputhiSF)$tot~rbind(outputloSF,outputhiSF)$passcutoff))
summary(lm(rbind(outputloSF,outputhiSF)$tot~rbind(outputloSF,outputhiSF)$passcutoff))

plot(rbind(outputloSF,outputhiSF)$richness,rbind(outputloSF,outputhiSF)$complexity)
abline(lm(rbind(outputloSF,outputhiSF)$complexity~rbind(outputloSF,outputhiSF)$richness))
summary(lm(rbind(outputloSF,outputhiSF)$complexity~rbind(outputloSF,outputhiSF)$richness))

plot(rbind(outputloSF,outputhiSF)$passcutoff,rbind(outputloSF,outputhiSF)$complexity)
abline(lm(rbind(outputloSF,outputhiSF)$complexity~rbind(outputloSF,outputhiSF)$passcutoff))
summary(lm(rbind(outputloSF,outputhiSF)$complexity~rbind(outputloSF,outputhiSF)$passcutoff))

outputSF2<-outputSF%>%
  filter(complexity<2)
plot(outputSF2$richness,outputSF2$complexity)
abline(lm(outputSF2$complexity~outputSF2$richness))
summary(lm(outputSF2$complexity~outputSF2$richness))

plot(outputSF2$passcutoff,outputSF2$complexity)
abline(lm(outputSF2$complexity~outputSF2$passcutoff))
summary(lm(outputSF2$complexity~outputSF2$passcutoff))



#### Run loop for long vectors ####

#I wrote this to test it out, but never used it really b/c the models take so long to run. If I try the short vectors code above with a theta=0.003, with 2000 reads, there is for example 334 average richness and 312 taxa pass my filter, my actual data with >7000 reads has a richness of 638 and 306 pass filter. Thus the numbers just seemed off, there is too much similarity in rare species comp among samples, the theta was not reproducing the variability in the samples well.

probabslo<-dmodlopi #the full 25 plot abundance vector
probabshi<-dmodhipi

#generate x if you want to
#simx<-cbind(rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1))

#generate output dataframes
outputlo<-data.frame(rep=rep(NA,10),richness=rep(NA,10),passcutoff=rep(NA,10),pos=rep(NA,10),neg=rep(NA,10),tot=rep(NA,10))
outputhi<-data.frame(rep=rep(NA,10),richness=rep(NA,10),passcutoff=rep(NA,10),pos=rep(NA,10),neg=rep(NA,10),tot=rep(NA,10))
#outputlomod<-list(NA)
#outputlorescor<-list(NA)

#lo
for(i in 1:10){
#  temp<-sort(t(hmscYlo2Bac[i,]),decreasing = T)[1:800]#800 I used 800 for my simulation runs
#  probabslo<-temp/sum(temp)  #take this and use it as probab
  
  set.seed(i+4)
  simlo<-t(simCounts(samples=25,counts=7987,pi=probabslo,theta=.00372,distrib="dm",maxcount = 100)) #counts=2000, theta=.01
  outputlo[i,"richness"]<-mean(rowSums(simlo>0)) #mean richness
  #sort(colSums(simlo>0),decreasing = T) #distribution of species frequencies
  outputlo[i,"passcutoff"]<-length(which(colSums(simlo>0)>11))  #passing filter
  
  ind<-which(colSums(simlo)>0);length(ind)
  simlo2<-simlo[,ind]#375 species
  #impute zeros
  simlo3 <- cmultRepl(simlo2,label=0, method="CZM")
  #Calculate clr
  simlo4 <- t(apply(simlo3, 1, function(x){log(x) - mean(log(x))}))
  ind<-which(colSums(simlo2>0)>11)
  simlo5<-simlo4[,ind]

  #generate random variables for the environment
  #simx<-cbind(rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1))
  mod.lo<- boral(y = simlo5, X = NULL, lv.control = list(num.lv = 3), family = "normal", save.model = TRUE, calc.ics = F, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
  rescor.lo <- get.residual.cor(mod.lo) 
  colMatlo<-rescor.lo$sig.correlaton
  colMatlo[which(rescor.lo$sig.correlaton>0)]<-1
  colMatlo[which(rescor.lo$sig.correlaton<0)]<- -1
  
  graphlo1<-graph_from_adjacency_matrix(colMatlo, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
  myedgelistlo<-data.frame(as_edgelist(graphlo1),weight=E(graphlo1)$weight) #just the edges
  
  outputlo[i,"rep"]<-i
  outputlo[i,"pos"]<-length(which(myedgelistlo$weight==1))
  outputlo[i,"neg"]<-length(which(myedgelistlo$weight==-1))
  outputlo[i,"tot"]<-length(which(myedgelistlo$weight==1))+length(which(myedgelistlo$weight==-1))
  print(outputlo)
}

#hi
for (j in 1:10){
#  temp<-sort(t(hmscYhi2Bac[j,]),decreasing = T)[1:1500]
#  probabshi<-temp/sum(temp)  #take this and use it as probab
  
  set.seed(j+4)
  simhi<-t(simCounts(samples=25,counts=7987,pi=probabshi,theta=.00309,distrib="dm",maxcount = 100))
  outputhi[j,"richness"]<-mean(rowSums(simhi>0)) #mean richness
  #sort(colSums(simhi>0),decreasing = T) #distribution of species frequencies
  outputhi[j,"passcutoff"]<-length(which(colSums(simhi>0)>11))  #passing filter
  
  ind<-which(colSums(simhi)>0);length(ind)
  simhi2<-simhi[,ind]#375 species
  #impute zeros
  simhi3 <- cmultRepl(simhi2,label=0, method="CZM")
  #Calculate clr
  simhi4 <- t(apply(simhi3, 1, function(x){log(x) - mean(log(x))}))
  ind<-which(colSums(simhi2>0)>11)
  simhi5<-simhi4[,ind]
  
  #generate random variables for the environment
  #simx<-cbind(rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1),rnorm(25,mean=0,sd=1))
  mod.hi<- boral(y = simhi5, X = NULL, lv.control = list(num.lv = 2), family = "normal", save.model = TRUE, calc.ics = F, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
  rescor.hi <- get.residual.cor(mod.hi) 
  colMathi<-rescor.hi$sig.correlaton
  colMathi[which(rescor.hi$sig.correlaton>0)]<-1
  colMathi[which(rescor.hi$sig.correlaton<0)]<- -1
  
  graphhi1<-graph_from_adjacency_matrix(colMathi, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
  myedgelisthi<-data.frame(as_edgelist(graphhi1),weight=E(graphhi1)$weight) #just the edges
  
  outputhi[j,"rep"]<-j
  outputhi[j,"pos"]<-length(which(myedgelisthi$weight==1))
  outputhi[j,"neg"]<-length(which(myedgelisthi$weight==-1))
  outputhi[j,"tot"]<-length(which(myedgelisthi$weight==1))+length(which(myedgelisthi$weight==-1))
  print(outputhi)
}

mean(outputlo$tot)
std.error(outputlo$tot)
mean(outputhi$tot)
std.error(outputhi$tot)


#1000 counts
outputlo1<-outputlo
outputhi1<-outputhi

#2000 counts
outputlo2<-outputlo
outputhi2<-outputhi

#3000 counts
outputlo3<-outputlo
outputhi3<-outputhi

#2000 counts each run using a different vector of probabilities and theta=0.01 and 10 reps
outputlo4<-outputlo
outputhi4<-outputhi

plot(rbind(outputlo4,outputhi4)$richness,rbind(outputlo4,outputhi4)$tot)
abline(lm(rbind(outputlo4,outputhi4)$tot~rbind(outputlo4,outputhi4)$richness))
summary(lm(rbind(outputlo4,outputhi4)$tot~rbind(outputlo4,outputhi4)$richness))
plot(rbind(outputlo4,outputhi4)$passcutoff,rbind(outputlo4,outputhi4)$tot)
abline(lm(rbind(outputlo4,outputhi4)$tot~rbind(outputlo4,outputhi4)$passcutoff))
summary(lm(rbind(outputlo4,outputhi4)$tot~rbind(outputlo4,outputhi4)$passcutoff))



ggplot(richdata,aes(x=log10(Plant_Dens+1),y=SR))+# as.numeric(fert),color=species
  labs(x="Plant density",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  #geom_smooth(method=lm,se=F,size=.8,color="black",formula = y ~ poly(x, 2)) +
  facet_wrap(~type,scales="free")





###### Old code trying different types of simulations #####

###### Using the default parameters from gamma distribution ####
#just to see what happens, doesn't work (hi has more species passing cutoff). 
#modifying taxa= does not really work b/c when there is a low number of taxa, all taxa pass filter
simlo<-t(simMat(taxa=150, samples=25, counts=1000, distrib="dm", maxcount=100, mode=6, k=0.08, theta=0.02, norm=F, shuffle.samples=F))
length(which(colSums(simlo>0)>11)) 
t(simlo)
t(simlo)[111:150,]
simhi<-t(simMat(taxa=150, samples=25, counts=1000, distrib="dm", maxcount=100, mode=6, k=0.04, theta=0.02, norm=F, shuffle.samples=F))#theta=0.0009


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


#using probabs1 and probabs2 below from the neg bin distribution, jsut messing with k doesn't make higher diveristy have fewer passing cutoff
#low diversity
sim<-t(simCounts(samples=25,counts=1000,pi=probabs1,theta=.02,distrib="dm",maxcount = 100))#theta=.0009 was from grassland
mean(rowSums(sim>0))
sort(colSums(sim>0),decreasing = T) #
length(which(colSums(sim>0)>11)) 

#high diversity
sim2<-t(simCounts(samples=25,counts=1000,pi=probabs2,theta=.02,distrib="dm",maxcount = 100))#theta=.0009 was from grassland
mean(rowSums(sim2>0))
sort(colSums(sim2>0),decreasing = T) #
length(which(colSums(sim2>0)>11)) 










##### Exploring the negative binomial #####

temp<-dnbinom(0:100, mu=2, size=0.02)
plot(0:100, temp)
temp<-dnbinom(0:100, mu=200, size=.02)
points(0:100, temp, col=2)

#numbers from low plots bacteria
temp<-sort(rnbinom(10000,mu=1.36,size=.011),decreasing = T)
hist(temp)
sum(temp>0)
temp2<-temp/sum(temp)
temp2[1:100]

#numbers from hi plots bacteria
temp<-sort(rnbinom(10000,mu=1.01,size=.014),decreasing = T)
hist(temp)
sum(temp>0)
temp2<-temp/sum(temp)
temp2[1:400]


##### Exploring my data structure #####

hmscYlo2Bac
hmscYlo2S

hmscYme2Bac
hmscYme2S

hmscYhi2Bac
hmscYhi2S

dim(hmscYlo2ITS)
colnames(hmscYlo2ITS)<-1:1325
dim(hmscYhi2ITS)
colnames(hmscYhi2ITS)<-1:1941
rarefied to 1023

#richness
mean(rowSums(hmscYlo2Bac>0))
mean(rowSums(hmscYhi2Bac>0))

length(which(colSums(hmscYlo2ITS)==1))
length(which(colSums(hmscYhi2ITS)==1))

#evenness
mean(vegan::diversity(hmscYlo2Bac)/log(specnumber(hmscYlo2Bac)))
mean(vegan::diversity(hmscYme2Bac)/log(specnumber(hmscYme2Bac)))

#more rich plots are slightly more even and have many more singletons

temp<-t(hmscYlo2ITS[2,])
sort(temp,decreasing = T)
cbind(t(hmscYlo2ITS[1,]),t(hmscYlo2ITS[4,]))
rowSums(hmscYlo2S)


#take one sample in lo and hi and calculate relative abundance and use that (moderately worked) OR fit a negbin model to it and use those parameters (didn't work hi plots had more passing cutoff)
#templo
sort(t(hmscYlo2Bac[1,]),decreasing = T)
temp<-sort(t(hmscYlo2Bac[1,]),decreasing = T)[1:800]#800 I used 800 for my simulation runs
hist(temp)
templo<-temp/sum(temp)  #take this and use it as probab
plot(temp)
mean(temp)
var(temp)#variance=mu + mu^2/k; k=(mu^2)/(var-mu)
mean(temp)^2/(var(temp)-mean(temp)) #(k)
rnbinom(n=5000,mu=1.36,size=.011)

#temphi
sort(t(hmscYhi2Bac[1,]),decreasing = T)
temp2<-sort(t(hmscYhi2Bac[1,]),decreasing = T)[1:1500]#1500 I used 1500 for my simulation runs
hist(temp2)
temphi<-temp2/sum(temp2)  #take this and use it as probab
mean(temp2)
var(temp2)#variance=mu + mu^2/k; k=(mu^2)/(var-mu)
mean(temp2)^2/(var(temp2)-mean(temp2)) #(k)
rnbinom(n=5000,mu=1.01,size=.014)

#take sum of a few or all samples in lo and hi and use relabun of that as probability vector
lovec<-sort(t(colSums(hmscYlo2Bac)),decreasing=T)
lovec<-lovec/sum(lovec)
hivec<-sort(t(colSums(hmscYhi2Bac)),decreasing=T)
hivec<-hivec/sum(hivec)
cbind(lovec[1:300],hivec[1:300])

#Using 10 samples from lo and 10 from hi to use in each randomization so that richness is different




####### Fitting the dirichlet multinomial to my data to determine theta parameter #######
#looks like it works, it took ~45 min to fit.
#dmodlo<-dirmult(hmscYlo2Bac,trace=T) #rows are subpopultions, columns are categories
dmodlo
#dirmult.summary(hmscYlo2Bac,dmodlo)#you should be able to run this but my matrix is too big and it crashes
dmodlopi<-sort(dmodlo$pi,decreasing = T) #a vector of the probabilities for each taxon
dmodlo$theta # theta=0.00372
#dmodhi<-dirmult(hmscYhi2Bac,trace=T) #
dmodhi
#dirmult.summary(hmscYhi2Bac,dmodhi)
dmodhipi<-sort(dmodhi$pi,decreasing = T) #a vector of the probabilities for each taxon
dmodhi$theta # theta=0.00309

#useful files for looking at my data structures
rm(list=setdiff(ls(), c("labelfile",
                        "richEukS2",
                        "richEukN2",
                        "richBac2",
                        "richITS2",
                        "datEukS5otu",
                        "datEukS5cotu",
                        "datEukN5otu",
                        "datEukN5cotu",
                        "datBacS5otu",
                        "datBacS5cotu",
                        "datITSS5otu",
                        "datITSS5cotu",
                        "datEukS5otu3",
                        "datEukS5cotu3",
                        "datEukN5otu3",
                        "datEukN5cotu3",
                        "datBacS5otu3",
                        "datBacS5cotu3",
                        "datITSS5otu3",
                        "datITSS5cotu3",
                        "datEukS5k2",
                        "datEukN5k2",
                        "datBacS5k2",
                        "datITSS5k2",
                        "plantcomp",
                        "plantcomp2",
                        "comm.bio",
                        "biogeo6")))






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
#            k = evenness parameter for mode 6
#
# Value:
# A taxon probability vector
#
##############################################################################

probabs1<-getProbabs(taxa=400,mode=7,k=.011);probabs1
probabs2<-getProbabs(taxa=400,mode=7,k=.021);probabs2
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
    pi=sort(rnbinom(n=taxa*20, mu=1.36, size = k),decreasing = T)[1:taxa]#sort(rnbinom(10000,mu=1.36,size=.011),decreasing = T)
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


