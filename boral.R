## Fitting a Latent Variable Model using the boral package ##
## Use boral to fit a LVM with 4 latent variables, a random site effect, and the effect of fixed environmental variables

## NOTE: boral uses JAGS so you need to first install JAGS-3.0.0.exe or higher from http://sourceforge.net/projects/mcmc-jags/files/JAGS/

#Workflow: 
#Rarefy each 16S/ITS/18SS/18SN data set (I need to because I'm dealing with differences in richness that are due to environment and I don't want read depth to affect richness, however Gloor does not rarefy data and DESeq won't work b/c it needs multiple samples within a treatment and tells you if taxa are more abundant in one vs the other treatment)
#Merge 16S/ITS/18SS/18SN/plants/biogeo just to get everything together
#Split into lo me hi
#Remove taxa that are zeros
#Split into 16S/ITS/18SS/18SN
#Do the cmultRepl sampling on each 16S/ITS/18SS/18SN dataset which imputes zeros. I went back and forth between first splitting into lo/me/hi vs. doing cmult and clr on whole 16S dataset for example Notes: Originally I was thinking of first splitting each into lo/me/hi (previous notes: not sure if I need to do this, it seems in some ways that I should b/c some taxa are not going to be present not due to incomplete sampling but b/c the environment is not right, but in some ways it probably doesnot matter b/c every sample will be treated the same and b/c a 1 is so similar to a 0 for linear type models). I might be able to do it on the whole datasets b/c I will need the full clr-ed datasets for ordination, so it makes the most sense to do clr on the full datasets. In the end I decided to do it on split lo/me/hi datasets b/c i feel like adding 1's will dilute the dataset and decrease differences between taxa that are actually present (you could be addign 1000 1's to the low dataset for example to acount for all the taxa that are present in me/hi plots but not low)
#Calculate clr on each data set 16S/ITS/18SS/18SN
#Replace back into lo/me/hi datasets
#Delete infrequent taxa
#Do modeling


####### Get species and environment data together #####

#microbes rarefied count data, not filtered
comm.dataEukS<-datEukS5cotu
comm.dataEukN<-datEukN5cotu
comm.dataBac<-datBacS5cotu
comm.dataITS<-datITSS5cotu

max(datEukS5otu[,34:dim(datEukS5otu)[2]])
[1] 0.4006889
max(datEukN5otu[,34:dim(datEukN5otu)[2]])
[1] 1
max(datBacS5otu[,34:dim(datBacS5otu)[2]])
[1] 0.1354701
max(datITSS5otu[,34:dim(datITSS5otu)[2]])
[1] 0.5356794

sort(matrix(as.matrix(datITSS5otu[,34:dim(datITSS5otu)[2]]),ncol=1),decreasing=T)

#plants
plantcomp2


#Merge things. all microbe datasets (not plants) should have the same 90 samples
#first merge comm.dataEuk with comm.data16S
#I need to remove all the description columns in one of the files, then merge
comm.dataEukSa<-cbind(Sample_name=comm.dataEukS$Sample_name,comm.dataEukS[,-c(1:33)])
comm.dataALL1<-merge(comm.dataEukN,comm.dataEukSa,"Sample_name",sort=F,all.y=F,all.x=F)

comm.dataBaca<-cbind(Sample_name=comm.dataBac$Sample_name,comm.dataBac[,-c(1:33)])
comm.dataALL2<-merge(comm.dataALL1,comm.dataBaca,"Sample_name",sort=F,all.y=F,all.x=F)

comm.dataITSa<-cbind(Sample_name=comm.dataITS$Sample_name,comm.dataITS[,-c(1:33)])
comm.dataALL3<-merge(comm.dataALL2,comm.dataITSa,"Sample_name",sort=F,all.y=F,all.x=F)

comm.dataALL3$Sample_name
dim(comm.dataALL3)[2]-33
3161+450+4234+16612 #matches, good, 24457 microbial taxa total

#then merge plants with microbes
comm.dataALL4<-merge(comm.dataALL3,plantcomp2,"Sample_name",sort=F,all.y=F)
comm.dataALL4$Sample_name

#substitute S for the N, since the mapping file was from the nematode dataset
comm.dataALL4$X.SampleID<-sub("N", "S", comm.dataALL4$X.SampleID) 

#delete mapping file data except for X.SampleID
comm.dataALL5<-comm.dataALL4[,-c(1,3:33)]
comm.dataALL5[1:10,1:10]



# biogeochemistry and plant density/cover data
biogeo6$X.SampleID

rcorr(as.matrix(data.frame(biogeo6$TC,biogeo6$snowdepth,biogeo6$pH,biogeo6$moisture)))

#Merge the biogeo6 with comm.dataALL, then split them to make sure the same samples and order are in each dataset
comm.bio<-merge(biogeo6,comm.dataALL5)
comm.bio[1:10,1:60]

#the comm.bio.csv file that is saved here is from the old bioinformatics. I saved it b/c the R environment file takes so long to load. however I might want to reinstate this, if this environment gets really big
#write.csv(comm.bio,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/comm.bio.csv",row.names=F)
#comm.bio<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/comm.bio.csv")





##### Split into lo/me/hi, calculate CLR #####

## lo ##
dim(comm.bio)

ind<-which(comm.bio$lomehi=="lo")
  
hmscYlo<-comm.bio[ind,55:24565]  #start at 54 if you want lomehi column
rownames(hmscYlo)<-comm.bio$X.SampleID[ind]
hmscYlo[1:10,1:10]

hmscXlo<-data.frame(snowdepth=comm.bio$snowdepth,TC=comm.bio$TC,pH=comm.bio$pH,moisture=comm.bio$moisture)[ind,] #,lomehi=comm.bio$lomehi,plantcov=comm.bio$plantcov,whc=comm.bio$WHC
rownames(hmscXlo)<-comm.bio$X.SampleID[ind]

rcorr(as.matrix(hmscXlo))
#plot(hmscX$TC,hmscX$moisture)

#take out species that are zeros
ind<-which(colSums(hmscYlo)>0);length(ind)
hmscYlo2<-hmscYlo[ind]

temp<-lapply(colnames(hmscYlo2),function(x){substr(x,1,1)})
head(which(temp=="I"))
tail(which(temp=="I"))

#separate N, S, 16S, ITS
hmscYlo2N<-hmscYlo2[,1:101]
hmscYlo2S<-hmscYlo2[,102:1247]
hmscYlo2Bac<-hmscYlo2[,1248:7101]
hmscYlo2ITS<-hmscYlo2[,7102:8426]
hmscYlo2Plant<-hmscYlo2[,8427:8442]

#Impute zeros
hmscYlo2N2 <- cmultRepl(hmscYlo2N,label=0, method="CZM")
hmscYlo2S2 <- cmultRepl(hmscYlo2S,label=0, method="CZM")
hmscYlo2Bac2 <- cmultRepl(hmscYlo2Bac,label=0, method="CZM")
hmscYlo2ITS2 <- cmultRepl(hmscYlo2ITS,label=0, method="CZM")

#Calculate clr
hmscYlo2N3 <- t(apply(hmscYlo2N2, 1, function(x){log(x) - mean(log(x))}))
hmscYlo2S3 <- t(apply(hmscYlo2S2, 1, function(x){log(x) - mean(log(x))}))
hmscYlo2Bac3 <- t(apply(hmscYlo2Bac2, 1, function(x){log(x) - mean(log(x))}))
hmscYlo2ITS3 <- t(apply(hmscYlo2ITS2, 1, function(x){log(x) - mean(log(x))}))

#I could rescale plant data essentially by dividing each cell by the max abundance, so that it is preserving all properties (e.g. I dont want to turn it into relative abundance), or by dividing each plant species by its max abundance

hmscYlo3<-cbind(hmscYlo2N3,hmscYlo2S3,hmscYlo2Bac3,hmscYlo2ITS3,hmscYlo2Plant)
hmscYlo2[1:5,1:5]
hmscYlo3[1:5,1:5]

hmscYlo2S2[1:5,1:5]
hmscYlo2S3[1:5,1:5]
colSums(hmscYlo2S2)
rowSums(hmscYlo2S2)
hist(hmscYlo2S2[,25])
hist(hmscYlo2S3[,25])
mean(apply(hmscYlo2S3,2,min))


## medium ##
ind<-which(comm.bio$lomehi=="me")

hmscYme<-comm.bio[ind,55:24565]  #start at 54 if you want lomehi column
rownames(hmscYme)<-comm.bio$X.SampleID[ind]
hmscYme[1:10,1:10]

hmscXme<-data.frame(snowdepth=comm.bio$snowdepth,TC=comm.bio$TC,pH=comm.bio$pH,moisture=comm.bio$moisture)[ind,] #,lomehi=comm.bio$lomehi,plantcov=comm.bio$plantcov,whc=comm.bio$WHC
rownames(hmscXme)<-comm.bio$X.SampleID[ind]

rcorr(as.matrix(hmscXme[,2:(dim(hmscXme)[2]-1)]))

#take out species that are zeros
ind<-which(colSums(hmscYme)>0);length(ind)
hmscYme2<-hmscYme[ind]

temp<-lapply(colnames(hmscYme2),function(x){substr(x,1,1)})
head(which(temp=="I"))
tail(which(temp=="I"))

#separate N, S, 16S, ITS
hmscYme2N<-hmscYme2[,1:167]
hmscYme2S<-hmscYme2[,168:1465]
hmscYme2Bac<-hmscYme2[,1466:9375]
hmscYme2ITS<-hmscYme2[,9376:10937]
hmscYme2Plant<-hmscYme2[,10938:10971]

#Impute zeros
hmscYme2N2 <- cmultRepl(hmscYme2N,label=0, method="CZM")
hmscYme2S2 <- cmultRepl(hmscYme2S,label=0, method="CZM")
hmscYme2Bac2 <- cmultRepl(hmscYme2Bac,label=0, method="CZM")
hmscYme2ITS2 <- cmultRepl(hmscYme2ITS,label=0, method="CZM")

#Calculate clr
hmscYme2N3 <- t(apply(hmscYme2N2, 1, function(x){log(x) - mean(log(x))}))
hmscYme2S3 <- t(apply(hmscYme2S2, 1, function(x){log(x) - mean(log(x))}))
hmscYme2Bac3 <- t(apply(hmscYme2Bac2, 1, function(x){log(x) - mean(log(x))}))
hmscYme2ITS3 <- t(apply(hmscYme2ITS2, 1, function(x){log(x) - mean(log(x))}))

#I could rescale plant data essentially by dividing each cell by the max abundance, so that it is preserving all properties (e.g. I dont want to turn it into relative abundance), or by dividing each plant species by its max abundance

hmscYme3<-cbind(hmscYme2N3,hmscYme2S3,hmscYme2Bac3,hmscYme2ITS3,hmscYme2Plant)
hmscYme2[1:5,1:5]
hmscYme3[1:5,1:5]


## high ##
ind<-which(comm.bio$lomehi=="hi")

hmscYhi<-comm.bio[ind,55:24565]  #start at 54 if you want lomehi column
rownames(hmscYhi)<-comm.bio$X.SampleID[ind]
hmscYhi[1:10,1:10]

hmscXhi<-data.frame(snowdepth=comm.bio$snowdepth,TC=comm.bio$TC,pH=comm.bio$pH,moisture=comm.bio$moisture)[ind,] #,lomehi=comm.bio$lomehi,plantcov=comm.bio$plantcov,whc=comm.bio$WHC
rownames(hmscXhi)<-comm.bio$X.SampleID[ind]

rcorr(as.matrix(hmscXhi[,2:(dim(hmscXhi)[2]-1)]))

#take out species that are zeros
ind<-which(colSums(hmscYhi)>0);length(ind)
hmscYhi2<-hmscYhi[ind]

temp<-lapply(colnames(hmscYhi2),function(x){substr(x,1,1)})
head(which(temp=="I"))
tail(which(temp=="I"))

#separate N, S, 16S, ITS
hmscYhi2N<-hmscYhi2[,1:282]
hmscYhi2S<-hmscYhi2[,283:2070]
hmscYhi2Bac<-hmscYhi2[,2071:9952]
hmscYhi2ITS<-hmscYhi2[,9953:11893]
hmscYhi2Plant<-hmscYhi2[,11894:11942]

#Impute zeros
hmscYhi2N2 <- cmultRepl(hmscYhi2N,label=0, method="CZM")
hmscYhi2S2 <- cmultRepl(hmscYhi2S,label=0, method="CZM")
hmscYhi2Bac2 <- cmultRepl(hmscYhi2Bac,label=0, method="CZM")
hmscYhi2ITS2 <- cmultRepl(hmscYhi2ITS,label=0, method="CZM")

#Calculate clr
hmscYhi2N3 <- t(apply(hmscYhi2N2, 1, function(x){log(x) - mean(log(x))}))
hmscYhi2S3 <- t(apply(hmscYhi2S2, 1, function(x){log(x) - mean(log(x))}))
hmscYhi2Bac3 <- t(apply(hmscYhi2Bac2, 1, function(x){log(x) - mean(log(x))}))
hmscYhi2ITS3 <- t(apply(hmscYhi2ITS2, 1, function(x){log(x) - mean(log(x))}))

#I could rescale plant data essentially by dividing each cell by the max abundance, so that it is preserving all properties (e.g. I dont want to turn it into relative abundance), or by dividing each plant species by its max abundance, but I will leave it as count data and model as neg.bin

hmscYhi3<-cbind(hmscYhi2N3,hmscYhi2S3,hmscYhi2Bac3,hmscYhi2ITS3,hmscYhi2Plant)
hmscYhi2[1:5,1:5]
hmscYhi3[1:5,1:5]



##### Select common species #####
#select species with greater than X (X+1 or more) occurrences

## low ##
ind<-which(colSums(hmscYlo2>0)>11)
#hmscYlo2sub<-hmscYlo2[,ind] #subset the raw data as well just for looking at histograms
length(ind)
hmscYlo4<-hmscYlo3[,ind]
hmscXlo
dim(hmscYlo4)
dim(hmscXlo)

## medium ##
ind<-which(colSums(hmscYme2>0)>11)
length(ind)
hmscYme4<-hmscYme3[,ind]
hmscXme
dim(hmscYme4)
dim(hmscXme)

## high ##
ind<-which(colSums(hmscYhi2>0)>11)
length(ind)
hmscYhi4<-hmscYhi3[,ind]
hmscXhi
dim(hmscYhi4)
dim(hmscXhi)

min(hmscYlo4)
min(hmscYme4)
min(hmscYhi4)
max(hmscYlo4[,1:519])#not including plants for >8
max(hmscYme4[,1:498])
max(hmscYhi4[,1:469])
hmscXlo
hmscXme
hmscXhi

#The y data are not normal (I have lots of options for distributions: normal, binary, poisson, overdispersed poisson (neg bin), log normal, gamma). however, the gloor paper and some other papers on clr seem to just plug clr values into ordinations etc without any additional transformations so I will not transform them here (yet)
hist(hmscYlo2sub[,320],breaks=20)
hist(hmscYlo4[,320],breaks=20)
hist(log(hmscYlo4[,320]+1),breaks=20)
qqnorm(log(hmscYlo4[,117]+1))
qqnorm(hmscYlo4[,117])
plot(hmscYlo2sub[,117],hmscYlo4[,117])

#check if the values are too low that some tolerance is messing up the CI estimates, yes important to scale y, but it is already fairly scaled here. I will scale x, since they differ so much in range
hmscXlo2<-scale(hmscXlo)
hmscXme2<-scale(hmscXme)
hmscXhi2<-scale(hmscXhi)

#Make them matrices
hmscYlo5<-as.matrix(hmscYlo4)
hmscYme5<-as.matrix(hmscYme4)
hmscYhi5<-as.matrix(hmscYhi4)

hmscXlo3<-as.matrix(hmscXlo2)
hmscXme3<-as.matrix(hmscXme2)
hmscXhi3<-as.matrix(hmscXhi2)


##### Fit the LVM using boral and calculate residual correlation matrix#####

#List of files produced:
mod.lo - 4 latent, fewer iterations to check things so it fits fast, >8 frequency
mod.hi
mod.lof - more iterations (f for final model)
mod.mef
mod.hif
mod.lo9 - >9 and neg.bin for plants
mod.me9 - >9 and neg.bin for plants
mod.hi9 - >9 and neg.bin for plants
mod.hi12 - >12 and neg.bin for plants
mod.hi16 - >16 and neg.bin for plants
mod.hi9lv0 - >9 and no latent, neg.bin for plants
mod.hi9fixed0 - >9 and no fixed, neg.bin for plants
mod.lo9fixed0 - >9 and no fixed, neg.bin for plants
mod.hi9pois - >9 and pois for plants
mod.mef9lv2 - testing how many latent variables to use
mod.mef9lv3
mod.mef9lv4
mod.mef9lv5
mod.mef9lv6
mod.lo9f - final model with >9 frequency, long chains, 4 latent variables, plants modeled as neg.bin
mod.me9f
mod.hi9f
mod.lo10f - final model with >10 frequency, long chains, 4 latent variables, plants modeled as neg.bin
mod.me10f
mod.hi10f
mod.lo11f - final model with >11 frequency, long chains, 4 latent variables, plants modeled as neg.bin
mod.me11f
mod.hi11f
mod.lo8f - final model with >8 frequency, long chains, 4 latent variables, plants modeled as neg.bin
mod.me8f
mod.hi8f
mod.hi10fs - s for swapping the order of deleting species then calculating clr, >10, 4 latent, plants neg.bin
mod.lo12f - final model with >12 frequency, long chains, 4 latent variables, plants modeled as neg.bin
mod.me12f
mod.hi12f
*mod.lo11flv3 - final model with >11 frequency, long chains, 3 latent variables, plants modeled as neg.bin
*mod.me11flv3
*mod.hi11flv3
mod.lo12flv3 - final model with >12 frequency, long chains, 3 latent variables, plants modeled as neg.bin

### *THE FINAL MODELS THAT I WILL USE FOR THE MANUSCRIPT: mod.lo11flv3, mod.me11flv3, mod.hi11flv3 (calculating clr then subsetting data, keeping species with >11 frequency, 3 latent variables, plants modeled as neg.big)

#Notes and decisions: 
#changing from >8 to >9 frequency has very little impact on networks, just a few fewer taxa and interactions. 
#changing from normal to neg.bin for plants should be done b/c of nonnormality of plant count data. neg.bin fits better than poisson - residuals with poisson are very large >6, with neg.bin ~3
#adding environmetnal data is super important, NOT for explaining variation per se or explaining co-variation in the models, but rather in defining what is a significant interactoin. There are much much fewer significant interactions after accounting for enviromental variables
#Changing the order of operations to remove >10, then calculate clr, does not change the nubmer of sig interactions but it dramatically changes the direction with many more negative interactions (50% vs. 80% positive interactions)
#using shorter chains for testing is ok to a certain extent but don't expect all the interactions to be included in the final model if you rerun it with longer chains, some shifting will occur

#Fitting notes
#Using the default mcmc parameters, the models take about 30 min to fit.  mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123)
#Using shorter chains, it takes about 12 min to fit.  mcmc.control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 6, seed = 123); I changed nthin from 3 to 6 so that the get.residual.cor would speed up.
#In the tutorial they add this to the model fitting code, hypparams = c(20,20,20,20), however, now if you wanted to change this you need to put it in a prior.control statement or something. I am just using the default here
#Be sure to remove row.eff = "random", a row effect standardizes across samples, but I don't want to do that b/c I already relativized my data and don't want to standardize across microbes and plants


#Using shorter chains, start 3:22pm, end 3:29 for model, 11 min for residual corr matrix
mod.lo9fixed0<- boral(y = hmscYlo5, X = NULL, lv.control = list(num.lv = 4), family = c(rep("normal",437),rep("negative.binomial",1)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 6, seed = 123))#
rescor.lo9fixed0 <- get.residual.cor(mod.lo9fixed0)

mod.me9<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 4), family = c(rep("normal",414),rep("negative.binomial",4)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 6, seed = 123))#
rescor.me9 <- get.residual.cor(mod.me9)

mod.hi9pois <- boral(y = hmscYhi5, X = hmscXhi3, lv.control = list(num.lv = 4), family = c(rep("normal",384),rep("poisson",12)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 6, seed = 123))#
rescor.hi9pois <- get.residual.cor(mod.hi9pois)


#Using longer chains - final models, start=5pm, end=5:30
mod.lo9f<- boral(y = hmscYlo5, X = hmscXlo3, lv.control = list(num.lv = 4), family = c(rep("normal",437),rep("negative.binomial",1)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.lo9f <- get.residual.cor(mod.lo9f) 

mod.me9f<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 4), family = c(rep("normal",414),rep("negative.binomial",4)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.me9f <- get.residual.cor(mod.me9f) 

mod.hi9f<- boral(y = hmscYhi5, X = hmscXhi3, lv.control = list(num.lv = 4), family = c(rep("normal",384),rep("negative.binomial",12)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.hi9f <- get.residual.cor(mod.hi9f) 

mod.lo10f<- boral(y = hmscYlo5, X = hmscXlo3, lv.control = list(num.lv = 4), family = c(rep("normal",357),rep("negative.binomial",1)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.lo10f <- get.residual.cor(mod.lo10f) 

mod.me10f<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 4), family = c(rep("normal",338),rep("negative.binomial",3)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.me10f <- get.residual.cor(mod.me10f) 

mod.hi10f<- boral(y = hmscYhi5, X = hmscXhi3, lv.control = list(num.lv = 4), family = c(rep("normal",311),rep("negative.binomial",10)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.hi10f <- get.residual.cor(mod.hi10f) 


mod.lo11f<- boral(y = hmscYlo5, X = hmscXlo3, lv.control = list(num.lv = 4), family = c(rep("normal",305),rep("negative.binomial",1)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.lo11f <- get.residual.cor(mod.lo11f) 

mod.me11f<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 4), family = c(rep("normal",298),rep("negative.binomial",3)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.me11f <- get.residual.cor(mod.me11f) 

mod.hi11f<- boral(y = hmscYhi5, X = hmscXhi3, lv.control = list(num.lv = 4), family = c(rep("normal",265),rep("negative.binomial",8)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.hi11f <- get.residual.cor(mod.hi11f) 


mod.lo8f<- boral(y = hmscYlo5, X = hmscXlo3, lv.control = list(num.lv = 4), family = c(rep("normal",519),rep("negative.binomial",1)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.lo8f <- get.residual.cor(mod.lo8f) 

mod.me8f<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 4), family = c(rep("normal",498),rep("negative.binomial",4)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.me8f <- get.residual.cor(mod.me8f) 

mod.hi8f<- boral(y = hmscYhi5, X = hmscXhi3, lv.control = list(num.lv = 4), family = c(rep("normal",469),rep("negative.binomial",14)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.hi8f <- get.residual.cor(mod.hi8f) 


#mod.lo12f<- boral(y = hmscYlo5, X = hmscXlo3, lv.control = list(num.lv = 4), family = c(rep("normal",305),rep("negative.binomial",1)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.lo12f <- get.residual.cor(mod.lo12f) 

#mod.me12f<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 4), family = c(rep("normal",298),rep("negative.binomial",3)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.me12f <- get.residual.cor(mod.me12f) 

mod.hi12f<- boral(y = hmscYhi5, X = hmscXhi3, lv.control = list(num.lv = 4), family = c(rep("normal",214),rep("negative.binomial",4)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.hi12f <- get.residual.cor(mod.hi12f) 


mod.lo11flv3<- boral(y = hmscYlo5, X = hmscXlo3, lv.control = list(num.lv = 3), family = c(rep("normal",305),rep("negative.binomial",1)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.lo11flv3 <- get.residual.cor(mod.lo11flv3) 

mod.me11flv3<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 3), family = c(rep("normal",298),rep("negative.binomial",3)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.me11flv3 <- get.residual.cor(mod.me11flv3) 

mod.hi11flv3<- boral(y = hmscYhi5, X = hmscXhi3, lv.control = list(num.lv = 3), family = c(rep("normal",265),rep("negative.binomial",8)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.hi11flv3 <- get.residual.cor(mod.hi11flv3) 


#mod.lo12flv3<- boral(y = hmscYlo5, X = hmscXlo3, lv.control = list(num.lv = 3), family = c(rep("normal",246),rep("negative.binomial",1)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.lo12flv3 <- get.residual.cor(mod.lo12flv3) 

mod.me12flv3<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 3), family = c(rep("normal",237),rep("negative.binomial",3)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.me12flv3 <- get.residual.cor(mod.me12flv3) 

#mod.hi12flv3<- boral(y = hmscYhi5, X = hmscXhi3, lv.control = list(num.lv = 3), family = c(rep("normal",214),rep("negative.binomial",4)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.hi12flv3 <- get.residual.cor(mod.hi12flv3) 



#Testing diffrent numbers of latent variables
mod.mef9lv2<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 2), family = c(rep("normal",414),rep("negative.binomial",4)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.mef9lv2 <- get.residual.cor(mod.mef9lv2) 
mod.mef9lv3<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 3), family = c(rep("normal",414),rep("negative.binomial",4)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.mef9lv3 <- get.residual.cor(mod.mef9lv3) 
mod.mef9lv4<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 4), family = c(rep("normal",414),rep("negative.binomial",4)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.mef9lv4 <- get.residual.cor(mod.mef9lv4) 
mod.mef9lv5<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 5), family = c(rep("normal",414),rep("negative.binomial",4)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.mef9lv5 <- get.residual.cor(mod.mef9lv5) 
mod.mef9lv6<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 6), family = c(rep("normal",414),rep("negative.binomial",4)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.mef9lv6 <- get.residual.cor(mod.mef9lv6) 

mod.mef9lv2$ics[1]
mod.mef9lv3$ics[1]
mod.mef9lv4$ics[1]
mod.mef9lv5$ics[1]
mod.mef9lv6$ics[1]
i<-1
plot(2:6,c(mod.mef9lv2$ics[i],
           mod.mef9lv3$ics[i],
           mod.mef9lv4$ics[i],
           mod.mef9lv5$ics[i],
           mod.mef9lv6$ics[i]),type = "b")

dic2<-get.dic(mod.mef9lv2)
#notes: the more latent variables you add, the fewer significant interactions you have. ex: with 2 latent variables there are 233 taxa with 2440 interactions, and with 4 latent variables there are 214 taxa and 5618 interactions. I think this is because with 4 latent variables you are giving the model more "space" or "dimentions" for the taxa to sort out, so there are fewer significant correlations among taxa. Whe I did this in MovingUphill3, I'm pretty sure the model with 4 latent varialbe had lowest DIC which is why I chose it. Here the model with 2 latent variables has lowest DIC. DIC is suspect though in these models and boral is actually not updating/supporting these functions anymore. I think I will stay with 4 latent variables because I have so many species and I think using only 2 latent variables actually gives spurious significant interactions



#Use modified function below to get 90% CI
#res.corshiocc8.90<-get.residual.cor2(fit.lvmhiocc8)



##### Look at results and check convergence/fit #####

#Model fit, information criteria
mod.lo11flv3$ics
#from when I tested the effect of different numbers of latent variables

summary(mod.lo11flv3) # To look at estimated parameter values
mod.lo$hpdintervals # 95% credible intervals for model parameters.

#Check convergence
#Geweke diagnostic - a z test testing whether the first 10% and the last 50% are diffrent (i think those are the fractions, doesn't really matter exactly), if it is significant, then the means are different and it didn't converge
plot(get.mcmcsamples(mod.lo11flv3)[,1])
plot(get.mcmcsamples(mod.lo11flv3)[,2])

#the order is effect of snowdepth for each of the 600 species, then effect of TC, then pH, then moisture 
mcmchi<-get.mcmcsamples(mod.lo11flv3)
dim(mcmchi)
colnames(mcmchi)[2500:3000]
mcmchi[1:10,1:5]

#TRUE means these did not converge
gew.pvals <- 2*pnorm(abs(unlist(mod.lo11flv3$geweke.diag[[1]])), lower.tail = FALSE)
length(gew.pvals)
gew.pvals[1:5]
gew.pvals[which(gew.pvals<.05)] #technically these did not converge, however, the trace plots look fine to me
p.adjust(gew.pvals, method = "holm")

mod.lo11flv3$geweke.diag
mod.me11flv3$geweke.diag
mod.hi11flv3$geweke.diag
mod.lo11flv3$geweke.diag$prop.exceed
mod.me11flv3$geweke.diag$prop.exceed
mod.hi11flv3$geweke.diag$prop.exceed

#example of one that did not converge
#(1st species) N6f914ead2160e51670d3dc70c25e107b for snowdepth did not converge, but looking at the trace plot, it seems fine
#geweke diagnostic
mod.lo11flv3$geweke.diag$geweke.diag$lv.coefs[1:5,]
#trace plot (it is the very first parameter)
plot(get.mcmcsamples(mod.lo11flv3)[,3])
#mean of the mcmc chain to make sure I'm looking at the right parameter
mean(get.mcmcsamples(fit.hilv4occ9exp4f)[,1]) #mean is -1.710607
#mean of the extracted model coefficients (to make sure I'm looking at the right parameter)
fit.hilv4occ9exp4f$X.coefs.mean  #-1.710607072, yes checks



##### Percent covariation explained by env #####

#One approach to quantify how much of the species co-occurrence is explained by covariates (that is, how well the predictor variables describe the species assemblage) is through differences in the trace of the estimated residual covariance matrix induced by the latent variables (Warton et al. 2015). From the above code, this can be obtained as rescors$trace. For the spider data set, when we compared a pure latent variable model (similar to the equation 1 but without site effects) to the correlated response model, the trace decreased from 178.92 to 107.92. This implies that environmental covariates accounted for approximately 40% of the covariation between species.
#The trace is interpreted as the amount of covariation explained by the latent variables. it is a random variable so will be different from run to run based on mcmc
(178.92-107.92)/178.92

#Need to fit a model with only latent variables
rescor.hi9fixed0$trace #expected to be higher than the model below (but it's not)
rescor.hi9$trace
(rescor.hi9$trace-rescor.hi9fixed0$trace)/rescor.hi9$trace*100
#this didn't work, not sure why


#### Percent variation explained by environment ####
#I don't konw how useful this is (b/c it is not R2 just partitioning the explained variation) - the results suggest that the vast majority of the variance is explained by the latent variables compared to the environment.  
varparthi<-calc.varpart(mod.mef9lv4)#,groupX=c(1,2,3,4,5)
varparthi$varpart.X[1:10]
varparthi$varpart.lv[1:10]
hist(varparthi$varpart.X)
mean(varparthi$varpart.X)

varparthi2<-calc.varpart(mod.hi9)
varparthi2$varpart.X[1:10]
varparthi2$varpart.lv[1:10]

hist(varparthi2$varpart.X)
mean(varparthi2$varpart.X)

#The above varpart is the contribution of fixed vs latent variables to explained variation, it is not an R2.
#I want to also try to fit some simple linear models to see how much environment matters to some of the species.
hmscXhi2
hmscYhi4

temp<-vector(length=396)
for (i in 1:396){
  mhi<-lm(hmscYhi4[,i]~hmscXhi2)
  temp[i]<-summary(mhi)$adj.r.squared
}
mean(temp)

mhi<-lm(hmscYhi4[,7]~hmscXhi2)
summary(mhi)
hist(resid(mhi))






##### Environmental correlations #####
envcor.lolv4occ9exp4f<-get.enviro.cor(fit.lolv4occ9exp4f)
corrplot(envcor.lolv4occ9exp4f$cor[1:100,1:100], type="lower", diag=F, title="Environmental correlations", mar=c(3,0.5,2,1), tl.srt=45,tl.pos="n")
corrplot(envcor.lolv4occ9exp4f$sig.cor[1:100,1:100], type="lower", diag=F, title="Environmental correlations", mar=c(3,0.5,2,1), tl.srt=45,tl.pos="n")




##Extract model coefficients
cbind(fit.lolv4occ10exp4$lv.coefs.mean,fit.lolv4occ10exp4$X.coefs.mean)
fit.hilv4occ9exp4f$X.coefs.mean

##Dunn-Smyth residual plots to check model assumption, outliers etc. The first plot should not have a funnel
plot(mod.lo11flv3)
plot(mod.me11flv3)
plot(mod.hi11flv3)
plot(mod.hi9pois)
plot(mod.hi9)

lof.resid<-ds.residuals(mod.lo11flv3)
str(lof.resid)
hist(lof.resid$residuals[,100])

mef.resid<-ds.residuals(mod.me11flv3)
str(mef.resid)
hist(mef.resid$residuals[,100])

hif.resid<-ds.residuals(mod.hi11flv3)
dim(hif.resid$residuals)
hist(hif.resid$residuals[,10])

## Residual ordination plot of the sites (please see Figure 2b in the main text for relevant explanation)
## Please note that boral currently does not automatically produce biplots, although species loadings can be obtained from fit.lvm$lv.coefs.median[,2:3]
lvsplot(mod.mef9lv4)



##### Extract number of significant correlations #####

(length(which(rescor.lo11flv3$sig.correlaton!=0))-dim(rescor.lo11flv3$sig.correlaton)[1])/2
(length(which(rescor.me11flv3$sig.correlaton!=0))-dim(rescor.me11flv3$sig.correlaton)[1])/2
(length(which(rescor.hi11flv3$sig.correlaton!=0))-dim(rescor.hi11flv3$sig.correlaton)[1])/2

#The analysis in MovingUpHill3 suggested that more mcmc iterations meant fewer significant interactions, but I'm not finding that to matter this time.

###### Use corrplot package to plot residual correlations between species #####

corrplot(rescor.lo$sig.correlaton[1:100,1:100], diag=F, type="lower", title="Residual correlations", mar=c(3,0.5,2,1), tl.srt = 45,method = "color")#
corrplot(rescor.hi$sig.correlaton[1:100,1:100], type="lower", diag=F, title="Residual correlations", mar=c(3,0.5,2,1), tl.srt=45,method = "color")
corrplot(rescor.melv5occ9exp4$sig.correlaton[1:100,1:100], type="lower", diag=F, title="Residual correlations", mar=c(3,0.5,2,1), tl.srt=45,method = "color")




##### Plotting with igraph #####

##### colors #####

labelcols<-data.frame(rbind(c("Bacteria","#7879BC"),
                            c("Eukaryota","#94BA3C"),
                            c("Mesofauna","#ff9c34"),
                            c("Fungi","#F6EC32"),
                            c("Plant","#E95275")))
colnames(labelcols)=c("group","color")


# including photosynthetic/not information
labelcols<-data.frame(rbind(c("PhotosyntheticEukaryota","#49874c"),# 466D24
                            c("HeterotrophicEukaryota","#673482"),
                            c("PhotosyntheticBacteria","#94BA3C"),
                            c("HeterotrophicBacteria","#7879BC"),
                            c("ChemoautotrophicBacteria","#6295cd"),
                            c("UnknownEukaryota","gray50"),
                            c("UnknownBacteria","gray70"),
                            c("Mesofauna","#ff9c34"),
                            c("Fungi","#F6EC32"),
                            c("Plant","#E95275")))
colnames(labelcols)=c("group2","color")


head(labelfile)
unique(labelfile$group2)

labelsall<-merge(labelfile,labelcols,"group2",all.x=F,all.y=F) #"labels"
labelsall$color<-as.character(labelsall$color)
head(labelsall)
labelsall$group2<-factor(labelsall$group2,levels=c("HeterotrophicBacteria","PhotosyntheticBacteria","ChemoautotrophicBacteria","UnknownBacteria","Fungi","HeterotrophicEukaryota","PhotosyntheticEukaryota","UnknownEukaryota","Mesofauna","Plant"))


pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/legend2.pdf")
plot(c(1,1),c(1,1))
#legend("topright",c("Bacteria","Small Eukaryota","Mesofauna","Fungi","Plants"),pt.bg=c("#7879BC","#94BA3C","#ff9c34","#F6EC32","#E95275"),bty="n",pch=21,cex=1.4)
legend("topright",c("Heterotrophic bacteria","Photosynthetic bacteria","Chemoautotrophic bacteria","Unknown bacteria","Heterotrophic Eukaryota","Photosynthetic Eukaryota","UnknownEukaryota","Mesofauna","Fungi","Plants"),pt.bg=c("#7879BC","#94BA3C","#6295cd","gray70","#673482","#466D24","gray50","#ff9c34","#F6EC32","#E95275"),bty="n",pch=21,cex=1.4)
legend("topleft",c("Positive","Negative"),col=c("#ce4d42","#687dcb"),lty=1,bty="n",cex=1.4)
#legend("top",as.character(1:10),col=c("#111110","#660011","#A80013","#118877","#4c3d3e","#118877","#7f783f","#aa8888","#aabbdd","#ff99a4"),lty=1,lwd=3)
dev.off()
 
#colors from nico: 111110,660011,112288,A80013,4c3d3e,118877,7f783f,aa8888,aabbdd,ff99a4, Ffccd1,ddd7d7,d8d3ad,e5001a



##### Network diagrams for lo #####
#creating sparse matrix
colMatlo<-rescor.lo11flv3$sig.correlaton
colMatlo[which(rescor.lo11flv3$sig.correlaton>0)]<-1
colMatlo[which(rescor.lo11flv3$sig.correlaton<0)]<- -1

# colMatlo<-rescor.lolv4occ9exp4nosite$sig.correlaton
# colMatlo[which(colMatlo>.85)]<-1
# colMatlo[which(colMatlo<(-.85))]<- -1
# colMatlo[which(colMatlo<.85&colMatlo>(-.85))]<-0

# temp<-colMatlo[upper.tri(colMatlo)]
# temp2<-temp[temp!=0]
# hist(temp2)
# sort(temp2)
# length(temp2)

graphlo1<-graph_from_adjacency_matrix(colMatlo, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelistlo<-data.frame(as_edgelist(graphlo1),weight=E(graphlo1)$weight) #just the edges

length(which(myedgelistlo$weight==1))
length(which(myedgelistlo$weight==-1))
length(which(myedgelistlo$weight==1))/(length(which(myedgelistlo$weight==1))+length(which(myedgelistlo$weight==-1)))

graphlo2<-graph.edgelist(as.matrix(myedgelistlo[,1:2]),directed=FALSE)
graphlo2

verticesgraphlo<-data.frame(otu=rownames(as.matrix(V(graphlo2))))
colorgraphlo<-merge(verticesgraphlo,labelsall,"otu",all.y=F,all.x=F,sort=F)
#shapesgraplo<-ifelse(colorgraphlo$group%in%c("Eukaryota"),"csquare",'circle')

##use colorgraphlo$group2 for ordering, if ther are ties it leaves them in thir original order, thus is still preserves some of the ordering that makes the lines look nice
orderlo<-order(colorgraphlo$group2)
#orderlo<-order(verticesgraphlo$otu)
graphlo2$layout <- layout_in_circle(graphlo2,order=orderlo)
#graphlo2$layout <- layout_in_circle

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/networklocircleposblue.pdf") 
#plot(graphlo2,vertex.size=4,edge.curved=F,edge.color=ifelse(myedgelistlo$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphlo$color,edge.width=.7,vertex.label=NA)#,vertex.shape=shapesgraplo  positive is red
plot(graphlo2,vertex.size=4,edge.curved=F,edge.color=ifelse(myedgelistlo$weight==1,"#687dcb","#ce4d42"),vertex.color=colorgraphlo$color,edge.width=.7,vertex.label=NA)#,vertex.shape=shapesgraplo  positive is blue
#dev.off()

colorgraphlo[which(colorgraphlo$group=="Mesofauna"),]
colorgraphlo[which(colorgraphlo$group=="Fungi"),]
colorgraphlo[which(colorgraphlo$group2=="PhotosyntheticBacteria"),]
colorgraphlo[which(colorgraphlo$group2=="PhotosyntheticEukaryota"),]

temp<-colorgraphlo[which(colorgraphlo$group=="Mesofauna"),"otu"]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp|myedgelistlo$X2%in%temp),]
dim(temp2)






##### Network diagrams for me #####
#creating sparse matrix
colMatme<-rescor.me11flv3$sig.correlaton
colMatme[which(rescor.me11flv3$sig.correlaton>0)]<-1
colMatme[which(rescor.me11flv3$sig.correlaton<0)]<- -1

# colMatme<-rescor.melv4occ9exp4f$sig.correlaton
# colMatme[which(colMatme>.6)]<-1
# colMatme[which(colMatme<(-.6))]<- -1
# colMatme[which(colMatme<.6&colMatme>(-.6))]<-0

graphme1<-graph_from_adjacency_matrix(colMatme, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelistme<-data.frame(as_edgelist(graphme1),weight=E(graphme1)$weight) #just the edges

length(which(myedgelistme$weight==1))
length(which(myedgelistme$weight==-1))
length(which(myedgelistme$weight==1))/(length(which(myedgelistme$weight==1))+length(which(myedgelistme$weight==-1)))

graphme2<-graph.edgelist(as.matrix(myedgelistme[,1:2]),directed=FALSE)
graphme2

verticesgraphme<-data.frame(otu=rownames(as.matrix(V(graphme2))))
colorgraphme<-merge(verticesgraphme,labelsall,"otu",all.y=F,all.x=F,sort=F)

##use colorgraphlo$group2 for ordering, if there are ties it leaves them in their original order, thus is still preserves some of the ordering that makes the lines look nice
orderme<-order(colorgraphme$group2)
#orderme<-order(verticesgraphme$otu)
graphme2$layout <- layout_in_circle(graphme2,order=orderme)
#graphme2$layout <- layout_in_circle

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/networkmecircleposblue.pdf") 
#plot(graphme2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelistme$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphme$color,edge.width=.7)#,layout=l3
plot(graphme2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelistme$weight==1,"#687dcb","#ce4d42"),vertex.color=colorgraphme$color,edge.width=.7)#,layout=l3
#dev.off()

temp<-colorgraphme[which(colorgraphme$group=="Plant"),"otu"]
temp2<-myedgelistme[which(myedgelistme$X1%in%temp|myedgelistme$X2%in%temp),]
temp2
dim(temp2)

temp<-colorgraphme[which(colorgraphme$group=="Mesofauna"),"otu"]
temp2<-myedgelistme[which(myedgelistme$X1%in%temp|myedgelistme$X2%in%temp),]
temp2
dim(temp2)

colorgraphme[which(colorgraphme$otu=="Bf8ab7e424f6976c81b44b3c809dc6ce7"),]
colorgraphme[which(colorgraphme$otu=="Bdc4a7fab972ac91dd37631d279420a08"),]




##### Network diagrams for hi #####
#creating sparse matrix
colMathi<-rescor.hi11flv3$sig.correlaton#rescor.hi9fixed0
colMathi[which(rescor.hi11flv3$sig.correlaton>0)]<-1
colMathi[which(rescor.hi11flv3$sig.correlaton<0)]<- -1

# colMathi<-rescor.hi9$correlation
# colMathi[which(colMathi>.5)]<-1
# colMathi[which(colMathi<(-.5))]<- -1
# colMathi[which(colMathi<.5&colMathi>(-.5))]<-0

graphhi1<-graph_from_adjacency_matrix(colMathi, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelisthi<-data.frame(as_edgelist(graphhi1),weight=E(graphhi1)$weight) #just the edges

length(which(myedgelisthi$weight==1))
length(which(myedgelisthi$weight==-1))
length(which(myedgelisthi$weight==1))/(length(which(myedgelisthi$weight==1))+length(which(myedgelisthi$weight==-1)))

graphhi2<-graph.edgelist(as.matrix(myedgelisthi[,1:2]),directed=FALSE)
graphhi2

verticesgraphhi<-data.frame(otu=rownames(as.matrix(V(graphhi2))))
colorgraphhi<-merge(verticesgraphhi,labelsall,"otu",all.y=F,all.x=F,sort=F)

#order starts at 3:00 and goes counterclockwise
##use colorgraphlo$group2 for ordering, if ther are ties it leaves them in thir original order, thus is still preserves some of the ordering that makes the lines look nice
orderhi<-order(colorgraphhi$group2)
#orderhi<-order(verticesgraphhi$otu)
graphhi2$layout <- layout_in_circle(graphhi2,order=orderhi)
#graphhi2$layout <- layout_in_circle

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/networkhicircleposblue.pdf") 
#plot(graphhi2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelisthi$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphhi$color,edge.width=.7)#,layout=l3  
plot(graphhi2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelisthi$weight==1,"#687dcb","#ce4d42"),vertex.color=colorgraphhi$color,edge.width=.7)#,layout=l3  
#dev.off()

colorgraphhi[which(colorgraphhi$group=="Mesofauna"),]
colorgraphhi[which(colorgraphhi$group=="Plant"),]

temp<-colorgraphhi[which(colorgraphhi$group=="Plant"),"otu"]
temp2<-myedgelisthi[which(myedgelisthi$X1%in%temp|myedgelisthi$X2%in%temp),]
temp2
dim(temp2)

#Creating a subgraph
colorgraphhi2<-colorgraphhi[which(colorgraphhi$group2=="Mesofauna"),]
myedgelisthi2<-myedgelisthi[which(myedgelisthi[,"X1"]%in%colorgraphhi2$otu|myedgelisthi[,"X2"]%in%colorgraphhi2$otu),]
graph3<-subgraph.edges(graphhi2, eids=which(myedgelisthi2[,"X1"]%in%colorgraphhi2$otu|myedgelisthi2[,"X2"]%in%colorgraphhi2$otu), delete.vertices = F)
plot(graph3,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelisthi2$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphhi$color)#,rescale=F,xlim=c(-1,1),ylim=c(-1,1)





###### Network statistics #######

length(E(graphlo2))/length(V(graphlo2))
length(E(graphme2))/length(V(graphme2))
length(E(graphhi2))/length(V(graphhi2))

temp<-colorgraphlo; temp$ones<-1
aggregate.data.frame(temp$ones,by=list(temp$group2),sum)

temp<-colorgraphme; temp$ones<-1
aggregate.data.frame(temp$ones,by=list(temp$group2),sum)

temp<-colorgraphhi; temp$ones<-1
aggregate.data.frame(temp$ones,by=list(temp$group2),sum)




##### Photosynthetic microbes #####
#connections with Ps bacteria/eukaryotes and bacterial heterotrophs

temp<-colorgraphlo[which(colorgraphlo$group2=="PhotosyntheticBacteria"|colorgraphlo$group2=="PhotosyntheticEukaryota"),"otu"]
length(temp)
outputPsinteractionlo<-myedgelistlo[which(myedgelistlo$X1%in%temp|myedgelistlo$X2%in%temp),]
head(outputPsinteractionlo)

outputPsinteractionlo$ID1<-NA
outputPsinteractionlo$ID2<-NA

for (i in 1:dim(outputPsinteractionlo)[1]){
  sp1<-outputPsinteractionlo$X1[i]
  outputPsinteractionlo$ID1[i]<-as.character(colorgraphlo$group2[which(colorgraphlo$otu==as.character(sp1))])
  sp2<-outputPsinteractionlo$X2[i]
  outputPsinteractionlo$ID2[i]<-as.character(colorgraphlo$group2[which(colorgraphlo$otu==as.character(sp2))])
}
ind1<-which(outputPsinteractionlo$ID1=="PhotosyntheticBacteria"&outputPsinteractionlo$ID2=="HeterotrophicBacteria")
ind2<-which(outputPsinteractionlo$ID1=="HeterotrophicBacteria"&outputPsinteractionlo$ID2=="PhotosyntheticBacteria")
ind3<-which(outputPsinteractionlo$ID1=="PhotosyntheticEukaryota"&outputPsinteractionlo$ID2=="HeterotrophicBacteria")
ind4<-which(outputPsinteractionlo$ID1=="HeterotrophicBacteria"&outputPsinteractionlo$ID2=="PhotosyntheticEukaryota")
length(ind1)+length(ind2)+length(ind3)+length(ind4)

outputPsinteractionlo[ind1,]

#ps with eachother
ind1<-which(outputPsinteractionlo$ID1=="PhotosyntheticBacteria"&outputPsinteractionlo$ID2=="PhotosyntheticEukaryota")
ind2<-which(outputPsinteractionlo$ID1=="PhotosyntheticEukaryota"&outputPsinteractionlo$ID2=="PhotosyntheticBacteria")


temp<-colorgraphme[which(colorgraphme$group2=="PhotosyntheticBacteria"|colorgraphme$group2=="PhotosyntheticEukaryota"),"otu"]

outputPsinteractionme<-myedgelistme[which(myedgelistme$X1%in%temp|myedgelistme$X2%in%temp),]
head(outputPsinteractionme)

outputPsinteractionme$ID1<-NA
outputPsinteractionme$ID2<-NA

for (i in 1:dim(outputPsinteractionme)[1]){
  sp1<-outputPsinteractionme$X1[i]
  outputPsinteractionme$ID1[i]<-as.character(colorgraphme$group2[which(colorgraphme$otu==as.character(sp1))])
  sp2<-outputPsinteractionme$X2[i]
  outputPsinteractionme$ID2[i]<-as.character(colorgraphme$group2[which(colorgraphme$otu==as.character(sp2))])
}
ind1<-which(outputPsinteractionme$ID1=="PhotosyntheticBacteria"&outputPsinteractionme$ID2=="HeterotrophicBacteria")
ind2<-which(outputPsinteractionme$ID1=="HeterotrophicBacteria"&outputPsinteractionme$ID2=="PhotosyntheticBacteria")
ind3<-which(outputPsinteractionme$ID1=="PhotosyntheticEukaryota"&outputPsinteractionme$ID2=="HeterotrophicBacteria")
ind4<-which(outputPsinteractionme$ID1=="HeterotrophicBacteria"&outputPsinteractionme$ID2=="PhotosyntheticEukaryota")
length(ind1)+length(ind2)+length(ind3)+length(ind4)

outputPsinteractionme[ind1,]


temp<-colorgraphhi[which(colorgraphhi$group2=="PhotosyntheticBacteria"|colorgraphhi$group2=="PhotosyntheticEukaryota"),"otu"]
outputPsinteractionhi<-myedgelisthi[which(myedgelisthi$X1%in%temp|myedgelisthi$X2%in%temp),]
head(outputPsinteractionhi)

outputPsinteractionhi$ID1<-NA
outputPsinteractionhi$ID2<-NA

for (i in 1:dim(outputPsinteractionhi)[1]){
  sp1<-outputPsinteractionhi$X1[i]
  outputPsinteractionhi$ID1[i]<-as.character(colorgraphhi$group2[which(colorgraphhi$otu==as.character(sp1))])
  sp2<-outputPsinteractionhi$X2[i]
  outputPsinteractionhi$ID2[i]<-as.character(colorgraphhi$group2[which(colorgraphhi$otu==as.character(sp2))])
}
ind1<-which(outputPsinteractionhi$ID1=="PhotosyntheticBacteria"&outputPsinteractionhi$ID2=="HeterotrophicBacteria")
ind2<-which(outputPsinteractionhi$ID1=="HeterotrophicBacteria"&outputPsinteractionhi$ID2=="PhotosyntheticBacteria")
ind3<-which(outputPsinteractionhi$ID1=="PhotosyntheticEukaryota"&outputPsinteractionhi$ID2=="HeterotrophicBacteria")
ind4<-which(outputPsinteractionhi$ID1=="HeterotrophicBacteria"&outputPsinteractionhi$ID2=="PhotosyntheticEukaryota")
length(ind1)+length(ind2)+length(ind3)+length(ind4)

outputPsinteractionhi[ind2,]



#counting total interations with Ps microbes, not knowing who the connection is with or pos/neg
temp<-colorgraphlo[which(colorgraphlo$group2=="PhotosyntheticBacteria"|colorgraphlo$group2=="PhotosyntheticEukaryota"),"otu"]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp|myedgelistlo$X2%in%temp),]
dim(temp2) #370 interactions

temp<-colorgraphme[which(colorgraphme$group2=="PhotosyntheticBacteria"|colorgraphme$group2=="PhotosyntheticEukaryota"),"otu"]
temp2<-myedgelistme[which(myedgelistme$X1%in%temp|myedgelistme$X2%in%temp),]
temp2
dim(temp2) #148 interactions
colorgraphme[which(colorgraphme$otu%in%temp2$X1|colorgraphme$otu%in%temp2$X2),]

temp<-colorgraphhi[which(colorgraphhi$group2=="PhotosyntheticBacteria"|colorgraphhi$group2=="PhotosyntheticEukaryota"),"otu"]
temp2<-myedgelisthi[which(myedgelisthi$X1%in%temp|myedgelisthi$X2%in%temp),]
dim(temp2) #2 interactions
temp2
colorgraphhi[which(colorgraphhi$otu%in%temp2$X1|colorgraphhi$otu%in%temp2$X2),]



#### Interactions with Ktedonobacteria #####
#Bacteria in networks (lo abundance, hi)
temp<-colorgraphlo[which(colorgraphlo$group=="Bacteria"),]
ind<-grep("Ktedonobacteria",temp$taxstring)
temp1<-temp[ind,]
#temp<-colorgraphlo[which(colorgraphlo$group=="Bacteria"),"otu"]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp1$otu|myedgelistlo$X2%in%temp1$otu),]
head(temp2)
sum(temp2$weight[temp2$weight==1])
dim(temp2) #210 interactions with Ktedonobacteria
colnames(temp2)[1]<-"otu"
temp2<-merge(temp2,colorgraphlo[,c(1,2)])
colnames(temp2)[4]<-"X1c"
colnames(temp2)[1]<-"X1"
colnames(temp2)[2]<-"otu"
temp2<-merge(temp2,colorgraphlo[,c(1,2)])
colnames(temp2)[5]<-"X2c"
colnames(temp2)[1]<-"X2"
head(temp2)
ind<-which(temp2$X1c!="HeterotrophicBacteria")
temp2$X2c[ind]<-temp2$X1c[ind]
temp2$X1c[ind]<-"HeterotrophicBacteria"
#now X2c is all the partner of the ktedonobacteria
temp3<-temp2[which(temp2$weight==1),]
aggregate.data.frame(temp3$weight,by=list(temp3$X2c),sum)
#subtract 7 positive interactions are between two Ktedonobacteria
myedgelistlo[which(myedgelistlo$X1%in%temp1$otu&myedgelistlo$X2%in%temp1$otu),]
dim(myedgelistlo[which(myedgelistlo$X1%in%temp1$otu&myedgelistlo$X2%in%temp1$otu),])
210-7

temp<-colorgraphme[which(colorgraphme$group=="Bacteria"),]
ind<-grep("Ktedonobacteria",temp$taxstring)
temp1<-temp[ind,]
#temp<-colorgraphlo[which(colorgraphlo$group=="Bacteria"),"otu"]
temp2<-myedgelistme[which(myedgelistme$X1%in%temp1$otu|myedgelistme$X2%in%temp1$otu),]
temp2
sum(temp2$weight[temp2$weight==1])
dim(temp2)
colnames(temp2)[1]<-"otu"
temp2<-merge(temp2,colorgraphme[,c(1,2)])
colnames(temp2)[4]<-"X1c"
colnames(temp2)[1]<-"X1"
colnames(temp2)[2]<-"otu"
temp2<-merge(temp2,colorgraphme[,c(1,2)])
colnames(temp2)[5]<-"X2c"
colnames(temp2)[1]<-"X2"
head(temp2)
ind<-which(temp2$X1c!="HeterotrophicBacteria")
temp2$X2c[ind]<-temp2$X1c[ind]
temp2$X1c[ind]<-"HeterotrophicBacteria"
#now X2c is all the partner of the ktedonobacteria
temp3<-temp2[which(temp2$weight==1),]
aggregate.data.frame(temp3$weight,by=list(temp3$X2c),sum)
#subtract 4 positive interactions are between two Ktedonobacteria
dim(myedgelistme[which(myedgelistme$X1%in%temp1$otu&myedgelistme$X2%in%temp1$otu),])
165-4

temp<-colorgraphhi[which(colorgraphhi$group=="Bacteria"),]
ind<-grep("Ktedonobacteria",temp$taxstring)
temp1<-temp[ind,]
#temp<-colorgraphlo[which(colorgraphlo$group=="Bacteria"),"otu"]
temp2<-myedgelisthi[which(myedgelisthi$X1%in%temp1$otu|myedgelisthi$X2%in%temp1$otu),]
temp2
sum(temp2$weight[temp2$weight==1])
dim(temp2)
colnames(temp2)[1]<-"otu"
temp2<-merge(temp2,colorgraphhi[,c(1,2)])
colnames(temp2)[4]<-"X1c"
colnames(temp2)[1]<-"X1"
colnames(temp2)[2]<-"otu"
temp2<-merge(temp2,colorgraphhi[,c(1,2)])
colnames(temp2)[5]<-"X2c"
colnames(temp2)[1]<-"X2"
head(temp2)
ind<-which(temp2$X1c!="HeterotrophicBacteria")
temp2$X2c[ind]<-temp2$X1c[ind]
temp2$X1c[ind]<-"HeterotrophicBacteria"
#now X2c is all the partner of the ktedonobacteria
temp3<-temp2[which(temp2$weight==1),]
aggregate.data.frame(temp3$weight,by=list(temp3$X2c),sum)
#subtract 1 positive interactions are between two Ktedonobacteria
dim(myedgelisthi[which(myedgelisthi$X1%in%temp1$otu&myedgelisthi$X2%in%temp1$otu),])
7-0



###### Relationship with known organisms #####

ind<-grep("Ciliophora",colorgraphlo$taxstring,perl=T,value=F)
colorgraphlo[ind,]
ind<-grep("Ciliophora",colorgraphme$taxstring,perl=T,value=F)
colorgraphme[ind,]
ind<-grep("Ciliophora",colorgraphhi$taxstring,perl=T,value=F)
colorgraphhi[ind,]

#lo

Ciliophora - 4 taxa
S7a97813269020725286a63989363ceef Litostomatea (class). when I blastd this, the closest identified hit was Foissnerides in the haptoria (subclass), Haptorida (an order of ciliate predators) but only 90% identity. The only other thing it could be in the litostomatea is from the group that consists of endosymbionts in th digestive tract of vertebrates which is it not probably, so otherwise it is a predator
S13a0576d1656e575cdc2dc94fd138719 Litostomatea/Isotricha (genus) uncultured_rumen_protozoa
  when I blasted this it came back as 100% idenity to Enchelyodon sp (haptorida). which is a freeliving ciliate predator (https://www.nies.go.jp/chiiki1/protoz/morpho/ciliopho/enchelyo.htm)
Sec88fea6bf21bef8baba7173b475de7d Scuticociliatia, Scuticociliates often feed on bacteria (subclass). "Scuticociliates often feed on bacteria, using complex morphological adaptations to create currents
and filters capable of capturing bacteria and other particles
from the water column or scraping them from hard surfaces"
Se4ae90d4c05ef00793fc6b93fb6a9af7 Halteria (genus) is a filter feeder also called a grazer on bacteria and other things non selective wide range of smallish sizes (halteria is quite small)

grep -C 2 "13a0576d1656e575cdc2dc94fd138719" dna-sequences.fasta

ind<-grep("Ciliophora",colorgraphlo$taxstring,perl=T,value=F)
temp<-colorgraphlo[ind,]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp$otu|myedgelistlo$X2%in%temp$otu),]
temp3<-temp2[temp2$weight==-1,]
colnames(temp3)[2]<-"otu"
temp4<-merge(temp3,colorgraphlo[,c(1,2,5)])
18ngative relationships with bacteria

#the Ochromonadales, a golden algae taht is a mixotroph and can be predatory/heterotrophic and photosynthetic
ind<-grep("Ochromo",colorgraphlo$taxstring,perl=T,value=F)
colorgraphlo[ind,]
temp<-colorgraphlo[which(colorgraphlo$otu%in%c("Sba9b295b95881ce85c632aeb3b4d3fd6","Sc45dd62d47c0618457b55cd98c0a28c1")),]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp$otu|myedgelistlo$X2%in%temp$otu),]
temp2[temp2$weight==-1,]  #7 negative relationships with bacteria


#mid 

#mesofauna
temp<-colorgraphme[which(colorgraphme$group=="Mesofauna"),]
myedgelistme[which(myedgelistme$X1%in%temp$otu|myedgelistme$X2%in%temp$otu),]
colorgraphme[which(colorgraphme$otu=="Bdc4a7fab972ac91dd37631d279420a08"),] #negative relatinoship acidobacteria
colorgraphme[which(colorgraphme$otu=="Be8accbbab64dba2fe98edb7072796cc5"),] #positive relationship

#plants
temp<-colorgraphme[which(colorgraphme$group=="Plant"),]
temp2<-myedgelistme[which(myedgelistme$X1%in%temp$otu|myedgelistme$X2%in%temp$otu),]
colnames(temp2)[1]<-"otu"
temp3<-merge(temp2,colorgraphme[,c(1,2,5)])
temp3[which(temp3$weight==1),]
temp3[which(temp3$weight==-1),]


#high
temp<-colorgraphhi[which(colorgraphhi$group=="Plant"),]
temp2<-myedgelisthi[which(myedgelisthi$X1%in%temp$otu|myedgelisthi$X2%in%temp$otu),]
colnames(temp2)[1]<-"otu"
temp3<-merge(temp2,colorgraphhi[,c(1,2,5)])
temp3[which(temp3$weight==1),]

#look up sequence for chitinophagaceeae
grep -C 2 "2d66d59115e5f8c6d716f72f4fc26f23" dna-sequences.fasta











##### try transforming the data with a clr #####

#started with comm.bio, took each data set apart, added 1 to the 0s, calculated the clr (on everything except plants), put the data back together, then subset into lo/hi and deleted rare taxa
colnames(hmscY)[7208:7262]

hmscYN<-hmscY[,1:142]; hmscYN[hmscYN==0]<-1
hmscYN2<-t(apply(hmscYN, 1, function(x){log(x) - mean(log(x))}))

cbind(hmscYN[,1],hmscYN2[,1])
hist(hmscYN[,9])
hist(hmscYN2[,9])

hmscYS<-hmscY[,143:1233]; hmscYS[hmscYS==0]<-1
hmscYS2<-t(apply(hmscYS, 1, function(x){log(x) - mean(log(x))}))

hmscYBac<-hmscY[,1234:6086]; hmscYBac[hmscYBac==0]<-1
hmscYBac2<-t(apply(hmscYBac, 1, function(x){log(x) - mean(log(x))}))

hmscYITS<-hmscY[,6087:7208]; hmscYITS[hmscYITS==0]<-1
hmscYITS2<-t(apply(hmscYITS, 1, function(x){log(x) - mean(log(x))}))

hmscYPlant<-hmscY[,7209:7262]

hmscY2<-cbind(hmscYN2,hmscYS2,hmscYBac2,hmscYITS2,hmscYPlant)
hmscY2[1:5,1:5]




##### to remove all the many trial files to make env smaller #####

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
                        "biogeo6",
                        "comm.bio",
                        "mod.lo11flv3", 
                        "mod.me11flv3", 
                        "mod.hi11flv3",
                        "rescor.lo11flv3",
                        "rescor.me11flv3",
                        "rescor.hi11flv3")))


