## Fitting a Latent Variable Model using the boral package ##
## Use boral to fit a LVM with 4 latent variables, a random site effect, and the effect of fixed environmental variables

## NOTE: boral uses JAGS so you need to first install JAGS-3.0.0.exe or higher from http://sourceforge.net/projects/mcmc-jags/files/JAGS/

library(boral) #Need version 0.7 or later, available on CRAN.
library(Matrix)


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

#Merge the biogeo6 with comm.dataALL, then split them to make sure the same samples and order are in each dataset
comm.bio<-merge(biogeo6,comm.dataALL5)
comm.bio[1:10,1:60]

#the comm.bio.csv file that is saved here is from the old bioinformatics. I saved it b/c the R environment file takes so long to load. however I might want to reinstate this, if this environment gets really big
#write.csv(comm.bio,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/comm.bio.csv",row.names=F)
#comm.bio<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/comm.bio.csv")





##### Split into lo/me/hi #####

## lo ##
dim(comm.bio)

ind<-which(comm.bio$lomehi=="lo")
  
hmscYlo<-comm.bio[ind,55:24565]  #start at 54 if you want lomehi column
rownames(hmscYlo)<-comm.bio$X.SampleID[ind]
hmscYlo[1:10,1:10]

hmscXlo<-data.frame(snowdepth=comm.bio$snowdepth,TC=comm.bio$TC,pH=comm.bio$pH,moisture=comm.bio$moisture,lomehi=comm.bio$lomehi)[ind,] #,plantcov=comm.bio$plantcov,whc=comm.bio$WHC
rownames(hmscXlo)<-comm.bio$X.SampleID[ind]

rcorr(as.matrix(hmscXlo[,2:(dim(hmscXlo)[2]-1)]))
#plot(hmscX$TC,hmscX$moisture)

hmscYct <- cmultRepl(hmscYc,label=0, method="CZM")
hmscYct[1:5,1:5]



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


#######




#select lo/me/hi
ind<-which(hmscX$lomehi=="hi")
hmscXb<-hmscX[ind,]
hmscYb<-hmscY[ind,]
hmscY2b<-hmscY2[ind,]

hmscYb[1:5,1:5]
hmscY2b[1:5,1:5]

#select species with greater than X (X+1 or more) occurrences and remove lo me hi (since you can't have text in a matrix or the whole matrix becomes character)
ind<-which(colSums(hmscYb>0)>8)
length(ind)
#hmscYc<-hmscYb[,ind]
hmscYc<-hmscY2b[,ind]
hmscXc<-hmscXb[,1:dim(hmscXb)[2]-1]#
dim(hmscYc)
dim(hmscXc)

hmscYc[1:5,1:5]

#the y data are not normal (the only options I have are normal, binary, poisson, overdispersed poisson), so I could do a sqrt transformation on Y (log(0) is -Inf). log(x+1) doesn't work since the proportions are so low, could do log(x*100+1) but the sqrt actually makes it more normal
#hmscYd<-sqrt(hmscYc*100)
#hmscYd<-hmscYc*100 #this is odd, 9 taxa have negative R2
#hmscYd<-log(hmscYc+1)
hist(hmscYc[,101])
#hist(hmscYd[,12])

#check if the values are too low that some tolerance is messing up the CI estimates, yes important to scale y, instead of scaling Y I will use the sqrt transform of the percent. I will scale x, since they differ so much in range
hmscXd<-scale(hmscXc)
hmscXd[,1]<-1

#hmscYd2<-scale(hmscYd)

####do randomization#####
hmscYr<-randomizeMatrix(hmscYc, null.model="independentswap",iterations=30000) #I upped to 30000 from 1000 bc it seems like things weren't getting randomized
hmscYc<-hmscYr
#######

#make them matrices
hmscXe<-as.matrix(hmscXd)

hmscYe<-as.matrix(hmscYc)
#hmscYe<-as.matrix(hmscYd2)

dim(hmscYe)
dim(hmscXe)


## Covariates need to be stored as a matrix for boral. no need for intercept
covX <- hmscXe[,-1]



#hi and lo random dataframes
hmscYehi<-hmscYe
covXhi<-covX
hmscYelo<-hmscYe
covXlo<-covX
hmscYehir<-hmscYe
covXhir<-covX
hmscYelor<-hmscYe
covXlor<-covX


##### Fit the LVM using boral and calculate residual correlation matrix#####

#List of files produced:
#fit.hilv4occ10exp4
#fit.lolv4occ10exp4
#fit.hilv4occ9exp4
#fit.lolv4occ9exp4
#fit.hilv4occ8exp4, do not have the corresponding rescor for this b/c my computer crashed
#fit.melv4occ9exp4
#fit.melv3occ9exp4
#fit.melv5occ9exp4
#fit.melv6occ9exp4
#fit.melv2occ9exp4
#fit.lolv4occ9exp4f - f means final model fitted with long mcmc chain, start 2:00pm, model finished 4:20pm, rescor finished 
#fit.melv4occ9exp4f - f means final model fitted with long mcmc chain
#fit.hilv4occ9exp4f - f means final model fitted with long mcmc chain
#fit.lolv4occ9exp4nosite - with no site random effect, fit with long chain
#fit.lolv4occ9exp0nosite - with no site random effct and only latent vaiables
#fit.melv4occ9exp4nosite
#fit.hilv4occ9exp4nosite


#Using the default mcmc parameters, the models take about 2 hrs to fit.  mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123)
#Using shorter chains, it takes about 12 min to fit.  mcmc.control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 3, seed = 123)
# in the tutorial they add this to the model fitting code, hypparams = c(20,20,20,20), however, now if you wanted to change this you need to put it in a prior.control statement or something. I am just using the default here

#hir <- boral(y = hmscYe, X = covX, lv.control = list(num.lv = 4), family = "negative.binomial", row.eff = "random", save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 3, seed = 123))

#fitting models for comparing percent covariance explained by env, take row.eff out
#it looks lik this is what I used for my final models in the second ISME submission. the results are very different from not using the site effect. In the Hui 2016 paper he says that the row effect adjusts for differnces in site total abundance so that, so if you want relative abundance include it. Except the other thing is that he only uses it in latent-only variable models, not models with environmental variables: "Note that we have chosen to remove the row effect, as we are now using the environmental covariates directly to also explain the differences in site total abundance". I *think* my thinking before was that I've already relativized my data, so I don't want it relativized again & also that I didn't want to standardize across different groups of organisms plants/bacteria/fungi. However, I am still very suprized that the differences between the models are so huge. it probably has to do with on sampl having a higher total than another so things are going to be mor positivly correlatd, but this is real, b/c if w included all taxa, they would have those abundances and be positively correlated.
#fit.hilv4occ9exp4nosite <- boral(y = hmscYe, X = covX, num.lv = 4, family = "negative.binomial", save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))
#fit.lolv4occ9exp0nosite <- boral(y = hmscYe, X = NULL, num.lv = 4, family = "negative.binomial", save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))
hinosite <- boral(y = hmscYehi, X = covXhi, lv.control = list(num.lv = 4), family = "normal", save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 30, seed = 123))
lonosite <- boral(y = hmscYelo, X = covXlo, lv.control = list(num.lv = 4), family = "normal", save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 30, seed = 123))
hirnosite <- boral(y = hmscYehir, X = covXhir, lv.control = list(num.lv = 4), family = "negative.binomial", save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 30, seed = 123))
lornosite <- boral(y = hmscYelor, X = covXlor, lv.control = list(num.lv = 4), family = "negative.binomial", save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 30, seed = 123)) #start 9:05, end 9:15


#Calculate residual correlation matrix
#Using
#Using shorter chains (see above comment), takes about 11 min to fit
#This will not calculate with a cuttoff of occ6, R studio crashes on my laptop

# rescor.lolv4occ9exp4nosite <- get.residual.cor(fit.lolv4occ9exp4nosite) 
# rescor.lolv4occ9exp0nosite <- get.residual.cor(fit.lolv4occ9exp0nosite) 
# rescor.melv4occ9exp4nosite <- get.residual.cor(fit.melv4occ9exp4nosite) 
# rescor.hilv4occ9exp4nosite <- get.residual.cor(fit.hilv4occ9exp4nosite) 
rescor.hirnosite <- get.residual.cor(hirnosite) 
rescor.lornosite <- get.residual.cor(lornosite) 
rescor.lonosite <- get.residual.cor(lonosite) #start 9:15, end 9:20
rescor.hinosite <- get.residual.cor(hinosite) 


#Use modified function below to get 90% CI
#res.corshiocc8.90<-get.residual.cor2(fit.lvmhiocc8)

save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill4_Workspace_Analysis3.Rdata")  


##### Look at results and check convergence/fit #####

#Model fit
fit.melv2occ9exp4$ics[1]
fit.melv3occ9exp4$ics[1]
fit.melv4occ9exp4$ics[1]
fit.melv5occ9exp4$ics[1]
fit.melv6occ9exp4$ics[1]
i<-1
plot(2:6,c(fit.melv2occ9exp4$ics[i],
           fit.melv3occ9exp4$ics[i],
           fit.melv4occ9exp4$ics[i],
           fit.melv5occ9exp4$ics[i],
           fit.melv6occ9exp4$ics[i]),type = "b")


summary(fit.hilv4occ9exp4f) # To look at estimated parameter values
fit.lolv4occ10exp4$hpdintervals # 95% credible intervals for model parameters.

#check information criteria
fit.lolv4occ10exp4$ics

#Check convergence
#Geweke diagnostic - a z test testing whether the first 10% and the last 50% are diffrent (i think those are the fractions, doesn't really matter exactly), if it is significant, then the means are different and it didn't converge
plot(get.mcmcsamples(fit.hilv4occ9exp4)[,1])
plot(get.mcmcsamples(fit.hilv4occ9exp4f)[,1])
plot(get.mcmcsamples(lornosite)[,2])
plot(get.mcmcsamples(hirnosite)[,2])

#the order is effct of snowdepth for each of th 600 species, then effect of TC, then pH, then moisture 
mcmchi<-get.mcmcsamples(fit.lolv4occ9exp4f)
dim(mcmchi)
colnames(mcmchi)[2500:3000]
mcmchi[1:10,1:5]

#TRUE means these did not converge
gew.pvals <- 2*pnorm(abs(unlist(fit.hilv4occ9exp4f$geweke.diag[[1]])), lower.tail = FALSE)
length(gew.pvals)
gew.pvals[1:5]
gew.pvals[which(gew.pvals<.05)] #technically these did not converge, however, the trace plots look fine to me
p.adjust(gew.pvals, method = "holm")

fit.hilv4occ9exp4f$geweke.diag
lornosite$geweke.diag
hirnosite$geweke.diag
lonosite$geweke.diag 
hinosite$geweke.diag 

#example of one that did not converge
#(1st species) N6f914ead2160e51670d3dc70c25e107b for snowdepth did not converge, but looking at the trace plot, it seems fine
#geweke diagnostic
fit.hilv4occ9exp4f$geweke.diag$geweke.diag$X.coefs[1:5,]
#trace plot (it is the very first parameter)
plot(get.mcmcsamples(fit.hilv4occ9exp4f)[,1])
#mean of the mcmc chain to make sure I'm looking at the right parameter
mean(get.mcmcsamples(fit.hilv4occ9exp4f)[,1]) #mean is -1.710607
#mean of the extractd model coefficients (to mak sure I'm looking at the right parameter)
fit.hilv4occ9exp4f$X.coefs.mean  #-1.710607072, yes checks



##### Percent covariation explained by env #####
#One approach to quantify how much of the species co-occurrence is explained by covariates (that is, how well the predictor variables describe the species assemblage) is through differences in the trace of the estimated residual covariance matrix induced by the latent variables (Warton et al. 2015). From the above code, this can be obtained as rescors$trace. For the spider data set, when we compared a pure latent variable model (similar to the equation 1 but without site effects) to the correlated response model, the trace decreased from 178.92 to 107.92. This implies that environmental covariates accounted for approximately 40% of the covariation between species.
(107.92-178.92)/178.92

#Need to fit a model with only latent variables
rescor.lolv4occ9exp4f$trace
rescor.lolv4occ9exp0f$trace






##### Environmental correlations #####
envcor.lolv4occ9exp4f<-get.enviro.cor(fit.lolv4occ9exp4f)
corrplot(envcor.lolv4occ9exp4f$cor[1:100,1:100], type="lower", diag=F, title="Environmental correlations", mar=c(3,0.5,2,1), tl.srt=45,tl.pos="n")
corrplot(envcor.lolv4occ9exp4f$sig.cor[1:100,1:100], type="lower", diag=F, title="Environmental correlations", mar=c(3,0.5,2,1), tl.srt=45,tl.pos="n")




##Extract model coefficients
cbind(fit.lolv4occ10exp4$lv.coefs.mean,fit.lolv4occ10exp4$X.coefs.mean)
fit.hilv4occ9exp4f$X.coefs.mean

##Dunn-Smyth residual plots to check model assumption, outliers etc. The first plot should not have a funnel
plot(fit.lolv4occ9exp4f)
plot(fit.melv4occ9exp4f)
plot(fit.hilv4occ9exp4f)
plot(lonosite)
plot(hinosite)

## Residual ordination plot of the sites (please see Figure 2b in the main text for relevant explanation)
## Please note that boral currently does not automatically produce biplots, although species loadings can be obtained from fit.lvm$lv.coefs.median[,2:3]
lvsplot(fit.lolv4occ10exp4)



##### Extract number of significant correlations #####

(length(which(rescor.melv3occ9exp4$sig.correlaton!=0))-dim(rescor.melv3occ9exp4$sig.correlaton)[1])/2
(length(which(rescor.melv4occ9exp4$sig.correlaton!=0))-dim(rescor.melv4occ9exp4$sig.correlaton)[1])/2
(length(which(rescor.melv5occ9exp4$sig.correlaton!=0))-dim(rescor.melv5occ9exp4$sig.correlaton)[1])/2

(length(which(rescor.lolv4occ9exp4f$sig.correlaton!=0))-dim(rescor.lolv4occ9exp4f$sig.correlaton)[1])/2
(length(which(rescor.melv4occ9exp4f$sig.correlaton!=0))-dim(rescor.melv4occ9exp4f$sig.correlaton)[1])/2
(length(which(rescor.hilv4occ9exp4f$sig.correlaton!=0))-dim(rescor.hilv4occ9exp4f$sig.correlaton)[1])/2
(length(which(rescor.hir$sig.correlaton!=0))-dim(rescor.hir$sig.correlaton)[1])/2

#strange - using more mcmc iterations changs how many significant interactions there are (ex: lo: 3000 iter - 3398 interactions, 30000 iter - 1806 interactions). the only way I can explain this is that thre is more error in th paramter estimates with the short chain (i.e. th mixing is wider and the density histogram is wider), thus if you are more confident in th location of the species in ordination space, thn the locations will not overlap as much, and you will say that thy are less correlated.


###### Use corrplot package to plot residual correlations between species #####

corrplot(rescor.melv3occ9exp4$sig.correlaton[1:100,1:100], diag = F, type = "lower", title = "Residual correlations from LVM", mar=c(3,0.5,2,1), tl.srt = 45,method = "color")#
corrplot(rescor.melv4occ9exp4$sig.correlaton[1:100,1:100], type="lower", diag=F, title="Residual correlations", mar=c(3,0.5,2,1), tl.srt=45,method = "color")
corrplot(rescor.melv5occ9exp4$sig.correlaton[1:100,1:100], type="lower", diag=F, title="Residual correlations", mar=c(3,0.5,2,1), tl.srt=45,method = "color")

corrplot(rescor.lolv4occ9exp4f$sig.correlaton[1:200,1:200], type="lower", diag=F, title="Residual correlations", mar=c(3,0.5,2,1), tl.srt=45,method = "color",tl.pos="n")




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

#labelsall<-merge(labelfile,labelcols,"group",all.x=F,all.y=F) #"labels"
labelsall<-merge(labelfile,labelcols,"group2",all.x=F,all.y=F) #"labels"
labelsall$color<-as.character(labelsall$color)
head(labelsall)
labelsall$group2<-factor(labelsall$group2,levels=c("HeterotrophicBacteria","PhotosyntheticBacteria","ChemoautotrophicBacteria","UnknownBacteria","Fungi","HeterotrophicEukaryota","PhotosyntheticEukaryota","UnknownEukaryota","Mesofauna","Plant"))

#old                            
# c("NonphotosyntheticEukaryota","#673482"),
# c("PhotosyntheticEukaryota","#466D24"),
# c("Fungi","#F6EC32"),
# c("Metazoa","#ff9c34"),
# c("Plant","#E95275"),
# c("PhotosyntheticBacteria","#94BA3C"),
# c("NonphotosyntheticBacteria","#7879BC"),
# c("OtherMetazoa","black"),
# c("Nematoda","gray30"),
# c("Tardigrada","gray55"),
# c("Rotifera","gray80"),
# c("Arthropoda","white"),
# c("AF","red"),
# c("AP","red"),
# c("BF","black"),
# c("FF","gray30"),
# c("OM","gray55"),
# c("PP","gray80"),
# c("RA","white"),
# c("unknown","red")))

#colnames(labelcols)=c("labels","color")

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/legend2.pdf")
plot(c(1,1),c(1,1))
#legend("topright",c("Bacteria","Small Eukaryota","Mesofauna","Fungi","Plants"),pt.bg=c("#7879BC","#94BA3C","#ff9c34","#F6EC32","#E95275"),bty="n",pch=21,cex=1.4)
legend("topright",c("Heterotrophic bacteria","Photosynthetic bacteria","Chemoautotophic bacteria","Unknown bacteria","Heterotrophic Eukaryota","Photosynthetic Eukaryota","UnknownEukaryota","Mesofauna","Fungi","Plants"),pt.bg=c("#7879BC","#94BA3C","#6295cd","gray70","#673482","#466D24","gray50","#ff9c34","#F6EC32","#E95275"),bty="n",pch=21,cex=1.4)
legend("topleft",c("Positive","Negative"),col=c("#ce4d42","#687dcb"),lty=1,bty="n",cex=1.4)
#legend("top",as.character(1:10),col=c("#111110","#660011","#A80013","#118877","#4c3d3e","#118877","#7f783f","#aa8888","#aabbdd","#ff99a4"),lty=1,lwd=3)
dev.off()
 
#colors from nico: 111110,660011,112288,A80013,4c3d3e,118877,7f783f,aa8888,aabbdd,ff99a4, Ffccd1,ddd7d7,d8d3ad,e5001a



##### lo #####
#creating sparse matrix
colMatlo<-rescor.lolv4occ9exp4f$sig.correlaton
colMatlo[which(rescor.lolv4occ9exp4f$sig.correlaton>0)]<-1
#colMatlo[which(rescor.lolv4occ9exp4f$sig.correlaton>0)]<-0
colMatlo[which(rescor.lolv4occ9exp4f$sig.correlaton<0)]<- -1
#colMatlo[which(rescor.lolv4occ9exp4f$sig.correlaton<0)]<- 0

colMatlo<-rescor.lolv4occ9exp4nosite$sig.correlaton
colMatlo[which(rescor.lolv4occ9exp4nosite$sig.correlaton>0)]<-1
colMatlo[which(rescor.lolv4occ9exp4nosite$sig.correlaton<0)]<- -1

colMatlo<-rescor.lolv4occ9exp4nosite$sig.correlaton
colMatlo[which(colMatlo>.85)]<-1
colMatlo[which(colMatlo<(-.85))]<- -1
colMatlo[which(colMatlo<.85&colMatlo>(-.85))]<-0

temp<-colMatlo[upper.tri(colMatlo)]
temp2<-temp[temp!=0]
hist(temp2)
sort(temp2)
length(temp2)

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

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/networklocircleposblue.pdf") 
plot(graphlo2,vertex.size=4,edge.curved=F,edge.color=ifelse(myedgelistlo$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphlo$color,edge.width=.7,vertex.label=NA)#,vertex.shape=shapesgraplo  
plot(graphlo2,vertex.size=4,edge.curved=F,edge.color=ifelse(myedgelistlo$weight==1,"#687dcb","#ce4d42"),vertex.color=colorgraphlo$color,edge.width=.7,vertex.label=NA)#,vertex.shape=shapesgraplo  
#dev.off()

colorgraphlo[which(colorgraphlo$group=="Mesofauna"),]
colorgraphlo[which(colorgraphlo$group2=="PhotosyntheticBacteria"),]

myedgelistlo[which(myedgelistlo$X1=="B783ef4ce2388b995de6b9b27b0c9209e"|myedgelistlo$X2=="B783ef4ce2388b995de6b9b27b0c9209e"),]

temp<-colorgraphlo[which(colorgraphlo$group=="Mesofauna"),"otu"]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp|myedgelistlo$X2%in%temp),]
dim(temp2)

#making sure it is plotting in the right order, yes. the nematode negative interactions are just covered up by all the red positive interactions in the main graph
dim(colMatlo)
rownames(colMatlo)<-1:644
colnames(colMatlo)<-1:644
graphlo1<-graph_from_adjacency_matrix(colMatlo, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
orderlo<-c(7,8,18,19,20,21,9,6,15,16,17,14,5,4,3,1,2,10,11,12,13)
colorgraphlo<-data.frame(color=as.character(c("#ff9c34",rep("gray80",20))));colorgraphlo$color<-as.character(colorgraphlo$color)

#subgraph lo
colorgraphlo2<-colorgraphlo[which(colorgraphlo$group2=="Mesofauna"),]
myedgelistlo2<-myedgelistlo[which(myedgelistlo[,"X1"]%in%colorgraphlo2$otu|myedgelistlo[,"X2"]%in%colorgraphlo2$otu),]
graph3<-subgraph.edges(graphlo2, eids=which(myedgelistlo2[,"X1"]%in%colorgraphlo2$otu|myedgelistlo2[,"X2"]%in%colorgraphlo2$otu), delete.vertices = F)
plot(graph3,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelistlo2$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphlo$color)#,rescale=F,xlim=c(-1,1),ylim=c(-1,1)





##### me #####
#creating sparse matrix
colMatme<-rescor.melv4occ9exp4f$sig.correlaton
colMatme[which(rescor.melv4occ9exp4f$sig.correlaton>0)]<-1
colMatme[which(rescor.melv4occ9exp4f$sig.correlaton<0)]<- -1

# colMatme<-rescor.melv4occ9exp4f$sig.correlaton
# colMatme[which(colMatme>.6)]<-1
# colMatme[which(colMatme<(-.6))]<- -1
# colMatme[which(colMatme<.6&colMatme>(-.6))]<-0

colMatme<-rescor.melv4occ9exp4nosite$sig.correlaton
colMatme[which(rescor.melv4occ9exp4nosite$sig.correlaton>0)]<-1
colMatme[which(rescor.melv4occ9exp4nosite$sig.correlaton<0)]<- -1

graphme1<-graph_from_adjacency_matrix(colMatme, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelistme<-data.frame(as_edgelist(graphme1),weight=E(graphme1)$weight) #just the edges

length(which(myedgelistme$weight==1))
length(which(myedgelistme$weight==-1))
length(which(myedgelistme$weight==1))/(length(which(myedgelistme$weight==1))+length(which(myedgelistme$weight==-1)))

graphme2<-graph.edgelist(as.matrix(myedgelistme[,1:2]),directed=FALSE)
graphme2

verticesgraphme<-data.frame(otu=rownames(as.matrix(V(graphme2))))
colorgraphme<-merge(verticesgraphme,labelsall,"otu",all.y=F,all.x=F,sort=F)

##use colorgraphlo$group2 for ordering, if ther are ties it leaves them in thir original order, thus is still preserves some of the ordering that makes the lines look nice
orderme<-order(colorgraphme$group2)
#orderme<-order(verticesgraphme$otu)
graphme2$layout <- layout_in_circle(graphme2,order=orderme)
#graphme2$layout <- layout_in_circle

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/networkmecircleposblue.pdf") 
plot(graphme2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelistme$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphme$color,edge.width=.7)#,layout=l3
plot(graphme2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelistme$weight==1,"#687dcb","#ce4d42"),vertex.color=colorgraphme$color,edge.width=.7)#,layout=l3
dev.off()

temp<-colorgraphme[which(colorgraphme$group=="Mesofauna"),"otu"]
temp2<-myedgelistme[which(myedgelistme$X1%in%temp|myedgelistme$X2%in%temp),]
dim(temp2)



##### hi #####
#creating sparse matrix
colMathi<-rescor.hilv4occ9exp4f$sig.correlaton
colMathi[which(rescor.hilv4occ9exp4f$sig.correlaton>0)]<-1
colMathi[which(rescor.hilv4occ9exp4f$sig.correlaton<0)]<- -1
#colMathi[which(rescor.hilv4occ9exp4f$sig.correlaton<0)]<- 0

# colMathi<-rescor.hilv4occ9exp4f$sig.correlaton
# colMathi[which(colMathi>.5)]<-1
# colMathi[which(colMathi<(-.5))]<- -1
# colMathi[which(colMathi<.5&colMathi>(-.5))]<-0

colMathi<-rescor.hilv4occ9exp4nosite$sig.correlaton
colMathi[which(rescor.hilv4occ9exp4nosite$sig.correlaton>0)]<-1
colMathi[which(rescor.hilv4occ9exp4nosite$sig.correlaton<0)]<- -1

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

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/networkhicircleposblue.pdf") 
plot(graphhi2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelisthi$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphhi$color,edge.width=.7)#,layout=l3  
plot(graphhi2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelisthi$weight==1,"#687dcb","#ce4d42"),vertex.color=colorgraphhi$color,edge.width=.7)#,layout=l3  
#dev.off()

colorgraphhi[which(colorgraphhi$group=="Mesofauna"),]

temp<-colorgraphhi[which(colorgraphhi$group=="Mesofauna"),"otu"]
temp2<-myedgelisthi[which(myedgelisthi$X1%in%temp|myedgelisthi$X2%in%temp),]
dim(temp2)

#Creating a subgraph
colorgraphhi2<-colorgraphhi[which(colorgraphhi$group2=="Mesofauna"),]
myedgelisthi2<-myedgelisthi[which(myedgelisthi[,"X1"]%in%colorgraphhi2$otu|myedgelisthi[,"X2"]%in%colorgraphhi2$otu),]
graph3<-subgraph.edges(graphhi2, eids=which(myedgelisthi2[,"X1"]%in%colorgraphhi2$otu|myedgelisthi2[,"X2"]%in%colorgraphhi2$otu), delete.vertices = F)
plot(graph3,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelisthi2$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphhi$color)#,rescale=F,xlim=c(-1,1),ylim=c(-1,1)



######random graphs#####
#creating sparse matrix

colMathi<-rescor.hirnosite$sig.correlaton
colMathi[which(rescor.hirnosite$sig.correlaton>0)]<-1
colMathi[which(rescor.hirnosite$sig.correlaton<0)]<- -1

colMathi<-rescor.lornosite$sig.correlaton
colMathi[which(rescor.lornosite$sig.correlaton>0)]<-1
colMathi[which(rescor.lornosite$sig.correlaton<0)]<- -1

colMathi<-rescor.lonosite$sig.correlaton
colMathi[which(rescor.lonosite$sig.correlaton>0)]<-1
colMathi[which(rescor.lonosite$sig.correlaton<0)]<- -1

colMathi<-rescor.hinosite$sig.correlaton
colMathi[which(rescor.hinosite$sig.correlaton>0)]<-1
colMathi[which(rescor.hinosite$sig.correlaton<0)]<- -1

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

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/networkhicircleposblue.pdf") 
#plot(graphhi2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelisthi$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphhi$color,edge.width=.7)#,layout=l3  
plot(graphhi2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelisthi$weight==1,"#687dcb","#ce4d42"),vertex.color=colorgraphhi$color,edge.width=.7)#,layout=l3  
#dev.off()

colorgraphhi[which(colorgraphhi$group=="Mesofauna"),]

temp<-colorgraphhi[which(colorgraphhi$group=="Mesofauna"),"otu"]
temp2<-myedgelisthi[which(myedgelisthi$X1%in%temp|myedgelisthi$X2%in%temp),]
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


#interactions involving mesofauna
temp<-colorgraphlo[which(colorgraphlo$group=="Mesofauna"),"otu"]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp|myedgelistlo$X2%in%temp),]
temp2
dim(temp2)
colorgraphlo[which(colorgraphlo$group=="Mesofauna"),]
colorgraphlo[which(colorgraphlo$otu%in%temp2$X1|colorgraphlo$otu%in%temp2$X2),]

temp<-colorgraphlo[which(colorgraphlo$otu=="Nb85db42310af5ddb08354eef2427cc8e"),"otu"]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp|myedgelistlo$X2%in%temp),]
colorgraphlo[which(colorgraphlo$otu%in%temp2$X2),]

Nb85db42310af5ddb08354eef2427cc8e predator nematode
N91502b60bb4cd795ba70dd05ed89e805 rotifer, one neg relationship with bacteroidetes
Nf607691c8f991641954b20294b81a9b7 rotifer, 3 neg relationship with heterotrophic bacteria, 1 neg with fungi

temp<-colorgraphme[which(colorgraphme$group=="Mesofauna"),"otu"]
temp2<-myedgelistme[which(myedgelistme$X1%in%temp|myedgelistme$X2%in%temp),]
dim(temp2)
colorgraphme[which(colorgraphme$group=="Mesofauna"),]

temp<-colorgraphhi[which(colorgraphhi$group=="Mesofauna"),"otu"]
temp2<-myedgelisthi[which(myedgelisthi$X1%in%temp|myedgelisthi$X2%in%temp),]
dim(temp2)
colorgraphhi[which(colorgraphhi$group=="Mesofauna"),]

temp<-colorgraphhi[which(colorgraphhi$otu=="N6f914ead2160e51670d3dc70c25e107b"),"otu"]
temp2<-myedgelisthi[which(myedgelisthi$X1%in%temp|myedgelisthi$X2%in%temp),]
colorgraphhi[which(colorgraphhi$otu%in%temp2$X2),]

N6f914ead2160e51670d3dc70c25e107b not well identified probably fungal feeder, Aphelenchoides fragariae, also positively correlated with 4 fungi
Nece270a38485359f4d1e2a33f964a983 Nagelus obscurus, plant parasite
     Leptosphaeria is the fungus it is positively correlated with, which is a plant pathogen
N1b65f0ff041f3f5b37e4dcb6607ea61b Pungentus sp, omnivore, positively correlated with two bacteria


#get sequences for all the mesofauna so Dorota can check them
grep -C 2 "bfe4b462f2c88679387207f098ee6220" dna-sequences.fasta


#photosynthetic microbes
temp<-colorgraphlo[which(colorgraphlo$group2=="PhotosyntheticBacteria"),"otu"]
temp<-colorgraphlo[which(colorgraphlo$group2=="PhotosyntheticEukaryota"),"otu"]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp|myedgelistlo$X2%in%temp),]
which(temp2$X1=="Ba6f5d1af8106e5c3abaf27735243cca4"|temp2$X2=="Ba6f5d1af8106e5c3abaf27735243cca4") #chlorobi (photoautotroph), 37 pos interactions with heterotrophic bacteria, 4 with unknown bacteria, 1 pos with algae, one pos photobac
which(temp2$X1=="Bc3082be027f0efa400acd22097db7422"|temp2$X2=="Bc3082be027f0efa400acd22097db7422") #Cyanobacteria, 12 interaction with het bacteria, 1 pos with photosynthetic bacteria
which(temp2$X1=="Bcbf23a259461788c7396423fe4fbf6de"|temp2$X2=="Bcbf23a259461788c7396423fe4fbf6de")#highest, chloroflexi (photoheterotroph), one pos photobac (didn't check for heterotrophs)
which(temp2$X1=="Bc4030127655425fbcde4a26e69b4859d"|temp2$X2=="Bc4030127655425fbcde4a26e69b4859d")#highest, chloroflexi (photoheterotroph), one pos photobac, one pos photoeuk (didnt check heterotrophs)
which(temp2$X1=="B3610973234fc5048dbab19dbef98d605"|temp2$X2=="B3610973234fc5048dbab19dbef98d605")#chloroflexi (photoheterotroph)
S8fcb977d1395a10dc83007a741b8a620 - Charophyta neg with fungi and bacteria
S0e23274231a1b7918be91d3630da6fe5 - Chlorophyta, 5 pos with heterotrophic bacteria, one pos with photosynthetic bactria
Sd4ea6eb20ce02a8c3e292634f92b2f10 - Chlorophyta, 2 pos with heterotrophic bacteria
Sdfc413b77592c93861526de8c906369c - chlorophyta, 2 positive with htertrophic bacteria
Sd2b04c59e082351f65cc351ef7d4ab04 - Chlorophyta 27 with het bac, 1 pos with photeuk and 1 with photobac
S6ae4196b4456d15855bbb3b088b90e1c - Archaeplastida, 12 with het bac, 1 pos with photeuk
S162acde9cd2e077ccdde140382ea48cd - Chlorophyta 2 with hetbac
Se84d5f50cb74d66781cd7b44849cccdd - chlorophyta 14 with het back
S87c0732ac6500f25336817c21ffc153e - chlorophyta 2 with het bac
Sa4906785d79edc355bd7c3048d0ed2a8 - chlorophyta 4 with het bac
Sd329e25e43835c63a4b27ccaf2403721 - Stramenopiles/Diatomea 4 pos with het bac
S9cfbc221cd3b50a2f192559d99142a13 - Stramenopiles/Diatomea 3 pos with het bac
Sc45dd62d47c0618457b55cd98c0a28c1 - Stramenopiles;__Chrysophyceae 4 pos het bac
S478bbfcf431d6dbd7438f33dc46a603e - chloroplastida (not to phylum), 1 with ht bac
37+12+5+2+2+27+12+2+14+2+4+4+3+4+1 #all
41+12+5+2+2+27+2+14+2+4+4+3+4 #this is how I got th original number I used
37+12+5+2+2+27+2+14+2+4+4+3+4 #taking out uncertain taxonomy

#2 pos b/t photobac, 2 pos b/t photobac and photo euk, 1 pos b/t photoeuk
dim(temp2)
temp2
colorgraphlo[which(colorgraphlo$otu%in%temp),]
colorgraphlo[which(colorgraphlo$otu%in%temp2$X1|colorgraphlo$otu%in%temp2$X2),]
colorgraphlo[which(colorgraphlo$otu=="S6ae4196b4456d15855bbb3b088b90e1c"),]
temp3<-myedgelistlo[which(myedgelistlo$X1=="S478bbfcf431d6dbd7438f33dc46a603e"|myedgelistlo$X2=="S478bbfcf431d6dbd7438f33dc46a603e"),]
names(temp3)[1]<-"otu"
temp3<-merge(temp3,colorgraphlo[,c(1,2)])
#temp3<-merge(temp3,colorgraphlo[,c(1,5)])
names(temp3)[1]<-"X1"
names(temp3)[2]<-"otu"
names(temp3)[4]<-"taxstring1"
temp3<-merge(temp3,colorgraphlo[,c(1,2)])
#temp3<-merge(temp3,colorgraphlo[,c(1,5)])
temp3
length(which(temp3$weight==1&temp3$group2=="HeterotrophicBacteria"|temp3$weight==1&temp3$taxstring1=="HeterotrophicBacteria"))


temp<-colorgraphme[which(colorgraphme$group2=="PhotosyntheticBacteria"),"otu"]
B783ef4ce2388b995de6b9b27b0c9209e cyanobacteria 6 het bacteria
B473e2ddb973674b822b0b8544a869f74 cyanobacteria 31 het bacteria
Bc4030127655425fbcde4a26e69b4859d chloroflexi 43 het bac, not used b/c photoheterotrohpic
Bcbf23a259461788c7396423fe4fbf6de chloroflexi 4 , not used b/c photoheterotrohpic
B23fcca0866f8adf23339be869a5c2f24 cyanobacteria 1
temp<-colorgraphme[which(colorgraphme$group2=="PhotosyntheticEukaryota"),"otu"]
S6fb7e847441f66fe0cd4eb9b74fb46eb Chlorophyta 14
Sd329e25e43835c63a4b27ccaf2403721 diatomea 9 
S1f590457d5e73c5eb225760a268c8fda chlorophyta 11
S21fd2c6d01f371f0d981893bc00ed79a chlorophyta  0
Sc14ca533c39ec57f91c0d5556164a531 chlorophyta 2
6+31+43+4+1+14+9+11+2 #all
6+31+1+14+9+11+2 #not using uncertain taxonomies

temp2<-myedgelistme[which(myedgelistme$X1%in%temp|myedgelistme$X2%in%temp),]
colorgraphme[which(colorgraphme$otu=="Sc14ca533c39ec57f91c0d5556164a531"),]
temp3<-myedgelistme[which(myedgelistme$X1=="Sc14ca533c39ec57f91c0d5556164a531"|myedgelistme$X2=="Sc14ca533c39ec57f91c0d5556164a531"),]
names(temp3)[1]<-"otu"
temp3<-merge(temp3,colorgraphme[,c(1,2)])
names(temp3)[1]<-"X1"
names(temp3)[2]<-"otu"
names(temp3)[4]<-"taxstring1"
temp3<-merge(temp3,colorgraphme[,c(1,2)])
temp3
length(which((temp3$weight==1&temp3$group2=="HeterotrophicBacteria")|(temp3$weight==1&temp3$taxstring1=="HeterotrophicBacteria")))
length(which(temp3$weight==1&temp3$group2=="HeterotrophicBacteria"))
length(which(temp3$weight==1&temp3$taxstring1=="HeterotrophicBacteria"))

temp<-colorgraphhi[which(colorgraphhi$group2=="PhotosyntheticBacteria"),"otu"]
Bc4030127655425fbcde4a26e69b4859d
temp2<-myedgelisthi[which(myedgelisthi$X1%in%temp|myedgelisthi$X2%in%temp),]
colorgraphhi[which(colorgraphhi$otu=="Bc4030127655425fbcde4a26e69b4859d"),]
temp3<-myedgelisthi[which(myedgelisthi$X1=="Bc4030127655425fbcde4a26e69b4859d"|myedgelisthi$X2=="Bc4030127655425fbcde4a26e69b4859d"),]
names(temp3)[1]<-"otu"
temp3<-merge(temp3,colorgraphhi[,c(1,2)])
names(temp3)[1]<-"X1"
names(temp3)[2]<-"otu"
names(temp3)[4]<-"taxstring1"
temp3<-merge(temp3,colorgraphhi[,c(1,2)])
temp3
length(which((temp3$weight==1&temp3$group2=="HeterotrophicBacteria")|(temp3$weight==1&temp3$taxstring1=="HeterotrophicBacteria")))


#counting total interations, not knowing who the connection is with or pos/neg
temp<-colorgraphlo[which(colorgraphlo$group2=="PhotosyntheticBacteria"|colorgraphlo$group2=="PhotosyntheticEukaryota"),"otu"]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp|myedgelistlo$X2%in%temp),]
dim(temp2) #307 interactions

temp<-colorgraphme[which(colorgraphme$group2=="PhotosyntheticBacteria"|colorgraphme$group2=="PhotosyntheticEukaryota"),"otu"]
temp2<-myedgelistme[which(myedgelistme$X1%in%temp|myedgelistme$X2%in%temp),]
temp2
dim(temp2) #153 interactions
colorgraphme[which(colorgraphme$otu%in%temp2$X1|colorgraphme$otu%in%temp2$X2),]

temp<-colorgraphhi[which(colorgraphhi$group2=="PhotosyntheticBacteria"|colorgraphhi$group2=="PhotosyntheticEukaryota"),"otu"]
temp2<-myedgelisthi[which(myedgelisthi$X1%in%temp|myedgelisthi$X2%in%temp),]
dim(temp2) #10 interactions
temp2
colorgraphhi[which(colorgraphhi$otu%in%temp2$X1|colorgraphhi$otu%in%temp2$X2),]



#plants
temp<-colorgraphlo[which(colorgraphlo$group=="Plant"),"otu"]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp|myedgelistlo$X2%in%temp),]
dim(temp2)

temp<-colorgraphme[which(colorgraphme$group=="Plant"),"otu"]
temp2<-myedgelistme[which(myedgelistme$X1%in%temp|myedgelistme$X2%in%temp),]
temp2
dim(temp2)
colorgraphme[which(colorgraphme$otu%in%temp2$X1|colorgraphme$otu%in%temp2$X2),]

temp<-colorgraphhi[which(colorgraphhi$group=="Plant"),"otu"]
temp2<-myedgelisthi[which(myedgelisthi$X1%in%temp|myedgelisthi$X2%in%temp),]
dim(temp2)
temp2
colorgraphhi[which(colorgraphhi$otu%in%temp2$X1|colorgraphhi$otu%in%temp2$X2),]



#Bacteria in networks (lo abundance, hi)
temp<-colorgraphlo[which(colorgraphlo$group=="Bacteria"),]
ind<-grep("Ktedonobacteria",temp$taxstring)
temp1<-temp[ind,]
#temp<-colorgraphlo[which(colorgraphlo$group=="Bacteria"),"otu"]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp1$otu|myedgelistlo$X2%in%temp1$otu),]
temp2
sum(temp2$weight[temp2$weight==1])
dim(temp2)
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
#subtract 26 positive interactions are between two Ktedonobacteria
dim(myedgelistlo[which(myedgelistlo$X1%in%temp1$otu&myedgelistlo$X2%in%temp1$otu),])
245-26

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
141-4

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
26-1


#high
temp<-colorgraphhi[which(colorgraphhi$group=="Fungi"),]
ind<-grep("Plenodomus_biglobosus",temp$taxstring)
temp[ind,]
Ide16ea0b9000538bf51c5d592e6003a0
temp2<-myedgelisthi[which(myedgelisthi$X1%in%c("Ide16ea0b9000538bf51c5d592e6003a0")|myedgelisthi$X2%in%c("Ide16ea0b9000538bf51c5d592e6003a0")),]
temp2




colMat<-res.corslolv4occ8exp5$sig.correlaton
colMat[which(res.corslolv4occ8exp5$sig.correlaton>0)]<-1
colMat[which(res.corslolv4occ8exp5$sig.correlaton<0)]<- -1

colMat<-res.corshilv4occ8exp4$sig.correlaton
colMat[which(res.corshilv4occ8exp4$sig.correlaton>0)]<-1
colMat[which(res.corshilv4occ8exp4$sig.correlaton<0)]<- -1

colMat<-res.corsmelv4occ8exp4$sig.correlaton
colMat[which(res.corsmelv4occ8exp4$sig.correlaton>0)]<-1
colMat[which(res.corsmelv4occ8exp4$sig.correlaton<0)]<- -1

colMat<-res.corslolv4occ8exp4$sig.correlaton
colMat[which(res.corslolv4occ8exp4$sig.correlaton>0)]<-1
colMat[which(res.corslolv4occ8exp4$sig.correlaton<0)]<- -1

colMat<-res.corshilv4occ8exp4$sig.correlaton
colMat[which(colMat>.6)]<-1
colMat[which(colMat<(-.6))]<- -1
colMat[which(colMat<.6&colMat>(-.6))]<-0

colMat<-res.corslolv4occ8exp4$sig.correlaton
colMat[which(colMat>.6)]<-1
colMat[which(colMat<(-.6))]<- -1
colMat[which(colMat<.6&colMat>(-.6))]<-0

#hist(colMat)
#dim(colMat)

#colMat2<-as(colMat, "dgCMatrix")
graph1<-graph_from_adjacency_matrix(colMat, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelist<-data.frame(as_edgelist(graph1),weight=E(graph1)$weight) #just the edges

length(which(myedgelist$weight==1))
length(which(myedgelist$weight==-1))
length(which(myedgelist$weight==1))/(length(which(myedgelist$weight==1))+length(which(myedgelist$weight==-1)))

#graph2<-simplify(graph.edgelist(as.matrix(myedgelist[,1:2]),directed=FALSE))
graph2<-graph.edgelist(as.matrix(myedgelist[,1:2]),directed=FALSE)
graph2

#(sum(colMat,na.rm=T)-322)/2 #(134 interactions)
#(24*24-24)/2 #every pairwise interaction (276)

graph2$layout <- layout_in_circle
verticesgraph<-data.frame(oldotu=rownames(as.matrix(V(graph2))))
colorgraph<-merge(verticesgraph,labelsall,"oldotu",all.y=F,all.x=F,sort=F)
plot(graph2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelist$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraph$color)#,layout=l3

unique(tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%colorgraph$oldotu])
temp<-myedgelist[which(myedgelist$X1%in%c("8f0fb68673a8b3c26f3fd210e29b7916")|myedgelist$X2%in%c("8f0fb68673a8b3c26f3fd210e29b7916")),]
8f0fb68673a8b3c26f3fd210e29b7916 #clostridium bowmanii, byproducts of acetic and butyric acid
temp2<-temp[which(temp$weight==1),]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]

temp<-myedgelist[which(myedgelist$X1%in%c("6472eb8b1e09f892aca2f23182962903")|myedgelist$X2%in%c("6472eb8b1e09f892aca2f23182962903")),]
8f0fb68673a8b3c26f3fd210e29b7916  #bradyrhizobium fixes N
temp2<-temp[which(temp$weight==1),]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datEukN3)[rownames(tax_table(datEukN3))%in%c(as.character(temp2$X1),as.character(temp2$X2))] #one chlorophyceae - possible syntrophic relationship
tax_table(datITSS3)[rownames(tax_table(datITSS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]


"631c2b290edb6b71d9e6b5577a0613ae"# syntrophobacteriaceae
"8a8895d160b2fda424ad5446e9bf8f4e"#Desulfosporosinus sulfate reducing
temp<-myedgelist[which(myedgelist$X1%in%c("8a8895d160b2fda424ad5446e9bf8f4e")|myedgelist$X2%in%c("8a8895d160b2fda424ad5446e9bf8f4e")),]
temp2<-temp[which(temp$weight==1),]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datITSS3)[rownames(tax_table(datITSS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]

"b8f64099026ffffde7cd2c6aa4642240" #nitrosomonadaceae, ammonia oxidizers
temp<-myedgelist[which(myedgelist$X1%in%c("b8f64099026ffffde7cd2c6aa4642240")|myedgelist$X2%in%c("b8f64099026ffffde7cd2c6aa4642240")),]
temp2<-temp[which(temp$weight==1),]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datITSS3)[rownames(tax_table(datITSS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]

myedgelist[myedgelist$X1%in%temp,]
myedgelist[myedgelist$X2%in%temp,]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(temp)]
colorgraph[which(colorgraph$oldotu=="854124b2d8e0badeb138eb9e0e1808c2"),]

#Creating a subgraph
colorgraph3<-colorgraph[which(colorgraph$labels=="Bacteria"),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Bacteria"),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Mesofauna"),]
tax_table(datEukN)[rownames(tax_table(datEukN))%in%colorgraph3$oldotu]

myedgelist2<-myedgelist[which(myedgelist[,"X1"]%in%colorgraph3$oldotu|myedgelist[,"X2"]%in%colorgraph3$oldotu),]
graph3<-subgraph.edges(graph2, eids=which(myedgelist[,"X1"]%in%colorgraph3$oldotu|myedgelist[,"X2"]%in%colorgraph3$oldotu), delete.vertices = F)
plot(graph3,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelist2$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraph$color)#,rescale=F,xlim=c(-1,1),ylim=c(-1,1)

#investigating nematode relationships
#lo
temp<-"b85db42310af5ddb08354eef2427cc8e" #omnivore nematode in lo
temp1<-myedgelist[which(myedgelist$X1=="b85db42310af5ddb08354eef2427cc8e"|myedgelist$X2=="b85db42310af5ddb08354eef2427cc8e"),]
temp2<-temp1[temp1$weight==-1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#5 bacteria
tax_table(datEukS)[rownames(tax_table(datEukS))%in%temp2$X2]#2 euks
tax_table(datITSS)[rownames(tax_table(datITSS))%in%temp2$X2]#0
temp2<-temp1[temp1$weight==1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#6
tax_table(datEukS)[rownames(tax_table(datEukS))%in%temp2$X2]#1
tax_table(datITSS)[rownames(tax_table(datITSS))%in%temp2$X2]#6
tax_table(datEukN)[rownames(tax_table(datEukN))%in%temp2$X2]#0

#hi
6f914ead2160e51670d3dc70c25e107b __Aphelenchida fungal feeder
e485c9c4bdda247bc813c36275d99dce __Plectidae bacterial feeder
temp<-"6f914ead2160e51670d3dc70c25e107b" #ff nematode in lo
temp1<-myedgelist[which(myedgelist$X1=="6f914ead2160e51670d3dc70c25e107b"|myedgelist$X2=="6f914ead2160e51670d3dc70c25e107b"),]
#only positive rels
temp2<-temp1[temp1$weight==1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#2 bacteria
temp<-"e485c9c4bdda247bc813c36275d99dce" #bf nematode in lo
temp1<-myedgelist[which(myedgelist$X1=="e485c9c4bdda247bc813c36275d99dce"|myedgelist$X2=="e485c9c4bdda247bc813c36275d99dce"),]
temp2<-temp1[temp1$weight==-1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#3 bacteria
tax_table(datEukS)[rownames(tax_table(datEukS))%in%temp2$X2]#0 euks
tax_table(datITSS)[rownames(tax_table(datITSS))%in%temp2$X2]#0
temp2<-temp1[temp1$weight==1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#3 bacteria
tax_table(datEukS)[rownames(tax_table(datEukS))%in%temp2$X2]#
tax_table(datITSS)[rownames(tax_table(datITSS))%in%temp2$X2]
# 2 plants 





#trying to keep the layout the same and plot only partial graph
l3<-layout_in_circle(graph2)#layout_with_fr(graph2) #I don't want to run this accidentally
rownames(l3) <- V(graph2)$name
l3b<-layout.norm(l3,xmin=-1,xmax=1,ymin=-1,ymax=1)
#plot(graph3,vertex.size=sizesgraph3,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA,vertex.shape=shapesgraph3,layout=l3b,rescale=F,xlim=c(-1,1),ylim=c(-1,1))




##Exploring nematode correlations in hi
#use colorgraph from below
colorgraph2<-colorgraph[which(colorgraph$labels=="Mesofauna"|colorgraph$labels=="Eukaryota"),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Mesofauna"|colorgraph$labels=="Eukaryota"),]
colorgraph2<-colorgraph[which(colorgraph$labels=="Mesofauna"),]
colorgraph2
tax_table(datEukN3)[rownames(tax_table(datEukN3))%in%colorgraph2$oldotu]
tax_table(datEukN)[rownames(tax_table(datEukN))%in%colorgraph2$oldotu]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%colorgraph2$oldotu]
temp<-c("a64f0e90f25d607e4c34fffa870f19e4","71f3f9709801ee03a20e33c46a8b797d","a7a7353d6231f298e75b69f4161926a1","ee45b87e6ac04257510a6d99644760cd","eafa1d1d561a47042f099e619cd25693","e4ae90d4c05ef00793fc6b93fb6a9af7","f7b389ab203fda4cc3cf1fffe23004a6","9bd3df7f714986d788a2cce65be56d10","7a97813269020725286a63989363ceef","d3898987a65cc50f9226529374935e74","854124b2d8e0badeb138eb9e0e1808c2")#lo 10
temp<-c("f4e7445ced25e7028c593551658c091a","f8d2629f619aea2ec1a4a2571c91f1af","01a1a6cfe5266fe585bed622ab7f64c4","124f8bbd2e10abcb816e84a974280321","3094f58213b725653eb04f35aea790f8","479660228dd139e20ad8ff42698a123e","65cb8aa63bed1f44c9a576ae2e94e0ab")#hi 8
temp<-c("f8d2629f619aea2ec1a4a2571c91f1af","01a1a6cfe5266fe585bed622ab7f64c4","124f8bbd2e10abcb816e84a974280321","3094f58213b725653eb04f35aea790f8")#hi 10
colorgraph3<-colorgraph[which(colorgraph$labels=="Mesofauna"|colorgraph$oldotu%in%temp),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Fungi"),]

#tax_table(datEukN3)[rownames(tax_table(datEukN3))=="6f914ead2160e51670d3dc70c25e107b"]

myedgelist2<-myedgelist[which(myedgelist[,"X1"]%in%colorgraph3$oldotu|myedgelist[,"X2"]%in%colorgraph3$oldotu),]
graph3<-graph.edgelist(as.matrix(myedgelist2[,1:2]),directed=FALSE)
graph3
graph3$layout <- layout_in_circle
verticesgraph<-data.frame(oldotu=rownames(as.matrix(V(graph3))))
colorgraph4<-merge(verticesgraph,labelsall,"oldotu",all.y=F,all.x=F,sort=F)
plot(graph3,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelist2$weight==1,"red","blue"),vertex.color=colorgraph$color,layout=l3b)#,rescale=F,xlim=c(-1,1),ylim=c(-1,1)



#using hi occ8 with correlation >.5 or <-.5, I get 10 mesofauna taxa
#using hi occ10 significnat correlations, I get 5 mesofauna and 4 heterotrophic eukaryotes. 67/51
#using hi occ8 significnat correlations, I get 5 mesofauna and 7 heterotrophic eukaryotes. 53/45
#lo 10, significant, 1 mesofuana, 10 heterotrophic euks 128/81
#using hi occ8 with 90%CI, 7 mesofauna (didn't look for euks)









##### Chord diagrams #####
(length(which(res.corslo$sig.correlaton!=0))-dim(res.corslo$sig.correlaton)[1])/2
colMat <- matrix(NA, nrow = nrow(rescor.lolv4occ10exp4$sig.correlaton), ncol = ncol(rescor.lolv4occ10exp4$sig.correlaton))
colMat[which(res.cors$sig.correlaton > 0.63, arr.ind = TRUE)] <- "red"
colMat[which(res.cors$sig.correlaton < -0.63, arr.ind = TRUE)] <- "blue"
chordDiagram(res.cors$sig.correlaton, symmetric = TRUE,annotationTrack = c("name", "grid"),grid.col = "grey",col=colMat)

colMat <- matrix(NA, nrow = nrow(res.corshiocc8$correlation), ncol = ncol(res.corshiocc8$correlation),dimnames = list(rownames(res.corshiocc8$correlation),rownames(res.corshiocc8$correlation)))
colMat[which(res.corshiocc8$correlation > 0.5, arr.ind = TRUE)] <- 1
colMat[which(res.corshiocc8$correlation < -0.5, arr.ind = TRUE)] <- -1
colMat[which(res.corshiocc8$correlation<.5&res.corshiocc8$correlation>(-.5))]<-0

length(which(res.cor90$correlation>.7&res.cor90$correlation<1,arr.ind = T))
length(which(res.corshi8$correlation>.7&res.corshi8$correlation<1,arr.ind = T))
hist(res.cors$correlation)
hist(res.corshi8$correlation)





#Trial runs were stored in MovingUphill3_Workspace_Analysis3b.Rdata
#fit.lvmhiocc6 
#load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill3_Workspace_Analysis3b.Rdata")
#fit.hilv4occ8exp5
#fit.lolv4occ8exp5
#fit.hilv4occ8exp4 (and b)
#fit.melv4occ8exp4

#res.cors (hi with 4latent)
#res.corslo (lo with 4 latent)
#res.corshi8 (hi with 8 latent) - hardly any correlations significant or strong
#res.corshiocc8 (hi with 8 occurrences)
#res.corshiocc8.90 (hi with 8 occurrences and 90CI)
#res.corshiocc6 (hi 6 occ)
#res.corshilv4occ8exp5 (hi with 5 environmetnal variables)
#res.corslolv4occ8exp5



##### Checking the models with multiple regressions ####
#there is no negative binomial option in glm, so using quasipoisson
#overall, I think I'm ok with th number of nonzero points (9 or greater) and having 4 explanatory variables
hmscXe
hmscYe

#n=10
sum(hmscYe[,17]>0)
m1<-glm(hmscYe[,17]~hmscXe[,2:5],family=quasipoisson)
m1<-glm(hmscYe[,17]~hmscXe[,c(2,3,4)],family=quasipoisson)
m1<-glm(hmscYe[,17]~hmscXe[,c(2,4,5)],family=quasipoisson)
summary(m1)
m1fitted<-predict(m1,type="response")
plot(hmscXe[,3],hmscYe[,17])
points(hmscXe[,3],m1fitted,col=2)
#

###### Checking whether I should include a random plot effect #####
#no, because a random plot effect standardizes the whole dataset



##### Functions #####
#90% CI
get.residual.cor2<-function (object, est = "median", prob = .9) 
{
  if (is.null(object$jags.model)) 
    stop("MCMC samples not found")
  fit.mcmc <- get.mcmcsamples(object)
  y <- object$y
  X <- object$X
  num.lv <- object$num.lv
  if (length(grep("lvs", colnames(fit.mcmc))) == 0) 
    stop("Cannot find MCMC samples corresponding to latent variables")
  n <- nrow(y)
  p <- ncol(y)
  sig_rescor_mat <- rescor_mat <- rescov_mat <- matrix(0, nrow = p, 
                                                       ncol = p)
  sig_respres_mat <- respres_mat <- matrix(0, nrow = p, ncol = p)
  if (is.null(colnames(y))) 
    colnames(y) <- 1:ncol(y)
  rownames(rescor_mat) <- colnames(rescor_mat) <- rownames(sig_rescor_mat) <- colnames(sig_rescor_mat) <- colnames(y)
  rownames(rescov_mat) <- colnames(rescov_mat) <- colnames(y)
  rownames(respres_mat) <- colnames(respres_mat) <- rownames(sig_respres_mat) <- colnames(sig_respres_mat) <- colnames(y)
  all_rescor_mat <- all.rescov_mat <- all.respres_mat <- array(0, 
                                                               dim = c(nrow(fit.mcmc), p, p))
  all_trace_rescor <- numeric(nrow(fit.mcmc))
  for (k0 in 1:nrow(fit.mcmc)) {
    lv.coefs <- matrix(fit.mcmc[k0, grep("lv.coefs", colnames(fit.mcmc))], 
                       nrow = p)
    lambdalambdaT <- tcrossprod(as.matrix(lv.coefs[, 2:(num.lv + 
                                                          1)]))
    all.rescov_mat[k0, , ] <- (lambdalambdaT)
    all_trace_rescor[k0] <- sum(diag(lambdalambdaT))
    all_rescor_mat[k0, , ] <- cov2cor(all.rescov_mat[k0, 
                                                     , ])
    all.respres_mat[k0, , ] <- ginv(all_rescor_mat[k0, , 
                                                   ])
  }
  for (j in 1:p) {
    for (j2 in 1:p) {
      if (est == "median") {
        rescov_mat[j, j2] <- median(all.rescov_mat[, 
                                                   j, j2])
        rescor_mat[j, j2] <- median(all_rescor_mat[, 
                                                   j, j2])
        respres_mat[j, j2] <- median(all.respres_mat[, 
                                                     j, j2])
      }
      if (est == "mean") {
        rescov_mat[j, j2] <- mean(all.rescov_mat[, j, 
                                                 j2])
        rescor_mat[j, j2] <- mean(all_rescor_mat[, j, 
                                                 j2])
        respres_mat[j, j2] <- mean(all.respres_mat[, 
                                                   j, j2])
      }
      sig_rescor_mat[j, j2] <- rescor_mat[j, j2]
      get.hpd.cors <- HPDinterval(as.mcmc(all_rescor_mat[, 
                                                         j, j2]), prob = .9)
      if (0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) 
        sig_rescor_mat[j, j2] <- 0
      sig_respres_mat[j, j2] <- respres_mat[j, j2]
      get.hpd.cors <- HPDinterval(as.mcmc(all.respres_mat[, 
                                                          j, j2]), prob = .9)
      if (0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) 
        sig_respres_mat[j, j2] <- 0
    }
  }
  if (est == "median") 
    final_trace <- median(all_trace_rescor)
  if (est == "mean") 
    final_trace <- mean(all_trace_rescor)
  return(list(correlation = rescor_mat, sig.correlaton = sig_rescor_mat, 
              covariance = rescov_mat, precision = respres_mat, sig.precision = sig_respres_mat, 
              trace = final_trace))
}
