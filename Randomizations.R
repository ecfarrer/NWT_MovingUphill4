#Randomizations

#Randomize my data, create networks to see false positive rate

#files needed:
comm.bio

#### lo ####
ind<-which(comm.bio$lomehi=="lo")
hmscYlo<-comm.bio[ind,55:24565]  #start at 54 if you want lomehi column
rownames(hmscYlo)<-comm.bio$X.SampleID[ind]

hmscXlo<-data.frame(snowdepth=comm.bio$snowdepth,TC=comm.bio$TC,pH=comm.bio$pH,moisture=comm.bio$moisture)[ind,] 
rownames(hmscXlo)<-comm.bio$X.SampleID[ind]

#take out species that are zeros
ind<-which(colSums(hmscYlo)>0);length(ind)
hmscYlo2<-hmscYlo[ind]

#separate N, S, 16S, ITS
hmscYlo2N<-hmscYlo2[,1:101]
hmscYlo2S<-hmscYlo2[,102:1247]
hmscYlo2Bac<-hmscYlo2[,1248:7101]
hmscYlo2ITS<-hmscYlo2[,7102:8426]
hmscYlo2Plant<-hmscYlo2[,8427:8442]

#randomize, this preserves species occurrence frequency and sample species richness
outputlo<-data.frame(rep=rep(NA,10),pos=rep(NA,10),neg=rep(NA,10),tot=rep(NA,10))
outputlomod<-list(NA)
outputlorescor<-list(NA)

for (i in 4:10){
hmscYlo2Nr<-randomizeMatrix(hmscYlo2N, null.model="independentswap",iterations=30000)#
hmscYlo2Sr<-randomizeMatrix(hmscYlo2S, null.model="independentswap",iterations=500000)#500000
hmscYlo2Bacr<-randomizeMatrix(hmscYlo2Bac, null.model="independentswap",iterations=5000000)#5000000
hmscYlo2ITSr<-randomizeMatrix(hmscYlo2ITS, null.model="independentswap",iterations=500000)
hmscYlo2Plantr<-randomizeMatrix(hmscYlo2Plant, null.model="independentswap",iterations=30000)

#this gives you a rough idea of how well it is randomizing - it tells you how many columns were not randomized at all. ("successful" randomization could still just be on number different, but it tells you something)
sum(colSums(hmscYlo2ITS[,]==hmscYlo2ITSr[,])==25)

#replace so I can carry on with same code
hmscYlo2N<-hmscYlo2Nr
hmscYlo2S<-hmscYlo2Sr
hmscYlo2Bac<-hmscYlo2Bacr
hmscYlo2ITS<-hmscYlo2ITSr
hmscYlo2Plant<-hmscYlo2Plantr
hmscYlo2<-cbind(hmscYlo2N,hmscYlo2S,hmscYlo2Bac,hmscYlo2ITS,hmscYlo2Plant)

#Impute zeros
hmscYlo2N2 <- cmultRepl(hmscYlo2N,label=0, method="CZM")
hmscYlo2N2[hmscYlo2N2<=0]<-min(hmscYlo2N2[hmscYlo2N2>0])
hmscYlo2S2 <- cmultRepl(hmscYlo2S,label=0, method="CZM")
hmscYlo2S2[hmscYlo2S2<=0]<-min(hmscYlo2S2[hmscYlo2S2>0])
hmscYlo2Bac2 <- cmultRepl(hmscYlo2Bac,label=0, method="CZM")
hmscYlo2Bac2[hmscYlo2Bac2<=0]<-min(hmscYlo2Bac2[hmscYlo2Bac2>0])
hmscYlo2ITS2 <- cmultRepl(hmscYlo2ITS,label=0, method="CZM")
hmscYlo2ITS2[hmscYlo2ITS2<=0]<-min(hmscYlo2ITS2[hmscYlo2ITS2>0])

#Calculate clr
hmscYlo2N3 <- t(apply(hmscYlo2N2, 1, function(x){log(x) - mean(log(x))}))
hmscYlo2S3 <- t(apply(hmscYlo2S2, 1, function(x){log(x) - mean(log(x))}))
hmscYlo2Bac3 <- t(apply(hmscYlo2Bac2, 1, function(x){log(x) - mean(log(x))}))
hmscYlo2ITS3 <- t(apply(hmscYlo2ITS2, 1, function(x){log(x) - mean(log(x))}))
hmscYlo3<-cbind(hmscYlo2N3,hmscYlo2S3,hmscYlo2Bac3,hmscYlo2ITS3,hmscYlo2Plant)

ind<-which(colSums(hmscYlo2>0)>11)
hmscYlo4<-hmscYlo3[,ind]
hmscXlo2<-scale(hmscXlo)
hmscYlo5<-as.matrix(hmscYlo4)
hmscXlo3<-as.matrix(hmscXlo2)

#run boral model
mod.lo11flv3r<- boral(y = hmscYlo5, X = hmscXlo3, lv.control = list(num.lv = 3), family = c(rep("normal",305),rep("negative.binomial",1)), save.model = T, calc.ics = F, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.lo11flv3r <- get.residual.cor(mod.lo11flv3r) 

outputlomod[[i]]<-mod.lo11flv3r #this saves the data too $X and $y
outputlorescor[[i]]<-rescor.lo11flv3r

colMatlo<-rescor.lo11flv3r$sig.correlaton
colMatlo[which(rescor.lo11flv3r$sig.correlaton>0)]<-1
colMatlo[which(rescor.lo11flv3r$sig.correlaton<0)]<- -1

graphlo1<-graph_from_adjacency_matrix(colMatlo, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelistlo<-data.frame(as_edgelist(graphlo1),weight=E(graphlo1)$weight) #just the edges

outputlo[i,"rep"]<-i
outputlo[i,"pos"]<-length(which(myedgelistlo$weight==1))
outputlo[i,"neg"]<-length(which(myedgelistlo$weight==-1))
outputlo[i,"tot"]<-length(which(myedgelistlo$weight==1))+length(which(myedgelistlo$weight==-1))
}

#to check that there were mcmc samples
plot(get.mcmcsamples(outputlomod[[3]])[,2])

#look at the network complexity in one of the randomizations
i=7
rescor.lo11flv3r<-outputlorescor[[i]]
colMatlo<-rescor.lo11flv3r$sig.correlaton
colMatlo[which(rescor.lo11flv3r$sig.correlaton>0)]<-1
colMatlo[which(rescor.lo11flv3r$sig.correlaton<0)]<- -1

graphlo1<-graph_from_adjacency_matrix(colMatlo, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelistlo<-data.frame(as_edgelist(graphlo1),weight=E(graphlo1)$weight) #just the edges

graphlo2<-graph.edgelist(as.matrix(myedgelistlo[,1:2]),directed=FALSE)
length(E(graphlo2))/length(V(graphlo2))



## medium ##
ind<-which(comm.bio$lomehi=="me")

hmscYme<-comm.bio[ind,55:24565]  #start at 54 if you want lomehi column
rownames(hmscYme)<-comm.bio$X.SampleID[ind]

hmscXme<-data.frame(snowdepth=comm.bio$snowdepth,TC=comm.bio$TC,pH=comm.bio$pH,moisture=comm.bio$moisture)[ind,] 
rownames(hmscXme)<-comm.bio$X.SampleID[ind]

#take out species that are zeros
ind<-which(colSums(hmscYme)>0);length(ind)
hmscYme2<-hmscYme[ind]

#separate N, S, 16S, ITS
hmscYme2N<-hmscYme2[,1:167]
hmscYme2S<-hmscYme2[,168:1465]
hmscYme2Bac<-hmscYme2[,1466:9375]
hmscYme2ITS<-hmscYme2[,9376:10937]
hmscYme2Plant<-hmscYme2[,10938:10971]

#randomize, this preserves species occurrence frequency and sample species richness
outputme<-data.frame(rep=rep(NA,10),pos=rep(NA,10),neg=rep(NA,10),tot=rep(NA,10))
outputmemod<-list(NA)
outputmerescor<-list(NA)

for (j in 2:10){
hmscYme2Nr<-randomizeMatrix(hmscYme2N, null.model="independentswap",iterations=30000)
hmscYme2Sr<-randomizeMatrix(hmscYme2S, null.model="independentswap",iterations=500000)
hmscYme2Bacr<-randomizeMatrix(hmscYme2Bac, null.model="independentswap",iterations=5000000)
hmscYme2ITSr<-randomizeMatrix(hmscYme2ITS, null.model="independentswap",iterations=500000)
hmscYme2Plantr<-randomizeMatrix(hmscYme2Plant, null.model="independentswap",iterations=30000)

sum(colSums(hmscYme2Plant[,]==hmscYme2Plantr[,])==25)

#replace so I can carry on with same code
hmscYme2N<-hmscYme2Nr
hmscYme2S<-hmscYme2Sr
hmscYme2Bac<-hmscYme2Bacr
hmscYme2ITS<-hmscYme2ITSr
hmscYme2Plant<-hmscYme2Plantr
hmscYme2<-cbind(hmscYme2N,hmscYme2S,hmscYme2Bac,hmscYme2ITS,hmscYme2Plant)

#Impute zeros
hmscYme2N2 <- cmultRepl(hmscYme2N,label=0, method="CZM")
hmscYme2N2[hmscYme2N2<=0]<-min(hmscYme2N2[hmscYme2N2>0])
hmscYme2S2 <- cmultRepl(hmscYme2S,label=0, method="CZM");min(hmscYme2S2)#
hmscYme2S2[hmscYme2S2<=0]<-min(hmscYme2S2[hmscYme2S2>0])
hmscYme2Bac2 <- cmultRepl(hmscYme2Bac,label=0, method="CZM")
hmscYme2Bac2[hmscYme2Bac2<=0]<-min(hmscYme2Bac2[hmscYme2Bac2>0])
hmscYme2ITS2 <- cmultRepl(hmscYme2ITS,label=0, method="CZM")
hmscYme2ITS2[hmscYme2ITS2<=0]<-min(hmscYme2ITS2[hmscYme2ITS2>0])

#Calculate clr
hmscYme2N3 <- t(apply(hmscYme2N2, 1, function(x){log(x) - mean(log(x))}))
hmscYme2S3 <- t(apply(hmscYme2S2, 1, function(x){log(x) - mean(log(x))}))
hmscYme2Bac3 <- t(apply(hmscYme2Bac2, 1, function(x){log(x) - mean(log(x))}))
hmscYme2ITS3 <- t(apply(hmscYme2ITS2, 1, function(x){log(x) - mean(log(x))}))
hmscYme3<-cbind(hmscYme2N3,hmscYme2S3,hmscYme2Bac3,hmscYme2ITS3,hmscYme2Plant)

ind<-which(colSums(hmscYme2>0)>11)
hmscYme4<-hmscYme3[,ind]
hmscXme2<-scale(hmscXme)
hmscYme5<-as.matrix(hmscYme4)
hmscXme3<-as.matrix(hmscXme2)

#run boral model
mod.me11flv3r<- boral(y = hmscYme5, X = hmscXme3, lv.control = list(num.lv = 3), family = c(rep("normal",298),rep("negative.binomial",3)), save.model = T, calc.ics = F, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.me11flv3r <- get.residual.cor(mod.me11flv3r) 

outputmemod[[j]]<-mod.me11flv3r #this saves the data too $X and $y
outputmerescor[[j]]<-rescor.me11flv3r

colMatme<-rescor.me11flv3r$sig.correlaton
colMatme[which(rescor.me11flv3r$sig.correlaton>0)]<-1
colMatme[which(rescor.me11flv3r$sig.correlaton<0)]<- -1

graphme1<-graph_from_adjacency_matrix(colMatme, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelistme<-data.frame(as_edgelist(graphme1),weight=E(graphme1)$weight) #just the edges

outputme[j,"rep"]<-j
outputme[j,"pos"]<-length(which(myedgelistme$weight==1))
outputme[j,"neg"]<-length(which(myedgelistme$weight==-1))
outputme[j,"tot"]<-length(which(myedgelistme$weight==1))+length(which(myedgelistme$weight==-1))
}




## high ##
ind<-which(comm.bio$lomehi=="hi")

hmscYhi<-comm.bio[ind,55:24565]  #start at 54 if you want lomehi column
rownames(hmscYhi)<-comm.bio$X.SampleID[ind]

hmscXhi<-data.frame(snowdepth=comm.bio$snowdepth,TC=comm.bio$TC,pH=comm.bio$pH,moisture=comm.bio$moisture)[ind,] #
rownames(hmscXhi)<-comm.bio$X.SampleID[ind]

#take out species that are zeros
ind<-which(colSums(hmscYhi)>0);length(ind)
hmscYhi2<-hmscYhi[ind]

#separate N, S, 16S, ITS
hmscYhi2N<-hmscYhi2[,1:282]
hmscYhi2S<-hmscYhi2[,283:2070]
hmscYhi2Bac<-hmscYhi2[,2071:9952]
hmscYhi2ITS<-hmscYhi2[,9953:11893]
hmscYhi2Plant<-hmscYhi2[,11894:11942]

#randomize, this preserves species occurrence frequency and sample species richness
outputhi<-data.frame(rep=rep(NA,10),pos=rep(NA,10),neg=rep(NA,10),tot=rep(NA,10))
outputhimod<-list(NA)
outputhirescor<-list(NA)

for (k in 3:10){
hmscYhi2Nr<-randomizeMatrix(hmscYhi2N, null.model="independentswap",iterations=30000)
hmscYhi2Sr<-randomizeMatrix(hmscYhi2S, null.model="independentswap",iterations=800000)
hmscYhi2Bacr<-randomizeMatrix(hmscYhi2Bac, null.model="independentswap",iterations=6000000)
hmscYhi2ITSr<-randomizeMatrix(hmscYhi2ITS, null.model="independentswap",iterations=800000)
hmscYhi2Plantr<-randomizeMatrix(hmscYhi2Plant, null.model="independentswap",iterations=30000)

#sum(colSums(hmscYhi2Plant[,]==hmscYhi2Plantr[,])==25)

#replace so I can carry on with same code
hmscYhi2N<-hmscYhi2Nr
hmscYhi2S<-hmscYhi2Sr
hmscYhi2Bac<-hmscYhi2Bacr
hmscYhi2ITS<-hmscYhi2ITSr
hmscYhi2Plant<-hmscYhi2Plantr
hmscYhi2<-cbind(hmscYhi2N,hmscYhi2S,hmscYhi2Bac,hmscYhi2ITS,hmscYhi2Plant)

#Impute zeros
hmscYhi2N2 <- cmultRepl(hmscYhi2N,label=0, method="CZM")
hmscYhi2N2[hmscYhi2N2<=0]<-min(hmscYhi2N2[hmscYhi2N2>0])
hmscYhi2S2 <- cmultRepl(hmscYhi2S,label=0, method="CZM");min(hmscYhi2S2)
hmscYhi2S2[hmscYhi2S2<=0]<-min(hmscYhi2S2[hmscYhi2S2>0])
hmscYhi2Bac2 <- cmultRepl(hmscYhi2Bac,label=0, method="CZM")
hmscYhi2Bac2[hmscYhi2Bac2<=0]<-min(hmscYhi2Bac2[hmscYhi2Bac2>0])
hmscYhi2ITS2 <- cmultRepl(hmscYhi2ITS,label=0, method="CZM");min(hmscYhi2ITS2)
hmscYhi2ITS2[hmscYhi2ITS2<=0]<-min(hmscYhi2ITS2[hmscYhi2ITS2>0])

#Calculate clr
hmscYhi2N3 <- t(apply(hmscYhi2N2, 1, function(x){log(x) - mean(log(x))}))
hmscYhi2S3 <- t(apply(hmscYhi2S2, 1, function(x){log(x) - mean(log(x))}))
hmscYhi2Bac3 <- t(apply(hmscYhi2Bac2, 1, function(x){log(x) - mean(log(x))}))
hmscYhi2ITS3 <- t(apply(hmscYhi2ITS2, 1, function(x){log(x) - mean(log(x))}))
hmscYhi3<-cbind(hmscYhi2N3,hmscYhi2S3,hmscYhi2Bac3,hmscYhi2ITS3,hmscYhi2Plant)

ind<-which(colSums(hmscYhi2>0)>11)
hmscYhi4<-hmscYhi3[,ind]
hmscXhi2<-scale(hmscXhi)
hmscYhi5<-as.matrix(hmscYhi4)
hmscXhi3<-as.matrix(hmscXhi2)

#run boral model
mod.hi11flv3r<- boral(y = hmscYhi5, X = hmscXhi3, lv.control = list(num.lv = 3), family = c(rep("normal",265),rep("negative.binomial",8)), save.model = T, calc.ics = F, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.hi11flv3r <- get.residual.cor(mod.hi11flv3r) 

outputhimod[[k]]<-mod.hi11flv3r #this saves the data too $X and $y
outputhirescor[[k]]<-rescor.hi11flv3r

colMathi<-rescor.hi11flv3r$sig.correlaton
colMathi[which(rescor.hi11flv3r$sig.correlaton>0)]<-1
colMathi[which(rescor.hi11flv3r$sig.correlaton<0)]<- -1

graphhi1<-graph_from_adjacency_matrix(colMathi, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelisthi<-data.frame(as_edgelist(graphhi1),weight=E(graphhi1)$weight) #just the edges

outputhi[k,"rep"]<-k
outputhi[k,"pos"]<-length(which(myedgelisthi$weight==1))
outputhi[k,"neg"]<-length(which(myedgelisthi$weight==-1))
outputhi[k,"tot"]<-length(which(myedgelisthi$weight==1))+length(which(myedgelisthi$weight==-1))
}

save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/MovingUphill4_WorkspaceRandomizations.Rdata")  # 


mean(outputlo[,"tot"]) #2.2
std.error(outputlo[,"tot"])# 1.7
mean(outputme[,"tot"]) #4.5
std.error(outputme[,"tot"]) #1.9
mean(outputhi[,"tot"]) #38.6
std.error(outputhi[,"tot"]) #11.6

alldat<-rbind(outputlo,outputme,outputhi)
alldat$Succession<-rep(c("Early","Mid","Late"),each=10)

alldat%>%
  group_by(Succession)%>%
  summarise(mean=mean(tot),se=std.error(tot))

alldat2<-alldat%>%
  mutate(tot=100*tot/rep(c(2829,2037,594),each=10))%>%
  group_by(Succession)%>%
  summarise(mean=mean(tot),se=std.error(tot))
alldat2$Succession<-factor(alldat2$Succession,levels=c("Early","Mid","Late"))

ggplot(alldat2,aes(x=Succession,y=mean))+
  labs(x = "",y="% False Positives")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.4)+#.5
  geom_point(size=1.5)+#2
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.4)+#.5
  guides(col = guide_legend(ncol = 1))



#I'm not sure why the high density plots have more false positives. there were 5 bacterialtaxa that seemed not to have been randomized but it was mainly eukaroytes that made up the networks and there weren't any euks that didn't get ranodmized . Plus, the high density data have fwere species going into the networks, and they have higher diversity so more species to mess up any negative correlations du to relative nature of the data. Possibly has to do with different structure of the environment??


##### Investigating models ####

#which bacteria got randomized
temp <- cmultRepl(hmscYhi2Bac,label=0, method="CZM")
temp2 <- t(apply(temp, 1, function(x){log(x) - mean(log(x))}))
temp3<-temp2[,which(colSums(hmscYhi2Bac>0)>11)]

temp<-outputhimod[[9]]$y[,39:247];colnames(temp)
#sum(colSums(temp[,]==temp3[,])==25)#doen't work b/c the zeros are imputed as a random process,so I just have to look through them or do some comparison with negative vs positive entries.
sum(colSums(sign(temp[,])==sign(temp3[,]))==25) # there are 5 whose signs are the same suggesting it is not or not well randomized
d=125;cbind(temp[,d],temp3[,d]);plot(temp[,d],temp3[,d])

#which euks got randomized
temp <- cmultRepl(hmscYhi2S,label=0, method="CZM")
temp2 <- t(apply(temp, 1, function(x){log(x) - mean(log(x))}))
temp3<-temp2[,which(colSums(hmscYhi2S>0)>11)]

temp<-outputhimod[[9]]$y[,10:38];colnames(temp)
sum(colSums(sign(temp[,])==sign(temp3[,]))==25) # there are 0 whose signs are the same suggesting it is not or not well randomized
d=15;cbind(temp[,d],temp3[,d]);plot(temp[,d],temp3[,d])


colMathi<-outputhirescor[[10]]$sig.correlaton
colMathi[which(outputhirescor[[10]]$sig.correlaton>0)]<-1
colMathi[which(outputhirescor[[10]]$sig.correlaton<0)]<- -1

graphhi1<-graph_from_adjacency_matrix(colMathi, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelisthi<-data.frame(as_edgelist(graphhi1),weight=E(graphhi1)$weight) #just the edges



#Code from when I deleted some of my initial run models by accident, then reimported them
rm(list=setdiff(ls(), c("outputlomod1",
                        "outputmemod1",
                        "outputhimod1",
                        "outputlorescor1",
                        "outputmerescor1",
                        "outputhirescor1",
                        "outputlo1",
                        "outputme1",
                        "outputhi1")))

outputlomod1<-outputlomod
outputlorescor1<-outputlorescor
outputlo1<-outputlo                
outputlomod[[1]]<-outputlomod1[[1]]
outputlomod[[2]]<-outputlomod1[[2]]
outputlomod[[3]]<-outputlomod1[[3]]
outputlorescor[[1]]<-outputlorescor1[[1]]
outputlorescor[[2]]<-outputlorescor1[[2]]
outputlorescor[[3]]<-outputlorescor1[[3]]
outputlo[1:3,]<-outputlo1[1:3,]

outputmemod1<-outputmemod
outputmerescor1<-outputmerescor
outputme1<-outputme
outputmemod[[1]]<-outputmemod1[[1]]
outputmerescor[[1]]<-outputmerescor1[[1]]
outputme[1,]<-outputme1[1,]

outputhimod1<-outputhimod
outputhirescor1<-outputhirescor
outputhi1<-outputhi       
outputhimod[[1]]<-outputhimod1[[1]]
outputhimod[[2]]<-outputhimod1[[2]]
outputhirescor[[1]]<-outputhirescor1[[1]]
outputhirescor[[2]]<-outputhirescor1[[2]]
outputhi[1:2,]<-outputhi1[1:2,]

rm(list=c("outputlomod1",
  "outputmemod1",
  "outputhimod1",
  "outputlorescor1",
  "outputmerescor1",
  "outputhirescor1",
  "outputlo1",
  "outputme1",
  "outputhi1"))








