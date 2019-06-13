
#Subsetting lo me and hi networks to the same number of species entering the network (boral) analysis

#load MovingUphill4_WorkspaceAnalysis2.Rdata or whatever the final workspace is from the boral analysis


dim(hmscYhi5)

## Lo ##

outputlos<-data.frame(rep=rep(NA,10),pos=rep(NA,10),neg=rep(NA,10),tot=rep(NA,10),nodes=rep(NA,10))
outputlomods<-list(NA)
outputlorescors<-list(NA)

for (i in 3:10){
  dim(hmscYlo5) #there is only one plant
  set.seed(i+4)
  ind<-c(sort(sample(1:305,272)),306)
  hmscYlo6<-hmscYlo5[,ind]
  dim(hmscYlo6)

mod.lo1<- boral(y = hmscYlo6, X = hmscXlo3, lv.control = list(num.lv = 3), family = c(rep("normal",272),rep("negative.binomial",1)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.lo1 <- get.residual.cor(mod.lo1) 

outputlomods[[i]]<-mod.lo1 #this saves the data too $X and $y
outputlorescors[[i]]<-rescor.lo1

colMatlo<-rescor.lo1$sig.correlaton
colMatlo[which(rescor.lo1$sig.correlaton>0)]<-1
colMatlo[which(rescor.lo1$sig.correlaton<0)]<- -1

graphlo1<-graph_from_adjacency_matrix(colMatlo, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelistlo<-data.frame(as_edgelist(graphlo1),weight=E(graphlo1)$weight) #just the edges

graphlo2<-graph.edgelist(as.matrix(myedgelistlo[,1:2]),directed=FALSE)

outputlos[i,"rep"]<-i
outputlos[i,"pos"]<-length(which(myedgelistlo$weight==1))
outputlos[i,"neg"]<-length(which(myedgelistlo$weight==-1))
outputlos[i,"tot"]<-length(which(myedgelistlo$weight==1))+length(which(myedgelistlo$weight==-1))
outputlos[i,"nodes"]<-length(V(graphlo2))
print(outputlos)
}

save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/MovingUphill4_WorkspaceSubsetting.Rdata")  # 
