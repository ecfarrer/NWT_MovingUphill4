
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

outputlos$complexity<-outputlos$tot/outputlos$nodes
mean(outputlos$tot)
std.error(outputlos$tot)
mean(outputlos$complexity)
std.error(outputlos$complexity)


## Med ##

outputmes<-data.frame(rep=rep(NA,10),pos=rep(NA,10),neg=rep(NA,10),tot=rep(NA,10),nodes=rep(NA,10))
outputmemods<-list(NA)
outputmerescors<-list(NA)

for (i in 1:10){
  dim(hmscYme5) #there are three plants, 301 total entering analysis
  set.seed(i+4)
  ind<-c(sort(sample(1:298,270)),299:301)
  hmscYme6<-hmscYme5[,ind]
  dim(hmscYme6)
  
  mod.me1<- boral(y = hmscYme6, X = hmscXme3, lv.control = list(num.lv = 3), family = c(rep("normal",272),rep("negative.binomial",1)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
  rescor.me1 <- get.residual.cor(mod.me1) 
  
  outputmemods[[i]]<-mod.me1 #this saves the data too $X and $y
  outputmerescors[[i]]<-rescor.me1
  
  colMatme<-rescor.me1$sig.correlaton
  colMatme[which(rescor.me1$sig.correlaton>0)]<-1
  colMatme[which(rescor.me1$sig.correlaton<0)]<- -1
  
  graphme1<-graph_from_adjacency_matrix(colMatme, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
  myedgelistme<-data.frame(as_edgelist(graphme1),weight=E(graphme1)$weight) #just the edges
  
  graphme2<-graph.edgelist(as.matrix(myedgelistme[,1:2]),directed=FALSE)
  
  outputmes[i,"rep"]<-i
  outputmes[i,"pos"]<-length(which(myedgelistme$weight==1))
  outputmes[i,"neg"]<-length(which(myedgelistme$weight==-1))
  outputmes[i,"tot"]<-length(which(myedgelistme$weight==1))+length(which(myedgelistme$weight==-1))
  outputmes[i,"nodes"]<-length(V(graphme2))
  print(outputmes)
}

outputmes$complexity<-outputmes$tot/outputmes$nodes
mean(outputmes$tot)
std.error(outputmes$tot)
mean(outputmes$complexity)
std.error(outputmes$complexity)


save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/MovingUphill4_WorkspaceSubsetting2.Rdata")  # 
