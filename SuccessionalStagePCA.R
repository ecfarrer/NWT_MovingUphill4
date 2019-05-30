#Successional stage PCA


##### biogeochemistry #####
biogeo<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Biogeochemistry/CN_enzymes_pH_moisture_whc_Niwot2015t.csv")
head(biogeo)
dim(biogeo)
#note on dataset - I updated this because for the data from microbial biomass (gravimetric moisture, IN, DOC, etc), sample 33 and 34 were mixed up. origianlly sample 33 was missing from the first dataset dorota gave me (CN_enzymes_pH_moisture_whc_Niwot2015.xlsx), and the nubmers fom 34 were really the values for 33. Then I got the All Data file from Dorota and this cleared it up b/c it had both samples 33 adn 34 in it, I also double checked with the Biogeochemistry_Niwot_2015.xlsx file where the calculations were done.

#merge with one of the mapping files
biogeo2<-merge(biogeo,datEukS5otu[,1:31])
#cbind(biogeo2$moisture,biogeo2$WHC)
head(biogeo2)

#I could replace all negative numbers with 0, not sure if I should do this, since they are relative, but it doesn't really make sense to have negative microbial biomass or inorganic N
#biogeo2[biogeo2<0]<-0



##### plant cover #####
#get plant cover data to use as "light" and merge with biogeo
plantcov<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Plants/Niwot_MovingUpHill_plots2015.csv")
names(plantcov)[1]<-"Sample_name"
plantcov$plantcover<-(plantcov$MOSS+plantcov$VEG)/plantcov$TOTAL
plantcov$plantcover[which(plantcov$Sample_name==64)]<-0 #this I'm 100% sure has zero plants, however bare and rock were recorded as 0s. I looked up the sample data sheet and it looks like Sam and I did not sample it because a soil smaple was not taken in 2007, however since it had 0 plants we added that plot and forgot to do cover on it.
head(plantcov)
plantcov$vascularplantcover<-(plantcov$VEG)/plantcov$TOTAL
plantcov$vascularplantcover[which(plantcov$Sample_name==64)]<-0

#take out some columns that I don't want to confuse with mapping file columns
plantcov2<-plantcov%>%select(Sample_name,X,Y,plantcover,vascularplantcover)

head(biogeo2)
biogeo3<-merge(biogeo2,plantcov2,"Sample_name") #with merge, they got put in the right order
rownames(biogeo3)<-biogeo3$X.SampleID
head(biogeo3)
#write.csv(biogeo3,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/biogeo3.csv",row.names=F)


##### PCA #####
#Use these plant variables: they include mosses/liverworts but not lichen
#Plant_Dens
#Plant_Div

#datITSS3otu3$X

biogeo4<-biogeo3%>%
  select(Plant_Dens,Plant_Div,TC,TN,NH4,NO3,MicC,MicN,pH,WHC,moisture,snowdepth,elevation,plantcover)#,IN,DOC,DON
ind<-which(is.na(rowSums(biogeo4))==F)
biogeo5<-biogeo4[ind,]
dim(biogeo5)

mypca<-rda(biogeo5,scale=T,na.action=na.omit)
plot(mypca)
summary(mypca)

scal=2
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/successionpca.pdf")
plot(0,0,type="n",xlab=paste("Axis 1 (",sprintf("%.1f",mypca$CA$eig["PC1"]/mypca$tot.chi*100,3),"%)",sep=""),ylab=paste("Axis 2 (",sprintf("%.1f",mypca$CA$eig["PC2"]/mypca$tot.chi*100,3),"%)",sep=""),cex.lab=1.4,xlim=c(-1,3),ylim=c(-1.5,1.5))
abline(h=0)
abline(v=0)
sorts<-order(scores(mypca,scaling=scal)$sites[,1],decreasing=T)
points(-scores(mypca,scaling=scal)$sites[sorts,1],scores(mypca,scaling=scal)$sites[sorts,2],cex=.9,pch=16,col=rep(c('gray10','gray50','gray80'),each=25))
arrows(0,0,-scores(mypca,scaling=scal)$species[,1],scores(mypca,scaling=scal)$species[,2],length=.05,col=1,lwd=2)
text(-scores(mypca,scaling=scal)$species[,1],scores(mypca,scaling=scal)$species[,2],labels=rownames(scores(mypca,scaling=scal)$species),cex=.9,col=1)
legend(2,1.5,c("Early succession","Mid succession","Late succession"),pch=16,col=c('gray10','gray50','gray80'),bty="n",cex=.9)
#dev.off()


succession<-(-scores(mypca)$sites[,1])

succession2<-data.frame(pca1=succession,X.SampleID=names(succession))
succession3<-succession2[order(succession2$pca1),]
succession3$lomehi<-rep(c('lo','me','hi'),each=25)

plot(biogeo2$snowdepth,biogeo2$Plant_Dens)
plot(biogeo2$snowdepth,biogeo2$pH)
plot(biogeo2$Plant_Dens,biogeo2$pH)

summary(lm(biogeo3$Plant_Dens~biogeo3$WHC))
summary(lm(biogeo3$Plant_Dens~biogeo3$moisture))


###### Merge back with biogeo3 #####
biogeo6<-merge(biogeo3[,-which(names(biogeo3)=="lomehi")],succession3)
head(biogeo6)
cbind(biogeo6$plantcover,biogeo6$Plant_Dens,biogeo6$lomehi)
#it is kind of strange that some of the plots with 0 plants were classified as medium, but there is not much I can do if that's how the pca shakes out. i'll just have to see what the networks look like

#write file for Cliff to make a map of the plots
write.csv(biogeo6,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/networkmsplots.csv")


#### Below analyses were from the previous iteration ####

###### Looking at variation in env characteristics ######

biogeom<-biogeo6noneg%>%
  select(pH:TC,Plant_Dens:Plant_Div,plantcover:pca1,lomehi)%>%
  group_by(lomehi)%>%
  summarise_at(vars(pH:pca1),funs(sd(.,na.rm=T)/mean(.,na.rm=T)))%>% #cv only makes sense for variables with only positive numbers
  #summarise_at(vars(pH:pca1),funs(sd(.,na.rm=T)))%>%
  #summarise_at(vars(pH:pca1),funs(mean(.,na.rm=T)))%>%
  mutate(lomehi=factor(lomehi,levels=c("lo","me","hi")))%>%
  gather(variable,range,pH:pca1)
as.data.frame(biogeom)

range(biogeo2$pH,na.rm=T)

ggplot(biogeom,aes(x=lomehi,y=range,color=variable))+
  labs(x = "",y="Range")+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  facet_wrap(~variable,scales="free")

#using sd
#variables with highest variation in late succession: MicC, MicN, moisture, NH4, pc1,plant density, plant diversity, plantcover,snowdepth,TC,TN,WHC
#variables with highest variation in lo/me or not much different from late: NO3, pH
#there is more variation in th high plots, in regular correlation networks you would think this would create more complex networks b/c there are more niches and you would capture species that cooccur due to shared environment (which is the opposite of what I see anyway). But taking out the effect of env you should be removing that hterogeneity, so I don't think there should be any lasting effect on the networks

#using CV (for positive variables)
range(biogeo2$moisture,na.rm=T)
sort(biogeo6$MicC)
sort(biogeo6noneg$MicC)

biogeo6noneg<-biogeo6%>%
  mutate(MicC=ifelse(MicC>0,MicC,0))%>%
  mutate(MicN=ifelse(MicN>0,MicN,0))%>%
  mutate(NO3=ifelse(NO3>0,NO3,0))%>%
  mutate(NH4=ifelse(NH4>0,NH4,0))

#negatives: MicC, MicN, NO3, NH4, pca1
#putting zero for all negatives from above (not pca1)
#highest variation in late: moisture, snowdepth, TC, TN, WHC, NO3
#highest in lo/me: pH, plant dens, plant div, plant cov,MicC, MicN, NH4 




###### Functional redundancy ######

comm.bio

#hmscY<-comm.bio[,54:7315] #for count data
frsp<-comm.bio[,1287:6139] #for count data, only bacteria

rownames(frsp)<-comm.bio$X.SampleID
frsp[1:10,1:10]

frenzymes<-comm.bio[,c(8:14,53)]
rownames(frenzymes)<-comm.bio$X.SampleID

#select lo/me/hi
ind<-which(frenzymes$lomehi=="hi")
frenzymesb<-frenzymes[ind,]
frspb<-frsp[ind,]

#select species with greater than X (X+1 or more) occurrences and remove lo me hi
ind<-which(colSums(frspb>0)>8)
length(ind)
frspc<-frspb[,ind]
frenzymesc<-frenzymesb[,1:dim(frenzymesb)[2]-1]#
dim(frspc)
dim(frenzymesc)

#remove NA plots (only for ordination)
#ind<-which(is.na(rowSums(frenzymesc))==F)
#frspd<-frspc[ind,]
#frenzymesd<-frenzymesc[ind,]

frenzymese<-scale(frenzymesc)

#myrda<-rda(frspd,frenzymese,scale=T,na.action=na.omit)
#plot(myrda,scaling=1)
#summary(myrda)

#fr<-cbind(frenzymese,frspd)

numspe<-dim(frspc)[2]

frout<-data.frame(enzyme=rep(colnames(frenzymesc),each=numspe),spearmanrho=rep(NA,7*dim(frspc)[2]),spearmanp.value=rep(NA,7*dim(frspc)[2]))

for (e in 1:7){
  for (i in 1:numspe){
    test<-cor.test(frenzymese[,e],frspc[,i],method="spearman",na.action=na.rm)
    frout[i+((e-1)*numspe),2]<-test$estimate
    frout[i+((e-1)*numspe),3]<-test$p.value
  }
}

for (e in 1:7){
  for (i in 1:numspe){
    test<-spearman_test(frenzymese[,e]~frspc[,i],distribution=approximate(B=9999))
    #spearman_test(frspc[,i]~frenzymese[,e],distribution="asymptotic",ties.method="mid-ranks")#this doesn't give the same asymptotic result as cor.test which is odd, but it is clos
    #frout[i+((e-1)*numspe),2]<-test$estimate #it does not seem to calculate rho
    test2<-cor.test(frenzymese[,e],frspc[,i],method="spearman",na.action=na.rm)
    frout[i+((e-1)*numspe),2]<-test2$estimate
    frout[i+((e-1)*numspe),3]<-pvalue(test)
  }
}

lofrout<-frout
lofrout$lomehi<-"lo"
head(lofrout)
lofrout$qval<-p.adjust(lofrout$spearmanp.value,method="fdr")
which(lofrout$spearmanp.value<0.05)
which(lofrout$qval<0.05)
hist(lofrout$qval)
min(lofrout$spearmanp.value)

mefrout<-frout
mefrout$lomehi<-"me"
mefrout$qval<-p.adjust(mefrout$spearmanp.value,method="fdr")

hifrout<-frout
hifrout$lomehi<-"hi"
hifrout$qval<-p.adjust(hifrout$spearmanp.value,method="fdr")
which(hifrout$spearmanp.value<0.05)
which(hifrout$qval<0.05)
min(hifrout$spearmanp.value)
24, 23, 71 = 118

frout2<-rbind(lofrout,mefrout,hifrout)
frout2$qvalall<-p.adjust(frout2$spearmanp.value,method="fdr")
head(frout2)
min(frout2$qvalall)

plot(frenzymese[,1],frspc[,6])

frout3<-frout2%>%
  filter(spearmanrho>0)%>%
  mutate(qvalpos=p.adjust(spearmanp.value,method="fdr"))%>%
  group_by(enzyme, lomehi)%>%
  filter(spearmanp.value<0.05)%>%
  summarise(number=n())

frout3[20,]<-c("LAP","lo",0)
frout3[21,]<-c("PHOS","lo",0)
frout3$number<-as.numeric(frout3$number)

frout4<-frout3%>%
  group_by(lomehi)%>%
  summarise(mean=mean(number),se=std.error(number))

ggplot(frout4,aes(x=lomehi,y=mean,group=lomehi,color=lomehi))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.5)+
  guides(col = guide_legend(ncol = 1))

t.test(frout3$number[frout3$lomehi=="lo"],frout3$number[frout3$lomehi=="hi"],paired=T)
