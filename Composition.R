###### Change in relative abundance ######
datBacS5k2
datITSS5k2
datEukS5k2
datEukN5k2

#change two colnames, they can't have a dash in it
colnames(datBacS5k2)[86]<-"WPS.2"
colnames(datBacS5k2)[37]<-"BHI80.139"
names(which(colSums(datBacS5k2[,34:89])>2))
relBac<-datBacS5k2 %>% 
  dplyr::select(Sample_name,Bacteroidetes,Cyanobacteria,Gemmatimonadetes,Heterotrophic_Acidobacteria,Heterotrophic_Actinobacteria,Heterotrophic_Chloroflexi,Heterotrophic_Planctomycetes,Heterotrophic_Proteobacteria,Heterotrophic_Verrucomicrobia,Photosynthetic_Acidobacteria) %>% #,WPS.2 is unknown function so I'm taking it out
  gather(Taxa,abun,Bacteroidetes:Photosynthetic_Acidobacteria) %>%
  mutate(Taxa = recode_factor(Taxa, Photosynthetic_Acidobacteria="Acidobacteria (P)",Bacteroidetes = "Bacteroidetes (H)",Cyanobacteria="Cyanobacteria (P)",Gemmatimonadetes="Gemmatimonadetes (H)",Heterotrophic_Acidobacteria="Acidobacteria (H)",Heterotrophic_Actinobacteria="Actinobacteria (H)",Heterotrophic_Chloroflexi="Chloroflexi (H)",Heterotrophic_Planctomycetes="Planctomycetes (H)",Heterotrophic_Proteobacteria="Proteobacteria (H)",Heterotrophic_Verrucomicrobia="Verrucomicrobia (H)"))%>%
  mutate(type="A. Bacteria")

names(which(colSums(datITSS5k2[,34:47])>.75))
relITS<-datITSS5k2 %>% 
  dplyr::select(Sample_name,Ascomycota,Basidiomycota,Glomeromycota,Mortierellomycota) %>%
  gather(Taxa,abun,Ascomycota:Mortierellomycota) %>%
  mutate(type="B. Fungi")

#if the label has the word "unknown" in it, then I don't want to plot it. it means that is it unknown at a level higher than phylum (even if I could tell hetero/photo)
sort(colSums(datEukS5k2[,34:62]))
names(which(colSums(datEukS5k2[,34:62])>1))#,Alveolata,Archaeplastida,Photosynthetic_Stramenopiles,Rhizaria
relEukS<-datEukS5k2 %>% 
  dplyr::select(Sample_name,Cercozoa,Charophyta,Chlorophyta,Ciliophora,Heterotrophic_Euglenozoa,Heterotrophic_Stramenopiles,Photosynthetic_Stramenopiles) %>%
  gather(Taxa,abun,Cercozoa,Charophyta,Chlorophyta,Ciliophora,Heterotrophic_Euglenozoa,Heterotrophic_Stramenopiles,Photosynthetic_Stramenopiles) %>%
  mutate(Taxa = recode_factor(Taxa, Cercozoa="Cercozoa (H)",Charophyta = "Charophyta (P)",Chlorophyta="Chlorophyta (P)",Ciliophora="Ciliophora (H)",Heterotrophic_Euglenozoa="Euglenozoa (H)",Heterotrophic_Stramenopiles="Stramenopiles (H)",Photosynthetic_Stramenopiles="Stramenopiles (P)"))%>%
  mutate(type="C. Small Eukaryotes")

names(which(colSums(datEukN5k2[,34:41])>1))
relEukN<-datEukN5k2 %>% 
  dplyr::select(Sample_name,Arthropoda,Nematoda,Rotifera,Tardigrada) %>%
  gather(Taxa,abun,Arthropoda,Nematoda,Rotifera,Tardigrada) %>%
  mutate(type="D. Soil mesofauna")

relALL1<-rbind(relBac,relITS,relEukS,relEukN)#
head(relALL1)

#merge with biogeo6 to get pca1
relALL<-merge(relALL1,biogeo6,"Sample_name")
head(relALL)

#plotdata<-relALL %>%
#  mutate(typeTaxa=paste(type,Taxa)) %>%
#  group_by(Taxa,lomehi,type,typeTaxa) %>%
#  summarise(mean_abun = mean(abun),se_abun=std.error(abun)) 
#  #%>%filter(mean_abun>.04)

#this was weird, maybe something changed in ggplot or dplyr because the colors were messing up and it was listing the legend in alfabetical order by taxa rather than the order in the "plotdata" dataframe. the workaroudn was to set the levels of plotdata$Taxa so they were correct
plotdata<-relALL %>%
  mutate(typeTaxa=paste(type,Taxa)) %>%
  group_by(typeTaxa,Taxa,lomehi,type) %>%
  summarise(mean_abun = mean(abun),se_abun=std.error(abun)) 
plotdata$Taxa<-factor(plotdata$Taxa,levels=unique(plotdata$Taxa))

as.data.frame(plotdata)
plotdata$lomehi<-factor(plotdata$lomehi,levels=c("lo","me","hi"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/relabuntaxavsplantdensitygroupsR2.pdf",width=6.5,height=4.3)#,width=4.3, height=5.3
ggplot(plotdata,aes(x=lomehi,y=mean_abun,group=typeTaxa,color=Taxa))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean_abun+se_abun, ymin=mean_abun-se_abun),width=.15,size=.5)+
  scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
#dev.off()

#10 bacteria, 4 fungi, 7 small euks, 4 large euks
mycols<-c("#D9A125",#yellow
          "#4BC366",#light green
          "#6F94DE",#light blue
          "#B4405E",#red
          "#D185E0",#light purple
          "#659125",#green
          "#ff99a4",#light pink
          "#cf6f23",#orange
          "#5C426C",#dark purple
          "#6768A3",#medium blue last bact

          "#cf6f23",#orange" 
          "#D9A125",#yellow
          "#B4405E",#red
          "#6768A3",#medium blue
          
          "#cf6f23",#orange
          "#659125",#green
          "#4BC366",#light green          
          "#D185E0",#light purple
          "#D9A125",#yellow
          "#5C426C",#dark purple 
          "#008B8B",#greenblue "#ff99a4",#light pink       

          "#cf6f23",#orange
          "#D9A125", #yellow          
          "#5C426C", #dark purple
          "#6F94DE")#light blue


#scatter plots - super messy
head(relALL)

plotdata<-relALL %>%
  mutate(typeTaxa=paste(type,Taxa))
plotdata$Taxa<-factor(plotdata$Taxa,levels=unique(plotdata$Taxa))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/relabuntaxavsplantdensitygroupsBFSLENscatter.pdf",width=6.5,height=6)#,width=4.3, height=5.3
ggplot(plotdata,aes(x=log10(Plant_Dens+1),y=abun,group=typeTaxa,color=Taxa))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  scale_color_manual(values=mycols) +
  geom_point(size=.2)+
  geom_smooth(method=lm,se=F,size=.8) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
dev.off()


#Doing anova on all of the above taxa groups
ind<-length(unique(relALL$Taxa))
anovaoutput<-data.frame(Taxa=rep(NA,ind),F=rep(NA,ind),P=rep(NA,ind))
for(i in 1:ind){
  current.taxa<-levels(relALL$Taxa)[i]
  temp<-relALL %>%
    filter(Taxa==current.taxa)
  mod<-anova(lm(abun~lomehi,data=temp))
  anovaoutput[i,1]<-as.character(current.taxa)
  anovaoutput[i,2]<-mod$`F value`[1]
  anovaoutput[i,3]<-mod$`Pr(>F)`[1]
}
anovaoutput$qval<-p.adjust(anovaoutput$P,method="fdr")
anovaoutput$qval<-format(anovaoutput$qval,scientific=F)

#anovaoutput$Taxa<-factor(anovaoutput$Taxa,levels=unique(plotdata$Taxa))
anovaoutput[order(anovaoutput$Taxa),]








##### Ordination #####
#input file
comm.ord<-comm.bio

#calculate CLR on full datasets for 16S, ITS, EukS, EukN

colnames(comm.ord)
temp<-lapply(colnames(comm.ord),function(x){substr(x,1,1)})
head(which(temp=="I"))
tail(which(temp=="I"))

#separate N, S, 16S, ITS
env<-comm.ord[,1:54]
commN<-comm.ord[,55:504]
commS<-comm.ord[,505:3665]
commB<-comm.ord[,3666:20277]
commI<-comm.ord[,20278:24511]

ind<-which(colSums(commN)>0);length(ind)
commN2<-commN[ind]






mynmds<-metaMDS(comm.dataBac[,32:3430],dist="bray",trymax = 1000)
#old lomehi
col=ifelse(comm.dataBac$lomehi=="lo","lightblue",NA)
col[which(comm.dataBac$lomehi=="me")]<-"dodgerblue"
col[which(comm.dataBac$lomehi=="hi")]<-"darkblue"

plot(scores(mynmds),col=col,pch=21,bg=col)#-scores(mynmds)[,1],scores(mynmds)[,2]
ordiellipse(mynmds,groups=col,col=c("darkblue","dodgerblue","lightblue"),conf=.99999,kind="se",lwd=2)#
legend("bottomright",c("Early","Mid","Late"),col=c("#ab3b57","#5268cb","#639e51"),pch=21,pt.bg=c("#ab3b57","#5268cb","#639e51"),lty=1,bty="n")


#contains all data plus biogeochemistry and plant cover, 75 samples
#biogeo info 53
#N 143
#S 1124
#Bact 3399
#ITS 1122
ordi.bact<-cbind(comm.bio[,1:53],comm.bio[,1321:4719])

mynmds<-metaMDS(ordi.bact[,54:3452],dist="bray",trymax = 1000)
#new lomehi
col=ifelse(ordi.bact$lomehi=="lo","lightblue",NA)
col[which(ordi.bact$lomehi=="me")]<-"dodgerblue"
col[which(ordi.bact$lomehi=="hi")]<-"darkblue"

plot(scores(mynmds),col=col,pch=21,bg=col)#-scores(mynmds)[,1],scores(mynmds)[,2]
ordiellipse(mynmds,groups=col,col=c("darkblue","dodgerblue","lightblue"),conf=.99999,kind="se",lwd=2)#




comm.dataBac[1:5,32:33]
ordi.bact[1:5,54:55]
