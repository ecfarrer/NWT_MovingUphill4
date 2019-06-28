
##### Phylogenetic Diversity #####
#Files:
richEukS2
richEukN2
richBac2
richITS2

biogeo6

richBac2$X.SampleID<-rownames(richBac2)
richBac3<-merge(richBac2,biogeo6,"X.SampleID")
richBac3$type<-"1Bacteria"

richITS2$X.SampleID<-rownames(richITS2)
richITS3<-merge(richITS2,biogeo6,"X.SampleID")
richITS3$type<-"2Fungi"

#which chao1
richITS3$PD<-richITS3$Chao1

richEukS2$X.SampleID<-rownames(richEukS2)
richEukS3<-merge(richEukS2,biogeo6,"X.SampleID")
richEukS3$type<-"3Small Eukaryotes"

richEukN2$X.SampleID<-rownames(richEukN2)
richEukN2$X.SampleID<-gsub(pattern = "N", replace = "S", x = richEukN2$X.SampleID)
richEukN3<-merge(richEukN2,biogeo6,"X.SampleID")
richEukN3$type<-"4Soil Mesofauna"


richdata<-rbind(richBac3,richITS3,richEukS3,richEukN3)

richmeans<-richdata%>%
  group_by(type,lomehi)%>%
  summarise(mean=mean(PD),se=std.error(PD))
richmeans$lomehi<-factor(richmeans$lomehi,levels=c("lo","me","hi"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/diversitybysuccessionalstagechao1.pdf",width=3.386,height=3.386) #width=3.386 or 7
ggplot(richmeans,aes(x=lomehi,y=mean,group=type))+
  labs(x = "",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.4)+#.5
  geom_point(size=1.5)+#2
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.4)+#.5
  #scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
#dev.off()

anova(mb<-lm(PD~lomehi,data=richBac3))
anova(mi<-lm(PD~lomehi,data=richITS3))
anova(ms<-lm(PD~lomehi,data=richEukS3))
anova(mn<-lm(PD~lomehi,data=richEukN3))

summary(mb <- aov(PD~lomehi, data = richBac3))
TukeyHSD(mb, "lomehi", ordered = TRUE)
plot(TukeyHSD(mb, "lomehi"))
summary(mb <- aov(PD~lomehi, data = richITS3))
TukeyHSD(mb, "lomehi", ordered = TRUE)
summary(mb <- aov(PD~lomehi, data = richEukS3))
TukeyHSD(mb, "lomehi", ordered = TRUE)
summary(mb <- aov(PD~lomehi, data = richEukN3))
TukeyHSD(mb, "lomehi", ordered = TRUE)

#Succession as continuous
ggplot(richdata,aes(x=log10(Plant_Dens+1),y=SR))+# as.numeric(fert),color=species
  labs(x="Plant density",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  #geom_smooth(method=lm,se=F,size=.8,color="black",formula = y ~ poly(x, 2)) +
  facet_wrap(~type,scales="free")

#could log transform if I wanted to, it looks less skewed
ggplot(richdata,aes(x=log10(pca1+1),y=PD))+# as.numeric(fert),color=species
  labs(x="Succession",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  #geom_smooth(method=lm,se=F,size=.8,color="black",formula = y ~ poly(x, 2)) +
  facet_wrap(~type,scales="free")


#looking into just bacteria pd vs richness
ggplot(pdBac2,aes(x=log10(Plant_Dens+1),y=PD))+# as.numeric(fert),color=species
  labs(x="Plant density",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black")

summary(lm(PD~log10(Plant_Dens+1),data=pdBac2))

ggplot(pdBac2,aes(x=log(pca1+1),y=PD))+# as.numeric(fert),color=species
  labs(x="Succession",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black")

summary(lm(PD~log(pca1+1),data=pdBac2))


richmeans<-pdBac2%>%
  group_by(lomehi)%>%
  summarise(mean=mean(PD),se=std.error(PD))
richmeans$lomehi<-factor(richmeans$lomehi,levels=c("lo","me","hi"))

anova(lm(PD~lomehi,data=richBac2))

ggplot(richmeans,aes(x=lomehi,y=mean,group=1))+
  labs(x = "",y="Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.5)







##### Species richness #####

richBac2$X.SampleID<-rownames(richBac2)
richBac3<-merge(richBac2,biogeo6,"X.SampleID")
richBac3$type<-"1Bacteria"

richITS2$X.SampleID<-rownames(richITS2)
richITS3<-merge(richITS2,biogeo6,"X.SampleID")
richITS3$type<-"2Fungi"

richEukS2$X.SampleID<-rownames(richEukS2)
richEukS3<-merge(richEukS2,biogeo6,"X.SampleID")
richEukS3$type<-"3Small Eukaryotes"

richEukN2$X.SampleID<-rownames(richEukN2)
richEukN2$X.SampleID<-gsub(pattern = "N", replace = "S", x = richEukN2$X.SampleID)
richEukN3<-merge(richEukN2,biogeo6,"X.SampleID")
richEukN3$type<-"4Soil Mesofauna"

richdata<-rbind(richBac3,richITS3,richEukS3,richEukN3)

richmeans<-richdata%>%
  group_by(type,lomehi)%>%
  summarise(mean=mean(Chao1),se=std.error(Chao1))
richmeans$lomehi<-factor(richmeans$lomehi,levels=c("lo","me","hi"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/richnessbysuccessionalstagechao1.pdf",width=3.386,height=3.386) #width=3.386 or 7
ggplot(richmeans,aes(x=lomehi,y=mean,group=type))+
  labs(x = "",y="Taxonomic Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.4)+#.5
  geom_point(size=1.5)+#2
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.4)+#.5
  #scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
dev.off()

anova(lm(Chao1~lomehi,data=richBac3))
anova(lm(Chao1~lomehi,data=richITS3))
anova(lm(Chao1~lomehi,data=richEukS3))
anova(lm(Chao1~lomehi,data=richEukN3))

summary(mb <- aov(Chao1~lomehi, data = richBac3))
TukeyHSD(mb, "lomehi", ordered = TRUE)
plot(TukeyHSD(mb, "lomehi"))
summary(mb <- aov(Chao1~lomehi, data = richITS3))
TukeyHSD(mb, "lomehi", ordered = TRUE)
summary(mb <- aov(Chao1~lomehi, data = richEukS3))
TukeyHSD(mb, "lomehi", ordered = TRUE)
summary(mb <- aov(Chao1~lomehi, data = richEukN3))
TukeyHSD(mb, "lomehi", ordered = TRUE)





##### Evenness #####

richBac2$Evenness<-vegan::diversity(datBacS5cotu[,-c(1:33)])/log(specnumber(datBacS5cotu[,-c(1:33)]))
richBac3<-merge(richBac2,biogeo6,"X.SampleID")
richBac3$type<-"1Bacteria"

richITS2$Evenness<-vegan::diversity(datITSS5cotu[,-c(1:33)])/log(specnumber(datITSS5cotu[,-c(1:33)]))
richITS3<-merge(richITS2,biogeo6,"X.SampleID")
richITS3$type<-"2Fungi"

richEukS2$Evenness<-vegan::diversity(datEukS5cotu[,-c(1:33)])/log(specnumber(datEukS5cotu[,-c(1:33)]))
richEukS3<-merge(richEukS2,biogeo6,"X.SampleID")
richEukS3$type<-"3Small Eukaryotes"

richEukN2$Evenness<-vegan::diversity(datEukN5cotu[,-c(1:33)])/log(specnumber(datEukN5cotu[,-c(1:33)]))
richEukN3<-merge(richEukN2,biogeo6,"X.SampleID")
richEukN3$type<-"4Soil Mesofauna"

richdata<-rbind(richBac3,richITS3,richEukS3,richEukN3)

richmeans<-richdata%>%
  group_by(type,lomehi)%>%
  summarise(mean=mean(Evenness,na.rm=T),se=std.error(Evenness,na.rm=T))
richmeans$lomehi<-factor(richmeans$lomehi,levels=c("lo","me","hi"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/evennessbysuccessionalstage.pdf",width=3.386,height=3.386) #width=3.386 or 7
ggplot(richmeans,aes(x=lomehi,y=mean,group=type))+
  labs(x = "",y="Evenness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.4)+#.5
  geom_point(size=1.5)+#2
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.4)+#.5
  #scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
dev.off()


anova(lm(Evenness~lomehi,data=richBac3))
anova(lm(Evenness~lomehi,data=richITS3))
anova(lm(Evenness~lomehi,data=richEukS3))
anova(lm(Evenness~lomehi,data=richEukN3))

summary(mb <- aov(Evenness~lomehi, data = richBac3))
TukeyHSD(mb, "lomehi", ordered = TRUE)
plot(TukeyHSD(mb, "lomehi"))
summary(mb <- aov(Evenness~lomehi, data = richITS3))
TukeyHSD(mb, "lomehi", ordered = TRUE)
summary(mb <- aov(Evenness~lomehi, data = richEukS3))
TukeyHSD(mb, "lomehi", ordered = TRUE)
summary(mb <- aov(Evenness~lomehi, data = richEukN3))
TukeyHSD(mb, "lomehi", ordered = TRUE)








# Old code when testing DADA2 with different parameters (none of them made a difference)

biogeo8<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/biogeo8.csv")

mapBac<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Bact/515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.txt")
mapITS<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITSsingle/ITS_Niwot_20072015_All_MapFilenewlomehi.txt")

ps <- phyloseq(otu_table(seqtab.nochimb, taxa_are_rows=FALSE), 
               sample_data(mapBac))
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(mapITS))

ps2<-ps%>%
  subset_samples(SampleType=="soil"&year==2015)%>%
  subset_samples(Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)%>%
  filter_taxa(function(x) sum(x) > (1), prune=T)%>%##even if there aren't singletons, it prunes them which is nice 
  rarefy_even_depth(sample.size=1023,rngseed=10,replace=F) #%>% for bac 6919 for its 1056
#transform_sample_counts(function(x) x/sum(x) )
sort(sample_sums(ps2))

plot_richness(ps2, x="VascPlant_Div", measures=c("Chao1", "Observed"))

ps3<-cbind(sample_data(ps2),otu_table(ps2))
ps3$Sample_name<-as.numeric(as.character(ps3$Sample_name))

ps3$richness<-rowSums(ps3[,32:dim(ps3)[2]]>0)
m1<-aggregate.data.frame(ps3$richness,by=list(ps3$lomehi),mean)
se<-aggregate.data.frame(ps3$richness,by=list(ps3$lomehi),std.error)
plot(log(ps3$VascPlant_Dens+1),ps3$richness)
abline(lm(ps3$richness~log(ps3$VascPlant_Dens+1)))
summary(lm(ps3$richness~log(ps3$VascPlant_Dens+1)))
plot(ps3$VascPlant_Div,ps3$richness)
abline(lm(ps3$richness~ps3$VascPlant_Div))

ps3means<-ps3%>%
  group_by(lomehi)%>%
  summarise(mean=mean(richness),se=std.error(richness))
ps3means$lomehi<-factor(ps3means$lomehi,levels=c("lo","me","hi"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/diversitybysuccessionalstage.pdf",width=3.386,height=3.386) #width=3.386 or 7
ggplot(ps3means,aes(x=lomehi,y=mean))+
  labs(x = "",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.4)+#.5
  geom_point(size=1.5)+#2
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.4)






