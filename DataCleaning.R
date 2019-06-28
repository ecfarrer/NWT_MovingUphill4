#Euks: joined paired reads, not truncated (b/c short read and primer should be successfully removed b/c primer comes before poor quality bps), used SILVA all taxa nontruncated database (I used release 111, I don't remember and don't see any notes why I originally used 128 and then switched to 111, an older release, maybe b/c I was already used to that release and the taxonomic assignments and labeling)
#Bact: joined reads, truncated aggressively. used greengenes
#ITS: joined reads using DADA2 within R b/c there were trimming options (not truncating), used UNITE database 


#######Read in OTU data#######

##### Read in euk files, not filtered #####

#first read in otu table and taxonomy file as txt, then merge (because the OTU table does not have taxonomy in it)

otufileEuk<-read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Euks/exported-table/otu_table2.txt",header=T)

head(otufileEuk)

#all ranks, consensus, not truncated, all taxa, 111 release
taxonomyfileEuk<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-taxonomy6/taxonomy2.csv",header=T)

head(taxonomyfileEuk)
colnames(taxonomyfileEuk)[1]<-"OTUID"
colnames(taxonomyfileEuk)[2]<-"taxonomy"

otufileEuk2<-merge(otufileEuk,taxonomyfileEuk,"OTUID")
otufileEuk2<-otufileEuk2[,-dim(otufileEuk2)[2]] #delete the last column which is confidence in taxonomy
head(otufileEuk2)

dim(otufileEuk)
dim(otufileEuk2)
dim(taxonomyfileEuk)

#Write it to a text file so that you can convert it to biom, which can then be import back into R using import_biom()

#all ranks and all taxa, consensus, not truncated, release 111
#write.table(otufileEuk2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table6.txt",sep="\t",row.names = F)

#open otu_table6 in excel and add '#OTU ID' as first cell name
biom convert -i otu_table6.txt -o otu_table6.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

otuEuk <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Euks/exported-table/otu_table6.biom",parseFunction = parse_taxonomy_default)
#I'm not sure what the differences is between parse tax greensgenes and default, the warnings told me to use default

head(otu_table(otuEuk))
head(tax_table(otuEuk))

#Import mapping and tree file
#On work computer I got an error that the hashtag at the beginning of the mapping file like was causing problems, so created a new mapping filel where I deleted it. ok then error stopped, not sure what's going on.
#mapEuk<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Euks/EukBr_Niwot_20072015_All_MapFilenohashtag.txt")

mapEuk<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Euks/EukBr_Niwot_20072015_All_MapFilenewlomehi.txt")

treeEuk<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Euks/exported-rooted-tree/tree.nwk")

datEuk<-merge_phyloseq(otuEuk,mapEuk,treeEuk)



####Filter and remove contaminant sequences with decontam()####

#### Euk soil, 111 release ####
##Filter, filter only soil samples, 2015, take out bacteria, unassigned, archaea, plants, fungi, metazoa, take out singletons
#note!!!! if you filter with subset_taxa and a !=, it will NOT return any rows that are NA, so you always have to do an "or NA" in the statement

datEukS<-datEuk%>%
  subset_samples(SampleType=="soil"|Sample_name%in%c("Control1","Control2","Control3","Control4"))%>%
  subset_samples(year==2015)%>%
  subset_taxa(is.na(Rank1)==T|Rank1!="Bacteria")%>%
  subset_taxa(is.na(Rank1)==T|Rank1!="Unassigned")%>%
  subset_taxa(is.na(Rank1)==T|Rank1!="Archaea")%>%
  subset_taxa(is.na(Rank3)==T|Rank3!="__Fungi")%>%
  subset_taxa(is.na(Rank3)==T|Rank3!="__Metazoa")%>%
  subset_taxa(is.na(Rank7)==T|Rank7!="__Embryophyta")%>%
  #subset_samples(Sample_name!=126&Sample_name!=5&Sample_name!=34)%>% none are super low
  filter_taxa(function(x) sum(x) > (1), prune=T) #there are no singletons

unique(tax_table(datEukS)[,"Rank3"])

min(sample_sums(datEukS))
sort(sample_sums(datEukS))

sample_data(datEukS)$Sample_or_Control<-c(rep("True Sample",97),rep("Control Sample",4))
tail(sample_data(datEukS))

dfEuk <- as.data.frame(sample_data(datEukS)) # Put sample_data into a ggplot-friendly data.frame
dfEuk$LibrarySize <- sample_sums(datEukS)
dfEuk <- dfEuk[order(dfEuk$LibrarySize),]
dfEuk$Index <- seq(nrow(dfEuk))
ggplot(data=dfEuk, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(datEukS)$is.neg <- sample_data(datEukS)$Sample_or_Control == "Control Sample"
contamdf.prevEuk <- isContaminant(datEukS, method="prevalence", neg="is.neg")
table(contamdf.prevEuk$contaminant)

datEukS.pa <- transform_sample_counts(datEukS, function(abund) 1*(abund>0))
datEukS.pa.neg <- prune_samples(sample_data(datEukS.pa)$Sample_or_Control == "Control Sample", datEukS.pa)
datEukS.pa.pos <- prune_samples(sample_data(datEukS.pa)$Sample_or_Control == "True Sample", datEukS.pa)
# Make data.frame of prevalence in positive and negative samples
datEukSdf.pa <- data.frame(pa.pos=taxa_sums(datEukS.pa.pos), pa.neg=taxa_sums(datEukS.pa.neg),
                           contaminant=contamdf.prevEuk$contaminant)
ggplot(data=datEukSdf.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#no contaminants! there are three otus found in control samples, but they are not found in any real samples

datEukS2 <- prune_taxa(!contamdf.prevEuk$contaminant, datEukS)
datEukS3 <-datEukS2 %>%
  subset_samples(SampleType=="soil")%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

#I checked that doing all this did not change the data for datEuk from before when I didn't do contaminant filtering, the only thing that changed was that 2 columns were added to the sample data

sort(unique(tax_table(datEukS3)[,"Rank7"]))

length(which(is.na(tax_table(datEukS3)[,"Rank2"])))
#1680 of the 3932 taxa are just assigned euk;__

#sample 61 has 680 reads, next sample76 has 871 reads, so could delete 61
min(sample_sums(datEukS3))
sort(sample_sums(datEukS3))
which(sample_data(datEukS3)$X.SampleID=="S.81.2015")
#sample 81 is not included



###### Euk Nematode Samples 111 release #####
##Filter, filter only nematode samples, filter sample 2A (duplicate), only metazoa, take out craniata (chordata here, chordata is a higher level than craniata but the only chordata here are vertebrates), take out singletons (there are none)

datEukN<-datEuk%>%
  subset_samples(SampleType=="nematode"|Sample_name%in%c("Control5","Control6"))%>%
  subset_samples(X.SampleID!="N.2A.2015")%>%
  subset_taxa(Rank3=="__Metazoa")%>%
  subset_taxa(is.na(Rank4)==T|Rank4!="__Craniata")%>%
  filter_taxa(function(x) sum(x) > 1, prune=T) #there are no singletons
datEukN

#subset_samples(Sample_name!=126&Sample_name!=5&Sample_name!=34)%>% none are super low
min(sample_sums(datEukN))
sort(sample_sums(datEukN))

sample_data(datEukN)$Sample_or_Control<-c(rep("True Sample",96),rep("Control Sample",2))
tail(sample_data(datEukN))

dfEukN <- as.data.frame(sample_data(datEukN)) # Put sample_data into a ggplot-friendly data.frame
dfEukN$LibrarySize <- sample_sums(datEukN)
dfEukN <- dfEukN[order(dfEukN$LibrarySize),]
dfEukN$Index <- seq(nrow(dfEukN))
ggplot(data=dfEukN, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(datEukN)$is.neg <- sample_data(datEukN)$Sample_or_Control == "Control Sample"
contamdf.prevEukN <- isContaminant(datEukN, method="prevalence", neg="is.neg")
table(contamdf.prevEukN$contaminant)

datEukN.pa <- transform_sample_counts(datEukN, function(abund) 1*(abund>0))
datEukN.pa.neg <- prune_samples(sample_data(datEukN.pa)$Sample_or_Control == "Control Sample", datEukN.pa)
datEukN.pa.pos <- prune_samples(sample_data(datEukN.pa)$Sample_or_Control == "True Sample", datEukN.pa)
# Make data.frame of prevalence in positive and negative samples
datEukNdf.pa <- data.frame(pa.pos=taxa_sums(datEukN.pa.pos), pa.neg=taxa_sums(datEukN.pa.neg),
                           contaminant=contamdf.prevEukN$contaminant)
ggplot(data=datEukNdf.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#one contaminant! 

datEukN2 <- prune_taxa(!contamdf.prevEukN$contaminant, datEukN)
datEukN3 <-datEukN2 %>%
  subset_samples(SampleType=="nematode")%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

head(tax_table(datEukN3))
sort(unique(tax_table(datEukN3)[,"Rank4"]))

#sample 78 has only 70 reads so delete; however this could be true b/c 78 is a zero plant plot so there are probably not a lot of nematodes, next sample 3 has 700 (3 is also a zero plant plot), next 83 has 2038 (also zero plants)
min(sample_sums(datEukN3))
sort(sample_sums(datEukN3))
which(sample_data(datEukN3)$X.SampleID=="S.33.2015")
which(sample_data(datEukN3)$X.SampleID=="S.56.2015")
#33 and 56 are missing







###### Bacteria #######
#paired ends
otufileBac<-read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Bact/exported-tableN/otu_table2.txt",header=T)

head(otufileBac)

#taxonomy from greengeenes
#paired ends
taxonomyfileBac<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Bact/exported-taxonomy_ggN/taxonomy2.csv",header=T)

head(taxonomyfileBac)
colnames(taxonomyfileBac)[1]<-"OTUID"
colnames(taxonomyfileBac)[2]<-"taxonomy"

otufileBac2<-merge(otufileBac,taxonomyfileBac,"OTUID")
otufileBac2<-otufileBac2[,-dim(otufileBac2)[2]] #delete the last column which is confidence in taxonomy
head(otufileBac2)

dim(otufileBac)
dim(otufileBac2)
dim(taxonomyfileBac)

#min(otufileEuk2$Confidence)

#Write it to a text file so that you can convert it to biom, which can then be import back into R using import_biom()
#paired end 
write.table(otufileBac2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Bact/exported-tableN/otu_table3.txt",sep="\t",row.names = F)

#open otu_table3 in excel and add '#OTU ID' as first cell name
biom convert -i otu_table3.txt -o otu_table3.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

#paired
otuBac <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Bact/exported-tableN/otu_table3.biom",parseFunction = parse_taxonomy_default)

#I'm not sure what the differences is between parse tax greensgenes and default, the warnings told me to use default

head(otu_table(otuBac))
head(tax_table(otuBac))


#Import mapping and tree file
mapBac<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Bact/515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.txt")

#paired
treeBac<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Bact/exported-rooted-treeN/tree.nwk")

datBac<-merge_phyloseq(otuBac,mapBac,treeBac)


####Filter and remove contaminant sequences with decontam()####
##Filter only soil samples, 2015, keep only bacteria and archaea (this filters unassigned), filter mitochondira and chloroplasts, take out singletons
#note!!!! if you filter with subset_taxa and a !=, it will NOT return any rows that are NA, so you always have to do an "or NA" in the statement

datBacS<-datBac%>%
  subset_samples(SampleType=="soil"|Sample_name%in%c("Control1","Control2","Control3","Control4"))%>%
  subset_samples(year==2015)%>%
  subset_samples(Sample_name!=126&Sample_name!=5&Sample_name!=34)%>%
  subset_taxa(Rank1=="k__Archaea"|Rank1=="k__Bacteria")%>%
  subset_taxa(is.na(Rank3)==T|Rank3!="c__Chloroplast")%>%
  subset_taxa(is.na(Rank5)==T|Rank5!="f__mitochondria")%>%
  filter_taxa(function(x) sum(x) > (1), prune=T) #there are no singletons, but this prunes off taxa that are 0

unique(tax_table(datBac)[,"Rank3"])

#sample 5, 34, 126 have low # reads, 2, 4, 17 respectively
min(sample_sums(datBacS))
sort(sample_sums(datBacS))

sample_data(datBacS)$Sample_or_Control<-c(rep("True Sample",95),rep("Control Sample",4))
tail(sample_data(datBacS))

dfBac <- as.data.frame(sample_data(datBacS)) # Put sample_data into a ggplot-friendly data.frame
dfBac$LibrarySize <- sample_sums(datBacS)
dfBac <- dfBac[order(dfBac$LibrarySize),]
dfBac$Index <- seq(nrow(dfBac))
ggplot(data=dfBac, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(datBacS)$is.neg <- sample_data(datBacS)$Sample_or_Control == "Control Sample"
contamdf.prevBac <- isContaminant(datBacS, method="prevalence", neg="is.neg")
table(contamdf.prevBac$contaminant)

datBacS.pa <- transform_sample_counts(datBacS, function(abund) 1*(abund>0))
datBacS.pa.neg <- prune_samples(sample_data(datBacS.pa)$Sample_or_Control == "Control Sample", datBacS.pa)
datBacS.pa.pos <- prune_samples(sample_data(datBacS.pa)$Sample_or_Control == "True Sample", datBacS.pa)
# Make data.frame of prevalence in positive and negative samples
datBacSdf.pa <- data.frame(pa.pos=taxa_sums(datBacS.pa.pos), pa.neg=taxa_sums(datBacS.pa.neg),
                    contaminant=contamdf.prevBac$contaminant)
ggplot(data=datBacSdf.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

datBacS2 <- prune_taxa(!contamdf.prevBac$contaminant, datBacS)
datBacS3 <-datBacS2 %>%
  subset_samples(SampleType=="soil")%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

tax_table(datBacS3)




###### ITS #######
#For some annoying reason, this takes about 15 minutes to read in (maybe the super long column names?)
otufileITS<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/exported-dada2/seqtab.nochim.csv",header=T, row.names=1)

head(otufileITS)

taxonomyfileITS<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/exported-dada2/taxa3.csv",header=T,row.names = 1)
taxonomyfileITS2<-as.matrix(taxonomyfileITS)
colnames(taxonomyfileITS2)<-c("Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7")

mapITS<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITSsingle/ITS_Niwot_20072015_All_MapFilenewlomehi.txt")

datITS <- phyloseq(otu_table(otufileITS, taxa_are_rows=FALSE), 
               sample_data(mapITS),tax_table(taxonomyfileITS2))

head(otu_table(datITS))
head(tax_table(datITS))
head(tax_table(otuBac))



####Filter and remove contaminant sequences with decontam()####
##Filter, filter only soil samples, 2015, there is only fungi here
#note!!!! if you filter with subset_taxa and a !=, it will NOT return any rows that are NA, so you always have to do an "or NA" in the statement

datITSS<-datITS%>%
  subset_samples(SampleType=="soil"|Sample_name%in%c("Control1","Control2","Control3","Control4"))%>%
  subset_samples(year==2015)%>%
  #subset_taxa(is.na(Rank2)==F)%>% #this takes out the Fungi;__ taxa, I could or not take them out. they could be real or they could be artifacts (i.e. primer errors). I'll leave them in for now. there are 3956 taxa like this - that's a lot!
  filter_taxa(function(x) sum(x) > (0), prune=T)  #there were 2 singletons, I'll keep them in
datITSS

print(as.data.frame(unique(tax_table(datITS)[,"Rank1"])),row.names=F)

#all have modest number of reads
min(sample_sums(datITSS))
sort(sample_sums(datITSS))

sample_data(datITSS)$Sample_or_Control<-c(rep("True Sample",97),rep("Control Sample",4))
tail(sample_data(datITSS))

dfITS <- as.data.frame(sample_data(datITSS)) # Put sample_data into a ggplot-friendly data.frame
dfITS$LibrarySize <- sample_sums(datITSS)
dfITS <- dfITS[order(dfITS$LibrarySize),]
dfITS$Index <- seq(nrow(dfITS))
ggplot(data=dfITS, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(datITSS)$is.neg <- sample_data(datITSS)$Sample_or_Control == "Control Sample"
contamdf.prevITS <- isContaminant(datITSS, method="prevalence", neg="is.neg")
table(contamdf.prevITS$contaminant)

#No contaminants!

datITSS.pa <- transform_sample_counts(datITSS, function(abund) 1*(abund>0))
datITSS.pa.neg <- prune_samples(sample_data(datITSS.pa)$Sample_or_Control == "Control Sample", datITSS.pa)
datITSS.pa.pos <- prune_samples(sample_data(datITSS.pa)$Sample_or_Control == "True Sample", datITSS.pa)
# Make data.frame of prevalence in positive and negative samples
datITSSdf.pa <- data.frame(pa.pos=taxa_sums(datITSS.pa.pos), pa.neg=taxa_sums(datITSS.pa.neg),
                           contaminant=contamdf.prevITS$contaminant)
ggplot(data=datITSSdf.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

datITSS2 <- prune_taxa(!contamdf.prevITS$contaminant, datITSS)
datITSS3 <-datITSS2 %>%
  subset_samples(SampleType=="soil")%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

#sample 112 has kind of low # reads, 13102301, the next lowest sample is 56 at 4145 reads but I think I'll keep it since I'm deleting too many samples from the other runs
min(sample_sums(datITSS3))
sort(sample_sums(datITSS3))
setdiff(sample_data(datBacS3)$Sample_name,sample_data(datITSS3)$Sample_name)
#ITS is missing 126, but we are already deleting that b/c low reads for bacteria








sort(sample_sums(datEukS3))
sort(sample_sums(datEukN3))
sort(sample_sums(datBacS3))
sort(sample_sums(datITSS3))
datEukS3
datEukN3
datBacS3
datITSS3




##### Remove samples across all datasets #####
#(from before) For euksS sample 81 did not amplify and 61 had low # reads. for euksN sample 33 and 56 did not have enough soil so were not done, and 78 had low # reads. for ITS 126 didnt amplify. for bacteria, samples 5,34,126 did not amplify. should have 90 samples left

#From EukS 61 (moderately low reads, could keep), 81 (no reads)
#From EukN 78 (very low reads), 33 (no reads), 56 (no reads)
#From Bac 5 (very low reads), 34 (very low reads), 126 (very low reads)
#From ITS 126 (no reads)

datEukS4<-datEukS3%>%
  subset_samples(Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)
sum(otu_table(datEukS2))

datEukN4<-datEukN3%>%
  subset_samples(Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

datBacS4<-datBacS3%>%
  subset_samples(Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

datITSS4<-datITSS3%>%
  subset_samples(Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)



##### Rarefy and transform to relative abundance #####
min(sample_sums(datEukS4))#rarefy to 871 
min(sample_sums(datEukN4))#rarefy to 700
min(sample_sums(datBacS4))#rarefy to 7987
min(sample_sums(datITSS4))#rarefy to 1023

datEukS5<-datEukS4%>%
  rarefy_even_depth(sample.size=min(sample_sums(datEukS4)),rngseed=10,replace=F) %>%
  transform_sample_counts(function(x) x/sum(x) )
#609 OTUs were removed because they are no longer present in any sample after random subsampling, 

datEukN5<-datEukN4%>%
  rarefy_even_depth(sample.size=min(sample_sums(datEukN4)),rngseed=10,replace=F)%>%
  transform_sample_counts(function(x) x/sum(x) )
#80 OTUs were removed because they are no longer present in any sample after random subsampling, 

datBacS5<-datBacS4%>%
  rarefy_even_depth(sample.size=min(sample_sums(datBacS4)),rngseed=10,replace=F)%>%
  transform_sample_counts(function(x) x/sum(x) )
#1524 OTUs were removed because they are no longer present in any sample after random subsampling, 

datITSS5<-datITSS4%>%
  rarefy_even_depth(sample.size=min(sample_sums(datITSS4)),rngseed=10,replace=F)%>%
  transform_sample_counts(function(x) x/sum(x) )
#2475 OTUs were removed because they are no longer present in any sample after random subsampling, 


#rarefying but not calculating relative abundance
datEukS5c<-datEukS4%>%
  rarefy_even_depth(sample.size=min(sample_sums(datEukS4)),rngseed=10,replace=F)
#609 OTUs were removed because they are no longer present in any sample after random subsampling, 

datEukN5c<-datEukN4%>%
  rarefy_even_depth(sample.size=min(sample_sums(datEukN4)),rngseed=10,replace=F)
#80 OTUs were removed because they are no longer present in any sample after random subsampling, 

datBacS5c<-datBacS4%>%
  rarefy_even_depth(sample.size=min(sample_sums(datBacS4)),rngseed=10,replace=F)
#1524 OTUs were removed because they are no longer present in any sample after random subsampling, 

datITSS5c<-datITSS4%>%
  rarefy_even_depth(sample.size=min(sample_sums(datITSS4)),rngseed=10,replace=F)
#2475 OTUs were removed because they are no longer present in any sample after random subsampling, 






##### Make a label column with kingdom/phylum level labels #####

##### Euks Soil #####

#labelsEukS<-data.frame(rep("Eukaryota",dim(tax_table(datEukS3))[1]))

#note - indented notes are from NWT_MovingUphill2
#Amoebozoa (phylum)
#Alveolata (Subkingdom, unranked clade, superphylum) - splitting into phyla, as I understand it the only possibly photosynthetic phylum is Dinoflagellata but they include both photo and heteros, and the 4 organisms I have could be either. So I'm calling Dinoflagellata unknown, all other phyla heterotrophs, and those things classified just to Alveolata unknown.
  #photosynthetic Alveolata (kingdom) (phylum Dinoflagellata: mostly photosynthetic; nonphotosynthetic Protoperidinium, both SL163A10 and Pfiesteria (can if it eats an alga), unknown D244),) only 4 taxa, not sure if they are Ps but I'll leave them here. 
  #nonphotosynthetic Alveolata (phyla Ciliophora(predators), protalveolata (predatory flagellates), apicomplexa (parasites)), BOLA553 is a near apicomplexan, not sure what SCM37C52 is so I put it here
#Archaeplastida (major_clade), phyla chlorophyta, charophyta, an unknowns/unclassified, technically a few chlorophyta are non-Ps, but it looks like all charophyta are Ps, they are such a large group and photosynthesis tends to be dominant so I will classify them all as Ps
#Rhizaria (subkingdom, unrankd clade, superphylum): unicellular amoeboid euks
#Holozoa(unranked: all animals, not fungi), includes metazoa (animals)
#Excavata: two phyla, Euglenozoa - not all Euglenida are Ps, but  one I have Euglena viridis is. the others Petalomonas are heterotrophic; all other Euglenozoa are Kinetoplastida which are heterotrophic. and phylum Heterolobosea which are heterotrophs
  #photosynthetic Excavata major clade (was Discoba kingdom) (class Euglenida: mostly but not all photosynthetic), 
  #nonphotosnthetic Excavata (was Discoba) (in discoba, phylum Heterolobosea: parasites, free living, symbiotic, amoeba-like, and phylum Jakoba heterotrophic), superphylum metamonada, doesn't have kingdom, parasites not Ps. Kinetoplastea: not Ps, Tetramitia: most likely not Ps
#we dont have any cryptophyceae
  #photosynthetic Cryptophyceae is a kingdom, they have one or two chloroplasts, except for Chilomonas, which has leucoplasts (not pigmented not for Ps) and Goniomonas (formerly Cyathomonas) which lacks plastids entirely.P1-31 (i can't find much info on this, some sites called them picoplankton and most of the cryptophyceae are photosynthetic). Cryptomonadales can't find much about it but assume Ps. __SA1-3C06 is at rank2 but notes in silva say "cryptophyte"
  #nonphotosynthetic cryptophyceae - Chilomonas and Goniomonas (formerly Cyathomonas)
#Haptophyta phylum, I think all are Ps
#Centrohelida order - I can't find it anymore but supposedly encyclopedia of life said some species have photosynthetic symbionts but I can't find any info on Ps taxa. I don't have many only like 8 taxa so I will leave them as centrohelida and assume they are nonPs. "free-living predatory protists" from Burke et al 2016
#Kathablepharidae class (but doesnt have a phylum) - heterotrophs
#NonPs Stramenopiles (no rank, phylum): MAST-12 (Karolina Kolodziej, Thorsten Stoeck 2007), Labyrinthulomycetes,Bicosoecida,Peronosporomycetes,Hyphochytriales
#Ps Stramenopiles: Diatomea, Eustigmatales, Xanthophyceae, Chrysophyceae, Raphidophyceae

#I separated out the telonema, breviatea and Apusomonadidae, even though these are very rare. the main contributer to "heterotrophic eukaryota" was unknown opithokonts. The Apusomonadidae is kind of abundant
#unknown opisthokonta. this is the huge group that includes (heterotrophic) holozoa, animals, fungi, and single celled euks like choanoflagellates. thus I shouldn't put this on a bargraph since it is too broad and unkown but I can still use it in networks
#Telonema - a genus sometimes assigned to a unique phylum or class, nonps protist phylum
#Breviatea class (amoeba like)
#Apusomonadidae (sister to breviatea, nonPs, in a phylum Apusozoa)
#heterotrophic eukaryota (things without a kingdom):Telonema (nonPs protist phylum);Breviatea ;Apusomonadidae sister to the breviatea nonPs

unique(tax_table(datEukS5)[,"Rank2"])
unique(tax_table(datEukS5)[,"Rank3"])

labelsEukS<-tax_table(datEukS5)[,"Rank3"]#

ind<-which(tax_table(datEukS5)[,"Rank2"]=="__Amoebozoa")
labelsEukS[ind]<-"Amoebozoa"

#splitting the Alveolata into phyla
#making dinoflagellates Ps Alveolata (from NWT_MovingUphill2)
#ind<-which(tax_table(datEukS3)[,"Rank4"]=="__Dinoflagellata")#
#labelsEukS[ind]<-"Photosynthetic_Alveolata"
#ind<-which(tax_table(datEukS3)[,"Rank4"]!="__Dinoflagellata"&tax_table(datEukS3)[,"Rank3"]=="__Alveolata")
#labelsEukS[ind]<-"Nonphotosynthetic_Alveolata"

tax_table(datEukS5)[ind,]

ind<-which(tax_table(datEukS5)[,"Rank4"]=="__Dinoflagellata");ind
labelsEukS[ind]<-"Dinoflagellata"
ind<-which(tax_table(datEukS5)[,"Rank4"]=="__Apicomplexa");ind
labelsEukS[ind]<-"Apicomplexa"
ind<-which(tax_table(datEukS5)[,"Rank4"]=="__Protalveolata");ind
labelsEukS[ind]<-"Protalveolata"
ind<-which(tax_table(datEukS5)[,"Rank4"]=="__Ciliophora");ind
labelsEukS[ind]<-"Ciliophora"
ind<-which(tax_table(datEukS5)[,"Rank3"]=="__Alveolata"&is.na(tax_table(datEukS5)[,"Rank4"])==T);ind#I'm calling these unknown
labelsEukS[ind]<-"Unknown_Alveolata"

#splitting the Archaeplastida into phyla
ind<-which(tax_table(datEukS5)[,"Rank4"]=="__Chlorophyta");ind
labelsEukS[ind]<-"Chlorophyta"
ind<-which(tax_table(datEukS5)[,"Rank4"]=="__Charophyta");ind
labelsEukS[ind]<-"Charophyta"
ind<-which(tax_table(datEukS5)[,"Rank2"]=="__Archaeplastida"&is.na(tax_table(datEukS5)[,"Rank3"])==T);ind
labelsEukS[ind]<-"Photosynthetic_Unknown_Archaeplastida"
ind<-which(tax_table(datEukS5)[,"Rank2"]=="__Archaeplastida"&is.na(tax_table(datEukS5)[,"Rank4"])==T);ind #this includes things only identified to Archaeplastida and the Chloroplastida which is not a phylum but unranked and includes land plants (streptophyta and Chlorophyta green algae)
labelsEukS[ind]<-"Photosynthetic_Unknown_Archaeplastida"

#splitting the Rhizaria (an infrakingdom) into phyla
#ind<-which(tax_table(datEukS3)[,"Rank3"]=="__Rhizaria");ind
#labelsEukS[ind]<-"Rhizaria"
ind<-which(tax_table(datEukS5)[,"Rank4"]=="__Cercozoa");ind
labelsEukS[ind]<-"Cercozoa"
ind<-which(tax_table(datEukS5)[,"Rank4"]=="__Foraminifera");ind
labelsEukS[ind]<-"Foraminifera"

#splitting the holozoa 
#ind<-which(tax_table(datEukS3)[,"Rank3"]=="__Holozoa");ind
#labelsEukS[ind]<-"Holozoa"

ind<-which(tax_table(datEukS5)[,"Rank4"]=="__Choanomonada");ind
labelsEukS[ind]<-"Choanomonada"
ind<-which(tax_table(datEukS5)[,"Rank4"]=="__Ichthyosporea");ind
labelsEukS[ind]<-"Ichthyosporea"
ind<-which(tax_table(datEukS5)[,"Rank4"]=="__Filasterea");ind
labelsEukS[ind]<-"Filasterea"
ind<-which(tax_table(datEukS5)[,"Rank3"]=="__Holozoa"&is.na(tax_table(datEukS5)[,"Rank4"])==T);ind
labelsEukS[ind]<-"Heterotrophic_Unknown_Holozoa"

#splitting the Excavata
#ind<-which(tax_table(datEukS3)[,"Rank6"]=="__Euglenida")
#labelsEukS[ind]<-"Photosynthetic_Excavata"
#ind<-which(tax_table(datEukS3)[,"Rank6"]!="__Euglenida"&tax_table(datEukS3)[,"Rank2"]=="__Excavata")
#labelsEukS[ind]<-"Nonphotosynthetic_Excavata"

ind<-which(tax_table(datEukS5)[,"Rank11"]=="__Euglena_viridis");ind
labelsEukS[ind]<-"Photosynthetic_Euglenozoa"
ind<-which(tax_table(datEukS5)[,"Rank8"]=="__Petalomonas");ind
labelsEukS[ind]<-"Heterotrophic_Euglenozoa"
ind<-which(tax_table(datEukS5)[,"Rank6"]!="__Euglenida"&tax_table(datEukS5)[,"Rank5"]=="__Euglenozoa");ind
labelsEukS[ind]<-"Heterotrophic_Euglenozoa"
ind<-which(tax_table(datEukS5)[,"Rank5"]=="__Heterolobosea");ind
labelsEukS[ind]<-"Heterolobosea"

ind<-which(tax_table(datEukS5)[,"Rank2"]=="__Haptophyta");ind
labelsEukS[ind]<-"Haptophyta"
ind<-which(tax_table(datEukS5)[,"Rank2"]=="__Centrohelida");ind
labelsEukS[ind]<-"Centrohelida"
ind<-which(tax_table(datEukS5)[,"Rank2"]=="__Kathablepharidae")
labelsEukS[ind]<-"Kathablepharidae"

ind<-which(tax_table(datEukS5)[,"Rank3"]=="__Stramenopiles"&tax_table(datEukS5)[,"Rank4"]%in%c("__MAST-12","__Labyrinthulomycetes","__Bicosoecida","__Peronosporomycetes","__Hyphochytriales"))
labelsEukS[ind]<-"Heterotrophic_Stramenopiles"
ind<-which(tax_table(datEukS5)[,"Rank3"]=="__Stramenopiles"&tax_table(datEukS5)[,"Rank4"]%in%c("__Diatomea","__Eustigmatales","__Xanthophyceae","__Chrysophyceae","__Raphidophyceae"))
labelsEukS[ind]<-"Photosynthetic_Stramenopiles"
ind<-which(tax_table(datEukS5)[,"Rank3"]=="__Stramenopiles"&is.na(tax_table(datEukS5)[,"Rank4"])==T);ind
labelsEukS[ind]<-"Unknown_Stramenopiles"

ind<-which(tax_table(datEukS5)[,"Rank2"]=="__SAR"&is.na(tax_table(datEukS5)[,"Rank3"])==T);ind
labelsEukS[ind]<-"Unknown_Eukaryota"
ind<-which(is.na(tax_table(datEukS5)[,"Rank2"])==T)
labelsEukS[ind]<-"Unknown_Eukaryota"

#Breaking up the Breviatea, Telomena, Apusomonadidae
#ind<-which(tax_table(datEukS3)[,"Rank3"]=="__Breviatea"|tax_table(datEukS3)[,"Rank3"]=="__Telonema"|tax_table(datEukS3)[,"Rank3"]=="__Apusomonadidae")
#labelsEukS[ind]<-"Heterotrophic_Eukarya"
ind<-which(tax_table(datEukS5)[,"Rank3"]=="__Breviatea");ind
labelsEukS[ind]<-"Breviatea"
ind<-which(tax_table(datEukS5)[,"Rank3"]=="__Telonema");ind
labelsEukS[ind]<-"Telonema"
ind<-which(tax_table(datEukS5)[,"Rank3"]=="__Apusomonadidae");ind
labelsEukS[ind]<-"Apusomonadidae"
ind<-which(tax_table(datEukS5)[,"Rank2"]=="__Opisthokonta"&is.na(tax_table(datEukS5)[,"Rank3"])==T);ind
labelsEukS[ind]<-"Heterotrophic_Unknown_Opisthokonta"

unique(tax_table(datEukS5)[which(tax_table(datEukS5)[,"Rank2"]=="__Archaeplastida"),])
#tax_table(datEukS3)[which(tax_table(datEukS3)[,"Rank3"]=="__Apusomonadidae"),]
#ind<-which(is.na(labelsEukS)==T)
#tax_table(datEukS3)[ind,]


#I could separate excavata into kingdoms/super phyla and also separate archaeplastida into its 3 kingdoms
unique(labelsEukS)
colnames(labelsEukS)<-"labels"
labelsEukS2<-as.data.frame(labelsEukS)
unique(labelsEukS2$labels)
labelsEukS2$group<-"Eukaryota"

labelsEukS2$group2<-NA
ind<-which(labelsEukS2$labels=="Photosynthetic_Unknown_Archaeplastida"|labelsEukS2$labels=="Photosynthetic_Stramenopiles"|labelsEukS2$labels=="Haptophyta"|labelsEukS2$labels=="Photosynthetic_Euglenozoa"|labelsEukS2$labels=="Chlorophyta"|labelsEukS2$labels=="Charophyta");ind
labelsEukS2$group2[ind]<-"PhotosyntheticEukaryota"                    
ind<-which(labelsEukS2$labels=="Unknown_Eukaryota"|labelsEukS2$labels=="Unknown_Stramenopiles"|labelsEukS2$labels=="Unknown_Alveolata"|labelsEukS2$labels=="Dinoflagellata");ind
labelsEukS2$group2[ind]<-"UnknownEukaryota"
labelsEukS2$group2[which(is.na(labelsEukS2$group2)==T)]<-"HeterotrophicEukaryota"

head(labelsEukS2)

dim(tax_table(datEukS5))
head(tax_table(datEukS5))
#tax_table(datEukS3)<-tax_table(datEukS3)[,-15]

#columns 12 and 13 are all NAs
unique(tax_table(datEukS5)[,11])
labelsEukS2$taxstring<-paste(tax_table(datEukS5)[,1],tax_table(datEukS5)[,2],tax_table(datEukS5)[,3],tax_table(datEukS5)[,4],tax_table(datEukS5)[,5],tax_table(datEukS5)[,6],tax_table(datEukS5)[,7],tax_table(datEukS5)[,8],tax_table(datEukS5)[,9],tax_table(datEukS5)[,10],tax_table(datEukS5)[,11],sep=";")

#replace tax table, I'm waiting on this and will only replace it if I find that I need to later on
#dim(tax_table(datEukS5))
#tax_table(datEukS5)<-cbind(tax_table(datEukS5),labelsEukS)

#look at tree of abundant taxa
#myTaxa<-c(names(sort(taxa_sums(datEukr2),decreasing=T)))[1:100]
#ex2 = prune_taxa(myTaxa, datEukr2)
#plot_tree(ex2, label.tips = "taxa_names",color="Rank3")






##### Euks N #####

#for lablefiles: column labels is for bargraphs, group is for a simple network, group2 is for photosynthetic/nonphotosynthetic network

#labelsEukN<-data.frame(rep("Mesofauna",dim(tax_table(datEukN3))[1]))

# ind<-which(is.na(labelsEukN))
# labelsEukN[ind]<-tax_table(datEukN3)[ind,"Rank5"]

unique(tax_table(datEukN5)[,"Rank4"])

labelsEukN<-substring(tax_table(datEukN5)[,"Rank4"],3)
ind<-which(is.na(labelsEukN));ind
labelsEukN[ind]<-"UnknownMetazoa"
 
colnames(labelsEukN)<-"labels"

labelsEukN2<-as.data.frame(labelsEukN)
unique(labelsEukN2$labels)
labelsEukN2$group<-"Mesofauna"

labelsEukN2$group2<-"Mesofauna"

head(labelsEukN2)

#9-13 are all NAs
unique(tax_table(datEukN5)[,8])
labelsEukN2$taxstring<-paste(tax_table(datEukN5)[,1],tax_table(datEukN5)[,2],tax_table(datEukN5)[,3],tax_table(datEukN5)[,4],tax_table(datEukN5)[,5],tax_table(datEukN5)[,6],tax_table(datEukN5)[,7],tax_table(datEukN5)[,8],sep=";")

#replace tax table
#tax_table(datEukN3)<-cbind(tax_table(datEukN3),labelsEukN)




##### Bacteria #####
#labelsBac<-data.frame(rep("Bacteria",dim(tax_table(datBacS3))[1]))
#I check the dataset for chemoautotrophs, no methanogens, yes ammonia oxidizers, 
#https://en.wikipedia.org/wiki/Microbial_metabolism#Sulfur_oxidation)
#https://link.springer.com/referenceworkentry/10.1007%2F978-3-642-38954-2_126

unique(tax_table(datBacS5)[,"Rank2"])

labelsBac<-substring(tax_table(datBacS5)[,"Rank2"],4)
sort(unique(labelsBac))

ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__");ind
labelsBac[ind]<-NA
ind<-which(is.na(labelsBac));ind
labelsBac[ind]<-"Unknown_Bacteria"

#Proteobacteria
#no info on o__Spirobacillales
ind<-which(labelsBac=="Proteobacteria");ind
labelsBac[ind]<-"Unknown_Proteobacteria"
#f__Nitrosomonadaceae (in Proteobacteria) from the Prokaryotes book "all of whose cultivated representatives are lithoautotrophic ammonia oxidizers". this is the only family in the nitrosomonadales order and everything I have is classified to genus, so no need to make any unknowns 
ind<-which(tax_table(datBacS5)[,"Rank5"]=="f__Nitrosomonadaceae");ind
labelsBac[ind]<-"Chemoautotrophic_Proteobacteria"
#The bradyrhizobiaceae (in proteobacteria) is complicated. the only species listed are bradyrhizobium (heterotroph that fixes N, although in papers some strains can be chemoautotrophs), Bosea (chemolithoheterotrophic) and another species that I dont know what it does (Balneimonas) or there is no genus/species listed. Some in this family however are photosynthetic and chemolithoautotrophs (nitrobacter which I don't have listed here). I will leave Bradyrhizobium and Bosea as heterotrophs and the rest as unknown
ind<-which(tax_table(datBacS5)[,"Rank5"]=="f__Bradyrhizobiaceae"&tax_table(datBacS5)[,"Rank6"]=="g__Bradyrhizobium");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
ind<-which(tax_table(datBacS5)[,"Rank5"]=="f__Bradyrhizobiaceae"&tax_table(datBacS5)[,"Rank6"]=="g__Bosea");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Rhodobiaceae the only genus we have afifella is heterotrophic
ind<-which(tax_table(datBacS5)[,"Rank6"]=="g__Afifella");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#f__Hyphomicrobiaceae unknown
#Rhizobiacea as far as I can tell are heterotrophs, Phyllobacteriaceae all heterotrophs except one that is facultatively chemolithoautotrophic so i will keep has heterotph, Methylocystaceae are methanotrophs use methane and methanol I will call them heterotrophs; Beijerinckiaceae on taxon is facultatively methanotrophic autotrophic so I will keep as a heterotroph; Xanthobacteraceae is chemoheterotrophs but many species are also facultatively chemolithoautotrophic soI will keep as heterotroph
ind<-which(tax_table(datBacS5)[,"Rank5"]%in%c("f__Rhizobiaceae","f__Phyllobacteriaceae","f__Methylocystaceae","f__Xanthobacteraceae"));ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Aurantimonadaceae too little is known
#Mycococcales are motile heterotrophs
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Myxococcales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Rickettiales are intracellular pathogens and maybe mutualists, so heterotrophs
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Rickettsiales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Syntrophobacterales, only one thing here in syntrophobacteraceae, most are heterotrophic but some are chemoautotrophs, so leave it at unknown
#o__Desulfuromonadales represented here by one genus geobacter, heterotroph
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Desulfuromonadales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Alteromonadales: main family is heterotrophic, two other unknown families there, but i will keep them all as heterotrophs
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Alteromonadales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Legionellales: two families, legionellaceae heterotrophs diseases, coxiellaceae also parasites
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Legionellales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Pseudomonadales, f__Moraxellaceae heterotrophs, pseudomonadaceae heterotrophs
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Pseudomonadales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Chromatilaes - (all taxa just to order) most of these are the photolithoautotrophic purple sulfur bacteria, however there is one family that is heterotrophic, so I'll leave at unknown
#Thiotrichales, one family Piscirickettsiaceae, both heterotrophs and chemolithoautotrophs
#Xanthomonadales: Sinobacteriaceae are heterotrophs, can't find much on Xanthomonadaceae but what i can find say they are heterotrophs
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Xanthomonadales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Enterobacteriales all appear to be heterotrophs
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Enterobacteriales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Procabacteriales not in the Prokaryotes book but looks like it is a heterotroph from google
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Procabacteriales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#burkholderiales are photoheterotrophs and other things but not autotrophs. one species is autotrophic (or can be maintained autotrophically in the lab) but I will keep it as a heterotroph group
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Burkholderiales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#o__Neisseriales, not too much info but mostly diseaes or microbiomes
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Neisseriales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Methylophilales, they are methylotrophs, I'll call them heterotrophs, many species soly use methanol CH3OH as the sources of carbon and entergy and some species use other small organic carbon compounds in addition to methanol.
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Methylophilales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Thiobacteraceae - unknown, cant find much info, might be a sulfur oxidizer but I don't knwo what that means
#Rhodocyclales - diverse metabolism, gemnus Uliginosibacterium and Propionivibrio and Sterolibacterium are heterotrophs, those with no genus leave as unknown 
ind<-which(tax_table(datBacS5)[,"Rank6"]%in%c("g__Uliginosibacterium","g__Propionivibrio","g__Sterolibacterium"));ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#o__Gallionellales chemolithoautotrophic
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Gallionellales");ind
labelsBac[ind]<-"Chemoautotrophic_Proteobacteria"
#[Entotheonellales] a sponge symbiont, it is an autotroph but also looks to have the capacity to be a heterotroph thus it is a mixotroph. leave it as unknown.
#The Rhodospirillaceae are purple non-sulfur bacteria but are hetero and autotrophs, alhtough g__Skermanella","g__Telmatospirillum are heteros. Acetobacteraceae all I can find is that they are heterotrophs, g__Roseomonas","g__Acidisphaera","g__Acidiphilium are all heterotrophs so I will jsut call them all heterotrophs
ind<-which(tax_table(datBacS5)[,"Rank6"]%in%c("g__Skermanella","g__Telmatospirillum"));ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
ind<-which(tax_table(datBacS5)[,"Rank5"]=="f__Acetobacteraceae");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#Sphingomonadales, heterotrophs
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Sphingomonadales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#o__Rhodobacterales" "f__Hyphomonadaceae" are heterotrophs, f__Rhodobacteraceae is unknown
ind<-which(tax_table(datBacS5)[,"Rank5"]=="f__Hyphomonadaceae");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
#caulobacterales chemoorganotrophic (as are many above, I assume this means heterotrophic but there might be some that are autotrophic in the archaea)
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Caulobacterales");ind
labelsBac[ind]<-"Heterotrophic_Proteobacteria"

#Spirochaetes are all heterotrophic

#OP3 are not well known, one is a chemolithoautotrophic but probably can't generalize

#p__Nitrospirae in bacteria. in Nitrospiraceae one genus nitrospira is chemoautotrophic, leptospirillacae is chemolithoautotrophic, but Thermodesulfovibrio heterotroph or a hydrogenotrophic sulfate reducer which I can't tell whether it is heterotroph or autotroph (might be autotroph, or at least I thin I read that some hydrogenotrophic sulfate reducers also have genes that fix co2), 0319-6A21 is unknown. and nitrospiracae is technically the only family so there are both autotrophs and heteotrophs within
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Nitrospirae");ind
labelsBac[ind]<-"Unknown_Nitrospirae"
ind<-which(tax_table(datBacS5)[,"Rank6"]=="g__Nitrospira"|tax_table(datBacS5)[,"Rank5"]=="f__[Leptospirillaceae]");ind
labelsBac[ind]<-"Chemoautotrophic_Nitrospirae"

#Verrucomicrobia, all except Methylacidiphilae are heterotrophic as far as I can tell, Methylacidiphilae are methanogens that are autotrophic and fix co2 (at least the genus that gives its name to the class is)
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Verrucomicrobia"&is.na(tax_table(datBacS5)[,"Rank3"]));ind
labelsBac[ind]<-"Unknown_Verrucomicrobia"
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Verrucomicrobia"&tax_table(datBacS5)[,"Rank3"]=="c__");ind
labelsBac[ind]<-"Unknown_Verrucomicrobia"
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Verrucomicrobia"&tax_table(datBacS5)[,"Rank3"]=="c__[Methylacidiphilae]");ind
labelsBac[ind]<-"Chemoautotrophic_Verrucomicrobia"
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Verrucomicrobia"&tax_table(datBacS5)[,"Rank3"]%in%c("c__[Spartobacteria]","c__Opitutae","c__[Pedosphaerae]","c__Verrucomicrobiae"));ind
labelsBac[ind]<-"Heterotrophic_Verrucomicrobia"

#Acidobacteria - all but the c__[Chloracidobacteria] are heterotrophs,Chloracidobacteria is photosynthetic
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Acidobacteria"&tax_table(datBacS5)[,"Rank3"]=="c__");ind
labelsBac[ind]<-"Unknown_Acidobacteria"
ind<-which(tax_table(datBacS5)[,"Rank3"]=="c__[Chloracidobacteria]");ind
labelsBac[ind]<-"Photoautotrophic_Acidobacteria"
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Acidobacteria"&tax_table(datBacS5)[,"Rank3"]!="c__[Chloracidobacteria]"&tax_table(datBacS5)[,"Rank3"]!="c__");ind
labelsBac[ind]<-"Heterotrophic_Acidobacteria"

#Elusimicrobia - heterotrophs

#Chloroflexi: Only the class Chloroflexi are photosynthetic, the other classes (Ktedontobacteria) are not photosynthetic. HOwever, there is one known recently discovered nitrifying bacterium (chemolithoautotroph) in the Thermomicrobia discovered in a reactor. but this is only one taxon not likely to be in our dataset so I am calling all other phyla heterotrophs.
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Chloroflexi"&tax_table(datBacS5)[,"Rank3"]=="c__");ind
labelsBac[ind]<-"Unknown_Chloroflexi"
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Chloroflexi"&is.na(tax_table(datBacS5)[,"Rank3"])==T);ind
labelsBac[ind]<-"Unknown_Chloroflexi"
ind<-which(tax_table(datBacS5)[,"Rank3"]=="c__Chloroflexi");ind
labelsBac[ind]<-"Photoautotrophic_Chloroflexi"
ind<-which(labelsBac=="Chloroflexi");ind
labelsBac[ind]<-"Heterotrophic_Chloroflexi"

#Armatimonadetes heterotrohpic
#p__[Thermi] hetrotrophic
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__[Thermi]");ind
labelsBac[ind]<-"Thermi"

#Gemmatimonadetes heterotrophic
#Fibrobacteres heterotrophic

#Actinobacteria a few can be chemoautotrophs
#groups with autotrophy and heterotrophy:
#o__Acidimicrobiales
#f__Nocardiaceae in the actinobactria (class and order) one species is autotrophic with hydrogen rhodococcus opacus, but most are saptroptrophs or pathogens, and four of my taxa match genus rhodococcus
#[c__Nitriliruptoria] not in dataset, [#o__Pseudonocardiales] not in dataset
#Start out with everything as a heterotroph 
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Actinobacteria");ind
labelsBac[ind]<-"Heterotrophic_Actinobacteria"
#label all organisms with no class as unknown
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Actinobacteria"&is.na(tax_table(datBacS5)[,"Rank3"]));ind
labelsBac[ind]<-"Unknown_Actinobacteria"
#Acidimicrobiales (only order in the class) is both hetero and auto, often the same organism can do both
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Acidimicrobiales");ind
labelsBac[ind]<-"Unknown_Actinobacteria"
#make nocardiaceae unknown and any actinobacteria class with no order unknown and any o__Actinomycetales order with no family unknown
ind<-which(tax_table(datBacS5)[,"Rank5"]=="f__Nocardiaceae");ind
labelsBac[ind]<-"Unknown_Actinobacteria"
ind<-which(tax_table(datBacS5)[,"Rank3"]=="c__Actinobacteria"&is.na(tax_table(datBacS5)[,"Rank4"]));ind
labelsBac[ind]<-"Unknown_Actinobacteria"
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Actinomycetales"&is.na(tax_table(datBacS5)[,"Rank5"]));ind
labelsBac[ind]<-"Unknown_Actinobacteria"
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Actinomycetales"&tax_table(datBacS5)[,"Rank5"]=="f__");ind
labelsBac[ind]<-"Unknown_Actinobacteria"

#BRC1 heterotrophic
#GAL15 unknown
#AD3 unknown
#NKB19 unknown
#FCPU426 unknown
#WS3 heterotrophic

#Planctomycetes
#The Planctomycetes are also complicated. there are five genera that do anammox (chemoautotrophic), all in the Brocadiales, none of that order or those genera are found in this dataset, however a few of the order slots are blanks. There was some confusion about a deeply branching group which is separate from the anammox group - a number of papers helped me clarify this (http://aem.asm.org/content/73/15/4707.full, https://aem.asm.org/content/67/2/623?ijkey=39ef9714e7624359db10d4584f438d28a65adb76&keytype2=tf_ipsecsha, https://aem.asm.org/content/69/12/7354?ijkey=b33ae38bdab6244d3695d980c8ce262530de592b&keytype2=tf_ipsecsha). Not much is known about the deeply branching group and one paper says that thy might hold unkonwn metabolic diversity, but as of now I will call them all heterotrophs, since they clearly do not cluster with brocadiales. also according to that paper and the prokaryotes book the genera in the o__Gemmatales and Pirellulales orders are all heterotrophs and genus Planctomyces is heterotrophic. [checking one other class OM190 is heterotrophic, https://www.frontiersin.org/articles/10.3389/fmicb.2014.00267/full] Overall, Brocadiales (in class Planctomycetia) are the only known anammox chemoautotrophs, so I will put anything without a class as unknown. All the planctomycetia have orders that are not brocadiales so I will call them all heterotrophs. everything else I will call heterotrophs

ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Planctomycetes");ind
labelsBac[ind]<-"Heterotrophic_Planctomycetes"
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Planctomycetes"&is.na(tax_table(datBacS5)[,"Rank3"]));ind
labelsBac[ind]<-"Unknown_Planctomycetes"
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Planctomycetes"&tax_table(datBacS5)[,"Rank3"]=="c__");ind
labelsBac[ind]<-"Unknown_Planctomycetes"

#BHI80-139 unknown
#p__Chlamydiae heterotrophic

#Firmicutes 
#some families in firmicutes are chemolithoautotrophs (or the taxa can do chemolithoautotrophy as well as heterotrophy). The class Bacilli has chemolithoautotrophy, as does some members of: Alicyclobacillaceae, Halanaerobiales, Heliobacteriaceae, Peptococcaceae, Desulfitobacteriaceae,Desulfotomaculaceae, Thermincolaceae,Thermodesulfobacteriaceae,Thermoanaerobacteraceae, Thermolithobacteriaceae. Bacillus - there may be one taxon that is capable of chemolithoautotrohpy but in general they are heterotrophs so I will leave them as such. Syntrophomonadaceae is often not a "real" autotroph b/c often depends on acetate (as is true of many H2 utilizing methanogens or sulfate reducers)
#we have Desulfosporosinus which is capable of chemolithoautotrophic growth
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Firmicutes");ind
labelsBac[ind]<-"Heterotrophic_Firmicutes"
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Firmicutes"&tax_table(datBacS5)[,"Rank6"]=="g__Desulfosporosinus");ind
labelsBac[ind]<-"Chemoautotrophic_Firmicutes"
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Firmicutes"&is.na(tax_table(datBacS5)[,"Rank5"]));ind
labelsBac[ind]<-"Unknown_Firmicutes"

#WS2 unknown
#WPS-2 unknown
#OD1 (or Parcubacteria) unknown. some are heterotrophs (The reduced genomes of Parcubacteria (OD1) contain signatures of a symbiotic lifestyle) but some can fix CO2 (Insights into the phylogeny and coding potential of microbial dark matter)
#SR1 unknown
#WS5 unknown
#TM7 (aka Saccharibacteria) definitely heterotrophic but found in weird reactors so could be both hetero and autotrophic, not much known, call unknown
#OP11 heterotrophic Partial genome assembly for a candidate division OP11 single cell from an anoxic spring (Zodletone Spring, Oklahoma)
# MVP-21 unknown
#Cyanobacteria ps
#GN02 (Gracilibacteria) unknown
#Tenericutes heterotrophs

##OD1, OD11, SR1, and GN02 are all related - are likely majority heterotrophs but so much is unknown and at least some OD1 can fix CO2 (Insights into the phylogeny and coding potential of microbial dark matter)

#Crenarchaeota in archaea: the o__Nitrososphaerales is an ammonia oxidizer, chemolithoautotroph fixes CO2 (all the abundant crenarchaeota are these), Nitrosotalea genus is also an chemolithoautotroph (Cultivation of an obligate acidophilic ammonia oxidizer from a nitrifying acid soil.) the other less abundant orders (like NRP-J and crenarcheales) I'm not sure about (I can't find information on them) but there are chemoorganotrophs in crenarchaeota so I can't assume they are autotrophs
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__Crenarchaeota");ind
labelsBac[ind]<-"Unknown_Crenarchaeota"
ind<-which(tax_table(datBacS5)[,"Rank4"]=="o__Nitrososphaerales");ind
labelsBac[ind]<-"Chemoautotrophic_Crenarchaeota"
ind<-which(tax_table(datBacS5)[,"Rank6"]=="g__Nitrosotalea");ind
labelsBac[ind]<-"Chemoautotrophic_Crenarchaeota"

#Euryarchaeota. all in class thermoplasmata and order E2 which is now the Methanoregulaceae. methanogens are confusing, some are autotrophs and do CO2 + 4 H2 → CH4 + 2 H2O, some are heterotrophs the most common relying on acetate CH3COOH → CH4 + CO2. autotrophic methanogens have genes for acetyl-CoA synthesis from CO2 (acetyl-CoA is the beginning of the krebs cycle, thus if you can make acetyl coa from co2 you are an autotroph, typically acetyl coa is from other organic compounds like glucose or at least acetate which is organic). this family in the prokaryotes says acetate is required for growht so I will call them heterotrophs

#Parvarchaeota - unknown
ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__[Parvarchaeota]");ind
labelsBac[ind]<-"Parvarchaeota"

#Bacteroidetes - heterotrophs
#Chlorobi - many are green sulfur bacteria (photolithoautotrophs), but now it also includes a few photoheterotrophic taxa. i will keep it as unknown
#FBP - unknown
#TM6 - many are symbionts (heterotroph),Comparative Genomics of Candidate Phylum TM6 Suggests That Parasitism Is Widespread and Ancestral in This Lineage - go with heterotrophs

ind<-which(tax_table(datBacS5)[,"Rank2"]=="p__TM6");ind
unique(as(tax_table(datBacS5),"matrix")[ind,])  #something weird happened on my work computer and unique() did not work on phyloseq objects
unique(tax_table(datBacS5)[ind,])

sort(unique(labelsBac))
colnames(labelsBac)<-"labels"

labelsBac2<-as.data.frame(labelsBac)
labelsBac2$group<-"Bacteria"

labelsBac2$group2<-NA
ind<-which(labelsBac2$labels=="Photoautotrophic_Acidobacteria"|labelsBac2$labels=="Cyanobacteria"|labelsBac2$labels=="Photoautotrophic_Chloroflexi") 
labelsBac2$group2[ind]<-"PhotosyntheticBacteria"
ind<-which(labelsBac2$labels%in%c("Heterotrophic_Acidobacteria","Heterotrophic_Actinobacteria","Heterotrophic_Chloroflexi","Heterotrophic_Firmicutes","Heterotrophic_Planctomycetes","Heterotrophic_Proteobacteria","Heterotrophic_Verrucomicrobia","Spirochaetes","Elusimicrobia","Armatimonadetes","Thermi","Gemmatimonadetes","Fibrobacteres","BRC1","WS3","Chlamydiae","OP11","Tenericutes","Euryarchaeota","Bacteroidetes","TM6"))
labelsBac2$group2[ind]<-"HeterotrophicBacteria"
ind<-which(labelsBac2$labels=="Chemoautotrophic_Crenarchaeota"|labelsBac2$labels=="Chemoautotrophic_Firmicutes"|labelsBac2$labels=="Chemoautotrophic_Proteobacteria"|labelsBac2$labels=="Chemoautotrophic_Nitrospirae"|labelsBac2$labels=="Chemoautotrophic_Verrucomicrobia") 
labelsBac2$group2[ind]<-"ChemoautotrophicBacteria"
ind<-which(is.na(labelsBac2$group2)==T)
labelsBac2$group2[ind]<-"UnknownBacteria"
head(labelsBac2)

#replace tax table
#tax_table(datBacS3)<-cbind(tax_table(datBacS3),labelsBac)

dim(tax_table(datBacS5))
unique(tax_table(datBacS5)[,7])
labelsBac2$taxstring<-paste(tax_table(datBacS5)[,1],tax_table(datBacS5)[,2],tax_table(datBacS5)[,3],tax_table(datBacS5)[,4],tax_table(datBacS5)[,5],tax_table(datBacS5)[,6],tax_table(datBacS5)[,7],sep=";")


##### Fungi #####

#labelsITS<-data.frame(rep("Fungi",dim(tax_table(datITSS3))[1]))

unique(tax_table(datITSS5)[,"Rank2"])

labelsITS<-substring(tax_table(datITSS5)[,"Rank2"],4)
ind<-which(is.na(labelsITS)==T);ind
labelsITS[ind]<-"Unknown_Fungi"
unique(labelsITS)
colnames(labelsITS)<-"labels"

labelsITS2<-as.data.frame(labelsITS)
#labelsITS2<-as.data.frame(as(labelsITS,"matrix"))
labelsITS2$group<-"Fungi"
labelsITS2$group2<-"Fungi"

print(labelsITS2[1:5,],row.names=F)

#tax_table(datBacS3)<-cbind(tax_table(datBacS3),labelsBac)

unique(tax_table(datITSS5)[,7])
labelsITS2$taxstring<-paste(tax_table(datITSS5)[,1],tax_table(datITSS5)[,2],tax_table(datITSS5)[,3],tax_table(datITSS5)[,4],tax_table(datITSS5)[,5],tax_table(datITSS5)[,6],tax_table(datITSS5)[,7],sep=";")

#replace tax table
#tax_table(datITSS3)<-cbind(tax_table(datITSS3),labelsITS)




##### Make otu tables (bact takes ~10 min) #####
datEukS5otu<-cbind(sample_data(datEukS5),t(otu_table(datEukS5)))
datEukS5otu$Sample_name<-as.numeric(as.character(datEukS5otu$Sample_name))

datEukN5otu<-cbind(sample_data(datEukN5),t(otu_table(datEukN5)))
datEukN5otu$Sample_name<-as.numeric(as.character(datEukN5otu$Sample_name))

datBacS5otu<-cbind(sample_data(datBacS5),t(otu_table(datBacS5)))
datBacS5otu$Sample_name<-as.numeric(as.character(datBacS5otu$Sample_name))

datITSS5otu<-cbind(sample_data(datITSS5),otu_table(datITSS5))
datITSS5otu$Sample_name<-as.numeric(as.character(datITSS5otu$Sample_name))

#otu for counts (not relative abundance)
datEukS5cotu<-cbind(sample_data(datEukS5c),t(otu_table(datEukS5c)))
datEukS5cotu$Sample_name<-as.numeric(as.character(datEukS5cotu$Sample_name))

datEukN5cotu<-cbind(sample_data(datEukN5c),t(otu_table(datEukN5c)))
datEukN5cotu$Sample_name<-as.numeric(as.character(datEukN5cotu$Sample_name))

datBacS5cotu<-cbind(sample_data(datBacS5c),t(otu_table(datBacS5c)))
datBacS5cotu$Sample_name<-as.numeric(as.character(datBacS5cotu$Sample_name))

datITSS5cotu<-cbind(sample_data(datITSS5c),otu_table(datITSS5c))
datITSS5cotu$Sample_name<-as.numeric(as.character(datITSS5cotu$Sample_name))





##### Phylogenetic diversity ######
# This needs to be done before any label changes (labels of taxa) below b/c it uses the tree data.
# I need to use root=T, if I use root=F it cannot calculate the diversity in a sample with only one taxon
pdEukS<-pd(as.matrix(datEukS5cotu[,-c(1:33)]),phy_tree(datEukS5c),include.root=TRUE) #took 5 minutes
pdEukN<-pd(as.matrix(datEukN5cotu[,-c(1:33)]),phy_tree(datEukN5c),include.root=TRUE) #took 3 minutes
pdBac<-pd(as.matrix(datBacS5cotu[,-c(1:33)]),phy_tree(datBacS5c),include.root=TRUE) #used to take an hour, now took 3 min
richITS1<-data.frame(PD=rowSums(datITSS5cotu[,-c(1:33)]>0),SR=rowSums(datITSS5cotu[,-c(1:33)]>0))


#add chao1 richness to all dfs
richEukS<-estimate_richness(datEukS5c, split = TRUE, measures = c("Observed", "Shannon", "Chao1","se.chao1"))
richEukN<-estimate_richness(datEukN5c, split = TRUE, measures = c("Observed", "Shannon", "Chao1","se.chao1"))
richBac<-estimate_richness(datBacS5c, split = TRUE, measures = c("Observed", "Shannon", "Chao1","se.chao1"))
richITS<-estimate_richness(datITSS5c, split = TRUE, measures = c("Observed", "Shannon", "Chao1","se.chao1"))

#merge the above two dfs
richEukS2<-cbind(richEukS,pdEukS)
richEukN2<-cbind(richEukN,pdEukN)
richBac2<-cbind(richBac,pdBac)
richITS2<-cbind(richITS,richITS1)

#Just to test
#richBac<-as.data.frame(rowSums(datBacS5cotu[,-c(1:33)]>0))



###### Grouping by kingdom/phylum #####
#for the bar graphs

datEukS5k<-aggregate.data.frame(otu_table(datEukS5),by=list(labels=labelsEukS2$labels),sum)
rownames(datEukS5k)<-datEukS5k$labels
datEukS5k$labels<-NULL
datEukS5k2<-cbind(sample_data(datEukS5),t(datEukS5k))
head(datEukS5k2)

datEukN5k<-aggregate.data.frame(otu_table(datEukN5),by=list(labels=labelsEukN2$labels),sum)
rownames(datEukN5k)<-datEukN5k$labels
datEukN5k$labels<-NULL
datEukN5k2<-cbind(sample_data(datEukN5),t(datEukN5k))
head(datEukN5k2)

datBacS5k<-aggregate.data.frame(otu_table(datBacS5),by=list(labels=labelsBac2$labels),sum)
rownames(datBacS5k)<-datBacS5k$labels
datBacS5k$labels<-NULL
datBacS5k2<-cbind(sample_data(datBacS5),t(datBacS5k))
head(datBacS5k2)

datITSS5k<-aggregate.data.frame(t(otu_table(datITSS5)),by=list(labels=labelsITS2$labels),sum)
rownames(datITSS5k)<-datITSS5k$labels
datITSS5k$labels<-NULL
datITSS5k2<-cbind(sample_data(datITSS5),t(datITSS5k))
head(datITSS5k2)
  
  
  
  


###### Filter data sets for network analysis ######
#Follwing Widder et al 2014 PNAS

#take out doubletons and singletons. this doesn't really matter b/c I will take out taxa with 7 or fewer occurrences later on, Im just doing it now to make the dataframe a little smaller and more manageable
ind<-(which(colSums(datEukS5otu[,34:dim(datEukS5otu)[2]]>0)>2))+33
datEukS5otu2<-cbind(datEukS5otu[,1:33],datEukS5otu[,ind])
datEukS5cotu2<-cbind(datEukS5cotu[,1:33],datEukS5cotu[,ind])

ind<-(which(colSums(datEukN5otu[,34:dim(datEukN5otu)[2]]>0)>2))+33
datEukN5otu2<-cbind(datEukN5otu[,1:33],datEukN5otu[,ind])
datEukN5cotu2<-cbind(datEukN5cotu[,1:33],datEukN5cotu[,ind])

ind<-(which(colSums(datBacS5otu[,34:dim(datBacS5otu)[2]]>0)>2))+33
datBacS5otu2<-cbind(datBacS5otu[,1:33],datBacS5otu[,ind])
datBacS5cotu2<-cbind(datBacS5cotu[,1:33],datBacS5cotu[,ind])

ind<-(which(colSums(datITSS5otu[,34:dim(datITSS5otu)[2]]>0)>2))+33
datITSS5otu2<-cbind(datITSS5otu[,1:33],datITSS5otu[,ind])
datITSS5cotu2<-cbind(datITSS5cotu[,1:33],datITSS5cotu[,ind])

#filter out taxa that have a summed relative abundance of <.002 (.2%)
#I tested this and found that there were about 35 bacteria who had 8 or more occurrences but still a summed rel abun of <.002, so this step does remove some taxa from the network independent of the occurrence removal
ind<-(which(colSums(datEukS5otu2[,34:dim(datEukS5otu2)[2]])>0.002))+33
datEukS5otu3<-cbind(datEukS5otu2[,1:33],datEukS5otu2[,ind]) #
datEukS5cotu3<-cbind(datEukS5cotu2[,1:33],datEukS5cotu2[,ind]) #1091 otu

ind<-(which(colSums(datEukN5otu2[,34:dim(datEukN5otu2)[2]])>0.002))+33
datEukN5otu3<-cbind(datEukN5otu2[,1:33],datEukN5otu2[,ind]) #142 otu
datEukN5cotu3<-cbind(datEukN5cotu2[,1:33],datEukN5cotu2[,ind]) #

ind<-(which(colSums(datBacS5otu2[,34:dim(datBacS5otu2)[2]])>0.002))+33
datBacS5otu3<-cbind(datBacS5otu2[,1:33],datBacS5otu2[,ind]) #4022 otu
datBacS5cotu3<-cbind(datBacS5cotu2[,1:33],datBacS5cotu2[,ind]) #

ind<-(which(colSums(datITSS5otu2[,34:dim(datITSS5otu2)[2]])>0.002))+33
datITSS5otu3<-cbind(datITSS5otu2[,1:33],datITSS5otu2[,ind]) #895 otu
datITSS5cotu3<-cbind(datITSS5cotu2[,1:33],datITSS5cotu2[,ind]) #

#order of doing things: filtered out unwanted taxa (chloroplasts, spiders) and samples (S.2015), rarefied, relativized, took out doubletons and singletons, took out samples <2% summed abundance. I am not going to rarefy or relativize here again b/c the doubletons/singletons/.2% otus that I removed are real, I could have included them in the network analyis if I wanted.


##### Plants #####

plantcomp<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Plants/Niwot_MovingUpHill_comp2015.csv")
head(plantcomp)
names(plantcomp)[1]<-"Sample_name"

#Remove plant species only present in one or two plots; there are some plots that have plant data but not microbe data. 69 70 71 77 81 108 117 118 147 148 149 151. This is because when we started doing the surveys we were going to all plots for plants and only sample some for microbes, then we realized that that was insane!
dim(plantcomp)
plantcomp2<-plantcomp[,colSums(plantcomp>0)>2]
plantcomp2$LICHEN<-NULL  #54 plants


###### Make labelfile ######

#Test if there is name overlap in full otu tables, no names overlap, no need to rename OTU columnson the full otu table files, however it is nice for figuring out what is going on in the modeling, so I will add B, N, S, I before all the otu names for the reduced otu table files
namesEukS<-names(d[,-c(1:33)])
namesEukN<-names(datEukN5otu[,-c(1:33)])
namesBac<-names(datBacS5otu[,-c(1:33)])
namesITS<-names(datITSS5otu[,-c(1:33)])
length(namesBac)+length(namesEukS)
length(union(namesBac,namesEukS))
intersect(namesEukS,namesEukN)
intersect(namesEukS,namesBac)
intersect(namesEukS,namesITS)
intersect(namesEukN,namesBac)
intersect(namesEukN,namesITS)
intersect(namesBac,namesITS)

# namesEukS2 <- sub("^", "S", namesEukS)
# namesEukN2 <- sub("^", "N", namesEukN)
# namesBac2 <- sub("^", "B", namesBac)
# namesITS2 <- sub("^", "I", namesITS)
# 
# names(datEukS3otu)[-c(1:31)]<-namesEukS2
# names(datEukN3otu)[-c(1:31)]<-namesEukN2
# names(datBacS3otu)[-c(1:31)]<-namesBac2
# names(datITSS3otu)[-c(1:31)]<-namesITS2


#Change colnames on the reduced datasets for the networks
#for widder et al reduced datasets
namesEukS<-names(datEukS5otu3[,-c(1:33)])
namesEukN<-names(datEukN5otu3[,-c(1:33)])
namesBac<-names(datBacS5otu3[,-c(1:33)])
namesITS<-names(datITSS5otu3[,-c(1:33)])

namesEukS2 <- sub("^", "S", namesEukS)
namesEukN2 <- sub("^", "N", namesEukN)
namesBac2 <- sub("^", "B", namesBac)
namesITS2 <- sub("^", "I", namesITS)

names(datEukS5otu3)[-c(1:33)]<-namesEukS2
names(datEukS5cotu3)[-c(1:33)]<-namesEukS2
names(datEukN5otu3)[-c(1:33)]<-namesEukN2
names(datEukN5cotu3)[-c(1:33)]<-namesEukN2
names(datBacS5otu3)[-c(1:33)]<-namesBac2
names(datBacS5cotu3)[-c(1:33)]<-namesBac2
names(datITSS5otu3)[-c(1:33)]<-namesITS2
names(datITSS5cotu3)[-c(1:33)]<-namesITS2

#for full datasets
namesEukS<-names(datEukS5otu[,-c(1:33)])
namesEukN<-names(datEukN5otu[,-c(1:33)])
namesBac<-names(datBacS5otu[,-c(1:33)])
namesITS<-names(datITSS5otu[,-c(1:33)])

namesEukS2 <- sub("^", "S", namesEukS)
namesEukN2 <- sub("^", "N", namesEukN)
namesBac2 <- sub("^", "B", namesBac)
namesITS2 <- sub("^", "I", namesITS)

names(datEukS5otu)[-c(1:33)]<-namesEukS2
names(datEukS5cotu)[-c(1:33)]<-namesEukS2
names(datEukN5otu)[-c(1:33)]<-namesEukN2
names(datEukN5cotu)[-c(1:33)]<-namesEukN2
names(datBacS5otu)[-c(1:33)]<-namesBac2
names(datBacS5cotu)[-c(1:33)]<-namesBac2
names(datITSS5otu)[-c(1:33)]<-namesITS2
names(datITSS5cotu)[-c(1:33)]<-namesITS2





#Make combined labelfile
labelsEukS2$otu<-sub("^", "S",rownames(labelsEukS2))
labelsEukN2$otu<-sub("^", "N",rownames(labelsEukN2))
labelsBac2$otu<-sub("^", "B",rownames(labelsBac2))
labelsITS2$otu<-sub("^", "I",rownames(labelsITS2))

labelfile1<-rbind(labelsEukS2,labelsEukN2,labelsBac2,labelsITS2)
head(labelfile1)
labelfile1$oldotu<-rownames(labelfile1)
#labelfile1$oldotu<-substring(rownames(labelfile1), 2)

#combine with plant labelfile
labelsPlant<-as.data.frame(cbind(labels="Plant",group="Plant",group2="Plant",taxstring=colnames(plantcomp2)[2:55],otu=colnames(plantcomp2)[2:55],oldotu=colnames(plantcomp2)[2:55]))
head(labelsPlant)

labelfile<-rbind(labelfile1,labelsPlant)
head(labelfile)
tail(labelfile)


#####Output files####
labelfile

#richness/diverity
richEukS2
richEukN2
richBac2
richITS2

#full rarefied count otu tables (for ordination and for networks)
datEukS5otu
datEukS5cotu
datEukN5otu
datEukN5cotu
datBacS5otu
datBacS5cotu
datITSS5otu
datITSS5cotu

#reduced rarefied rel abun and count otu tables (originally for networks, but I didn't actually use this, I used above)
datEukS5otu3
datEukS5cotu3
datEukN5otu3
datEukN5cotu3
datBacS5otu3
datBacS5cotu3
datITSS5otu3
datITSS5cotu3

#grouped by kingdom (for bargraphs)
datEukS5k2
datEukN5k2
datBacS5k2
datITSS5k2

#plants
plantcomp
plantcomp2

#create an env with only the above objects
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
                        "plantcomp2")))

