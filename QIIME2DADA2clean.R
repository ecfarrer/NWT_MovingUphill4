##QIIME2

#notes: this is for installing an old version of qiime2, be sure to find the current release and modify the code accordingly

##First I needed to uninstall and then reinstall miniconda3. it would not update by itself with the update code
https://conda.io/docs/user-guide/install/macos.html

#The second time I tried (when reinstalling for the patch) it worked
conda update conda
conda install wget

##How to install qiime2
https://docs.qiime2.org/2017.12/install/native/

wget https://data.qiime2.org/distro/core/qiime2-2018.2-py35-osx-conda.yml
conda env create -n qiime2-2018.2 --file qiime2-2018.2-py35-osx-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2018.2-py35-osx-conda.yml

wget https://data.qiime2.org/distro/core/qiime2-2019.1-py36-osx-conda.yml
conda env create -n qiime2-2019.1 --file qiime2-2019.1-py36-osx-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2019.1-py36-osx-conda.yml



#How to uninstll qiime2
#I needed to do this to download the new patch for the ITS dada2 error
conda env remove -n qiime2-2017.12

#Activate Qiime env, must do this in every new terminal tab that is opened
source activate qiime2-2018.2 #(this was used for ITS demultiplexing and Euk everything)
source activate qiime2-2019.1 #(this was used for bacteria paired end)

#Test the new environment
qiime --help

#Paired end tutorial
https://docs.qiime2.org/2017.12/tutorials/atacama-soils/
#General tutorial
https://docs.qiime2.org/2017.12/tutorials/moving-pictures/


#discussion about whether it is worth it to join paired end reads when the reads completely overlap
https://forum.qiime2.org/t/question-about-dada2-denoise-paired-analysis/464
#antoher discussion of relaxing the maxEE filtering parameter when using full reads with poor quality scores near the ends (in qiime dada2 denoise-paired). however the default does delete everyting after the first instance of a bp with quality score of 2 (this can be changed as well)
https://github.com/qiime2/q2-dada2/issues/48






##### Taxonomic databases #####

##### Silva #####
#notes about the silva release: https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_128_notes.txt
#workflow from this thread: https://forum.qiime2.org/t/18s-classifier-using-silva-database-and-emb-primers/361

#With all taxa and the 111 release, start 1:31pm, end 1:33
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path Silva_111_post/rep_set/Silva_111_full_unique.fasta \
--output-path all_SILVA111_unique_otus.qza


#Import taxonomy, only with 18S.
# There are lots of options here, there is a file that says "no ambiguous" but I dont know if they mean, I assume no ambiguous taxa. there is also a file with consistent 6 ranks (the help said that RDP requires all taxa to have consistent number of ranks. I can try this if the full taxonomy doesn't work). there are also 99 clustered files, but I will use full database
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path Silva_111_post/taxonomy/Silva_111_taxa_map_full.tsv \
--output-path all_SILVA111_unique_taxonomy.qza


##### Extract EMBP variable region for 18S #####
#with single end data you should set a length. When I was playing around I used 150 b/c I was following the person writing the workflow. I could use the max length from my euk data below (which is generally less than 150, except there are about 2000 reads per sample that are 301)
#with paired end data (what I ended up useing), trunc-len should not be used (set at default 0). post saying that you shouldn't truncate with paired end data: https://forum.qiime2.org/t/how-can-i-train-classifier-for-paired-end-reads/1512/7, I think truncation is for when you have only a forward read and it is a particular length, start 1:32pm, end 2:14pm

#Without the 150 truncation, start 2:30pm, end 3:17pm
qiime feature-classifier extract-reads \
--i-sequences all_SILVA111_unique_otus.qza \
--p-f-primer GTACACACCGCCCGTC \
--p-r-primer TGATCCTTCTGCAGGTTCACCTAC \
--o-reads ref-seqs_all_unique_SILVA111.qza

#export to see what it trimmed. 
#Most lengths are less than 140, there are some longer though
qiime tools export \
ref-seqs_all_unique_SILVA111.qza \
--output-dir exported-ref-seqs_all_unique_SILVA111

awk '{print length}' dna-sequences.fasta | sort | uniq -c


#Train classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_all_unique_SILVA111.qza \
--i-reference-taxonomy all_SILVA111_unique_taxonomy.qza \
--o-classifier all_EMB_SILVA111_classifier.qza


#talking about what the different outputs (truncated taxonomy vs blank genus/species levels means: 
https://forum.qiime2.org/t/consensus-blast-taxonomy-strings/586/2




##### Greengenes database ######
#I redid this, need to if you update qiime
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path gg_13_8_otus/rep_set/99_otus.fasta \
--output-path gg_13_8_otus_99_otus.qza

#Import taxonomy
#Need to duplicate and save the txt file as tsv
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path gg_13_8_otus/taxonomy/99_otu_taxonomy.tsv \
--output-path gg_13_8_otus_99_taxonomy.qza

#start 12:28 pm, end 12:48 pm
#note: if you have paired reads, do not set --p-trunc-len. I added the min and max - my reads are between 235-394, this was used by the moving pictures tutorial which use the same primers
qiime feature-classifier extract-reads \
--i-sequences gg_13_8_otus_99_otus.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer GGACTACNVGGGTWTCTAAT \
--p-min-length 100 \
--p-max-length 400 \
--o-reads ref-seqs_all_99_gg_13_8.qza

#export to see what it trimmed
#qiime tools export \
#ref-seqs_all_99_gg_13_8.qza \
#--output-dir exported-ref-seqs_all_99_gg_13_8

#cd exported-ref-seqs_all_99_gg_13_8
#awk '{print length}' dna-sequences.fasta | sort | uniq -c
#Most lengths are 253, there are some longer and shorter from about 251-257

##Train classifier
#start 3:35pm, end 3:42pm, then move into Bact folder
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_all_99_gg_13_8.qza \
--i-reference-taxonomy gg_13_8_otus_99_taxonomy.qza \
--o-classifier all_EMB_gg_13_8_classifier.qza




##### UNITE database #####
#I downloaded the file 12_11 alpha release from the qiime page here: http://qiime.org/home_static/dataFiles.html I didn't see it on the unite website...
#tutorial https://github.com/gregcaporaso/2017.06.23-q2-fungal-tutorial
#I tried downloading the UNITE release and processing it but I kept getting errors, I think there were a bunch of formatting errors in the file, so I downloaded an already trained classifir from unite here (UNITE version 7.2): https://forum.qiime2.org/t/unite-ver-7-2-2017-12-01-classifiers-for-qiime2-ver-2017-12-available-here/3020
#She followed the emp protocol for training it

unite-ver7-99-classifier-01.12.2017.qza



##### Euk #####

##### Demultiplexing #####

#Import the data into an artifact
#Navigate to: /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figures\&Stats/kingdata/QIIME2/Euks

#first rename your files to (R1) forward.fastq, (R2) reverse.fastq, (I) and barcodes.fastq
#then gzip them, takes a couple minutes
gzip barcodes.fastq
gzip reverse.fastq 
gzip forward.fastq

#start 4:08pm, end 4:15
qiime tools import \
--type EMPPairedEndSequences \
--input-path emp-paired-end-sequences \
--output-path emp-paired-end-sequences.qza

#move mapping file into directory, and delete .txt and replace with .tsv

#barcode is on the forward read (see word doc showing how the R1 has the reverse primers in the read, while the R2 has the forward primers, this is b/c the sequence is so short that it kept reading the other half of the primer. However the mapping file is correct: "reverse primer" matches to the R2). looking at the barcodes and the barcode file, it looks like it is the reverse complement
#start 5:20pm, end 6:30
qiime demux emp-paired \
--m-barcodes-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--m-barcodes-category BarcodeSequence \
--i-seqs emp-paired-end-sequences.qza \
--o-per-sample-sequences demux \
--p-rev-comp-mapping-barcodes

#summarize, start 8:44pm, end 9:03pm
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

qiime tools view demux.qzv

#Export it so you can look inside at the files
#started 8:58pm, end 9:00pm
qiime tools export \
demux.qza \
--output-dir exported-demux
  
##### Trim primers/adapters #####
#To get these adapters I looked in the raw sequnce data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
#start 9:44pm, end 10:02
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux.qza \
--p-cores 4 \
--p-adapter-f GTAGGTGAACCTGCAGAAGGATCA \
--p-adapter-r GACGGGCGGTGTGTAC \
--o-trimmed-sequences trimmed-seqs.qza \
--verbose
#note that it said that for the reverse read, the primer was preceeded by C extremely often, which might mean it is part of the primer. it shouldn't be though, so I'm leaving it in.
#note I could have used --p-error-rate 0 \  meaning no errors in the adapter, the default used above is 10%, not sure what difference it would make

#summarize
qiime demux summarize \
--i-data trimmed-seqs.qza \
--o-visualization trimmed-seqs.qzv

qiime tools view trimmed-seqs.qzv

#export
qiime tools export \
trimmed-seqs.qza \
--output-dir exported-trimmed-seqs
#look at this tomorrow and make sure it trimmed things correctly!!! Yes they are, looks good.
awk '{print length}' N.0.2015_49_L001_R1_001copy.fastq | sort | uniq -c

# note that for some of the forward reads start with N for a basepair, this is actually "true" when you look at the reverse and forward reads, the N represents a certain bp that was apparently unknown in the sequencing. (also when it starts with a g rather than an n (gctac) the g blasts to something so the g is correct). The DADA2 tutorial https://benjjneb.github.io/dada2/tutorial.html states that no N are allowable in DADA2, so I need to get rid of that first basepair.

grep --color -n "^N"  N.0.2015_49_L001_R1_001copy.fastq # ^ means at the beginning of the line. there actually aren't that many, only like 30 reads start with N
grep --color -n "^N"  N.0.2015_49_L001_R2_001copy.fastq #not any Ns at the beginning of reads like in R1

##### Denoising with DADA2 ##### 
#220 and 210 respectively are when the median quality scores hit 25, however reading the documentation it says reads that are shorter than the truncation numbers are discarded (this seems silly but it is what it is), thus I should keep everything (I tried doing truncation prior to dada2 but it didn't work). Also, I don't think it matters as much for euks b/c the read is so short, if the primer is successfully removed, there is no issue with quality. reads are ~130bp, primer is 24 or 16bp, so the primers are well within the good quality read.
#note for n-threads specifying 0 means use all cores
#note for trunc-len specifying 0 means don't truncate
#the p-trim-left-f is set at 1 to get rid of the N basepair
#start 8:11pm, end 9:22
qiime dada2 denoise-paired \
--i-demultiplexed-seqs trimmed-seqs.qza \
--o-table table \
--o-representative-sequences rep-seqs \
--p-n-threads 6 \
--p-trim-left-f 1 \
--p-trim-left-r 0 \
--p-trunc-len-f 0 \
--p-trunc-len-r 0

#export rep seqs just to take a look
qiime tools export \
rep-seqs.qza \
--output-dir exported-rep-seqs

#to print the line with the longest number of characters:
cat filename|awk '{print length, $0}'|sort -nr|head -1
#the longest line in the rep set is 244bp, and blasts to a fungus

#export table to get otu table
qiime tools export \
table.qza \
--output-dir exported-table

#convert biom to otu table text file!
biom convert -i feature-table.biom -o otu_table.txt --to-tsv

#delete the space and # from '#OTU ID' on line 2 
sed '2s/[ ]//' otu_table.txt | sed '2s/.//' > otu_table2.txt

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv

qiime tools view table.qzv

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

qiime tools view rep-seqs.qzv


###### Create a phylogenetic tree #####

#do the alignment
#start 9:56pm, end  10:00. could add --p-n-threads to the code
qiime alignment mafft \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza

#Mask (or filter) the alignment to remove positions that are highly variable. These positions are generally considered to add noise to a resulting phylogenetic tree.
#start 10:30pm, end 10:38
qiime alignment mask \
--i-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza

#generate tree from masked alignment
#start 10:39, end 10:45
qiime phylogeny fasttree \
--i-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza

#apply midpoint rooting to place the root of the tree at the midpoint of the longest tip-to-tip distance in the unrooted tree
#start 2:04pm, and 2:05pm 
qiime phylogeny midpoint-root \
--i-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

#export
qiime tools export \
rooted-tree.qza \
--output-dir exported-rooted-tree

###### Assign taxonomy #####

#note: I tried using blast just to see what I would get. I accepted a lot of defaults here because I couldn't figure out what they were in the last 97% analysis. These defaults might actually be wrong b/c thre is one about percent identity whose default is 80%. There is a post saying blast is terrible and will classify anything, no matter how bad the fit is: https://github.com/benjjneb/dada2/issues/323
#I don't think I'll use blast or uclust, I feel like they work with OTU clustering methods, not with amplicon sequence variant analysis like DADA2.

#start 3:25, end 3:36
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA111_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy6.qza


#explanation of sklearn and training 
#https://forum.qiime2.org/t/classifier-training-questions/1162/3

qiime tools export \
taxonomy6.qza \
--output-dir exported-taxonomy6

#navigate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy6.qza \
--o-visualization taxonomy6.qzv

qiime tools view taxonomy6.qzv


#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy6.qza \
--m-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots6.qzv

qiime tools view taxa-bar-plots6.qzv
#you can download a csv file from this visualization, might be useful

#I did a bunch of trials with assigning taxonomy using different databases and classifiers (code not shown above) here are my notes: there are unassigned taxa in the file that was blasted against the all database, there are no unassigned taxa in the file that blasted agains the Euk only database. this might be because when you blast against euks it just leavs out anything that doesnot have a hit (b/c they are also likely bacteria so it doesn't make sense to call them unassigned ??)
#looking only at eukaroyte classified reads, many of the numbers in the groups (level2) are identical, a few are a little different, with the Euk only database always having more reads (if they are not the same). the euk dataset also has many many more reads for the Eukaryota;__ classification. That is where the main difference is. This is because it looks like any read that is unclassified or classified as bacteria in the dataset (from the all database) is calssified as a Eukaryota;__ in the euk only database. (i.e. the total number of reads is the same). Looking at some of the taxonomic groups at Level 2, some of the numbers of reads are higher in the all data base, some higher in the euk-only database. I can't explain that, except it is a different training so there is variability
#comparing concensus and majority at level 2 they are nearly identical, does not affect the number of reads in Eukaryota;__
#comparing truncated vs not, they look almost identical
#comparing release 128 to 111, 111 has fewer unassigned, about the same number of Eukaryota;__





##### Bact #####
##### Demultiplexing #####

#Import the data into an artifact
#Navigate to: /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figures\&Stats/kingdata/QIIME2/Bact

#first rename your files to (R1) forward.fastq, (R2) reverse.fastq, (I) and barcodes.fastq
#then gzip them, takes a couple minutes
gzip barcodes.fastq
gzip reverse.fastq 
gzip forward.fastq

#start 
qiime tools import \
--type EMPPairedEndSequences \
--input-path emp-paired-end-sequences \
--output-path emp-paired-end-sequences.qza

#move mapping file into directory, and delete .txt and replace with .tsv

#deleted N.47.2015 from the mapping file b/c it was giving errors for trimming
#barcode is on the forward read, it is not the reverse complement
#start 3:04pm, end 4:59
qiime demux emp-paired \
--m-barcodes-file 515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.tsv \
--m-barcodes-category BarcodeSequence \
--i-seqs emp-paired-end-sequences.qza \
--o-per-sample-sequences demux 

#summarize, start 6:21pm, end
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

qiime tools view demux.qzv

#export it so you can look inside at the files
#started 6:22pm, end
qiime tools export \
demux.qza \
--output-dir exported-demux

##### Trim primers/adapters ##### 
#To get these adapters I looked in the raw sequnce data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
#start 5:45pm, end 6:20pm, ran fine, no error
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux.qza \
--p-cores 4 \
--p-adapter-f ATTAGAWACCCBNGTAGTCC \
--p-adapter-r TTACCGCGGCKGCTGRCAC \
--o-trimmed-sequences trimmed-seqs.qza \
--verbose

#summarize
qiime demux summarize \
--i-data trimmed-seqs.qza \
--o-visualization trimmed-seqs.qzv

qiime tools view trimmed-seqs.qzv
#nice, the samples with low # sequences here matches that in the DataCleaning file (samples 5, 34, 126)

#export
qiime tools export \
trimmed-seqs.qza \
--output-dir exported-trimmed-seqs
#R1 almost all got trimmed to 253 bp. R2 some got trimmed and some didn't probably b/c there are lots of errors at the end of the read and it did not recognize the primer.

#Getting a histogram of line lengths
awk '{print length}' S.0.2015_5_L001_R1_001copy.fastq | sort | uniq -c
#read 1 80036 at 253, 572 at 252, 84 at 251, 2 at 250

awk '{print length}' S.99.2015_101_L001_R1_001copy.fastq | sort | uniq -c
#read 1 96184 at 253, 1066 at 252, 224 at 251, 34 at 250

awk '{print length}' S.0.2015_5_L001_R2_001copy.fastq | sort | uniq -c
#read 2 56518 at 253, 300 at 252, 54 at 251, 28 at 242

awk '{print length}' S.99.2015_101_L001_R2_001copy.fastq | sort | uniq -c
#read 2 65670 at 253, 550 at 252, 142 at 251, 24 at 250

##### Denoising with DADA2 #####
#this is where I start using qiime 19.1 and redoing the analysis
#the vast majority of reads are 253bp
#Before I truncated at p-trunc-len-f 251, and p-trunc-len-r 233. This was because if the primers are 20 and 19 bp, then for R1 the primer should be within the high quality region, and I chose 251 to keep some of the shorter fragments. For R2, the primer will be within the poor quality region, and quality really got bad around 233bp
#after playing with DADA2 in R and looking at the moving pictures tutorial, they truncated more and had higher standards for quality, before I was looking at 25 but now I was going to 30-31. So I decided to truncate at 235 for the forward read and 180 for the reverse.
#note for n-threads specifying 0 means use all cores
#note for trunc-len specifying 0 means don't truncate
#start 8:41pm, end 7:25am; start 2:34pm and 11:04
qiime dada2 denoise-paired \
--i-demultiplexed-seqs trimmed-seqs.qza \
--o-table tableN \
--o-representative-sequences rep-seqsN \
--o-denoising-stats denoising-statsN \
--p-n-threads 0 \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 235 \
--p-trunc-len-r 180

#number of seqs out for S.99.2015 from DADA2 R script is 36180
#from this analyis: 37782, but thre were sequences that wre only 33 bp

#export rep seqs just to take a look at. especially to look at how the long sequences got identified, I don't think I can look at this, since there isn't a file that has the ID DNA header, the sequence, and the feature ID in it

qiime tools export \
--input-path rep-seqsN.qza \
--output-path exported-rep-seqsN

#to print the line with the longest number of characters:
cat dna-sequences.fasta|awk '{print length, $0}'|sort -nr|head -1
#the longest line in the rep set is 460bp, and blasts to a fungus, odd

#26277 are 253 bp
awk '{print length}' dna-sequences.fasta | sort | uniq -c

#export table to get otu table
qiime tools export \
--input-path tableN.qza \
--output-path exported-tableN

#convert biom to otu table text file!
biom convert -i feature-table.biom -o otu_table.txt --to-tsv

#delete the space and # from '#OTU ID' on line 2 
sed '2s/[ ]//' otu_table.txt | sed '2s/.//' > otu_table2.txt

qiime feature-table summarize \
--i-table tableN.qza \
--o-visualization tableN.qzv \
--m-sample-metadata-file 515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.tsv

qiime tools view tableN.qzv

qiime feature-table tabulate-seqs \
--i-data rep-seqsN.qza \
--o-visualization rep-seqsN.qzv

qiime tools view rep-seqsN.qzv


##### Create a phylogenetic tree #####

#do the alignment, mask (filter) positions that are highly variable and will add noise to a phylogenetic tree, make tree from masked alignment, make an unrooted and rooted tree
#start 9:52am, end 11:58am

qiime phylogeny align-to-tree-mafft-fasttree \
--p-n-threads 0 \
--i-sequences rep-seqsN.qza \
--o-alignment aligned-rep-seqsN.qza \
--o-masked-alignment masked-aligned-rep-seqsN.qza \
--o-tree unrooted-treeN.qza \
--o-rooted-tree rooted-treeN.qza

#export
qiime tools export \
--input-path rooted-treeN.qza \
--output-path exported-rooted-treeN


##### Assign taxonomy #####

#using greengenes
#start 2:45pm, end 3:03pm

qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_gg_13_8_classifier.qza \
--i-reads rep-seqsN.qza \
--p-n-jobs 2 \
--o-classification taxonomy_ggN.qza

qiime tools export \
--input-path taxonomy_ggN.qza \
--output-path exported-taxonomy_ggN



#navigate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy_ggN.qza \
--o-visualization taxonomy_ggN.qzv

qiime tools view taxonomy_ggN.qzv

#visualize barplots
qiime taxa barplot \
--i-table tableN.qza \
--i-taxonomy taxonomy_ggN.qza \
--m-metadata-file 515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.tsv \
--o-visualization taxa-bar-plots_ggN.qzv

qiime tools view taxa-bar-plots_ggN.qzv
#nice! there are very few unassigned and very few assigned to bacteria;__











##### ITS #####

##### Demultiplexing #####
# (from MovingUphill3)

#Import the data into an artifact
#Navigate to: /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS

#first rename your files to (R1) forward.fastq, (R2) reverse.fastq, (I) and barcodes.fastq
#then gzip them, takes a couple minutes
gzip barcodes.fastq
gzip reverse.fastq 
gzip forward.fastq

#import paired end reads
qiime tools import \
--type EMPPairedEndSequences \
--input-path emp-paired-end-sequences \
--output-path emp-paired-end-sequences.qza

#move mapping file into directory, and delete .txt and replace with .tsv
grep "GCATCGATGAAGAACGCAGC" Undetermined_S0_L001_R1_001.fastq
grep "GTGTAGATCTCGGTGGTCGCCGTATCATT" Undetermined_S0_L001_R2_001.fastq
grep "TTACTTCCTCTAAATGACCAAG" Undetermined_S0_L001_R2_001.fastq
grep "CGCAAATTCGAC" Undetermined_S0_L001_R1_001.fastq

#it is the reverse complement of the barcode in the index file
reverseComplement(DNAString("GTCGAATTTGCG"))
grep "CGCAAATTCGAC" Undetermined_S0_L001_I1_001.fastq

#barcode is on the forward read (see word doc showing how the R1 has the reverse primers in the read, while the R2 has the forward primers. The initial primers got removed, but the sequence is so short that it kept reading the other half of the primer. However the mapping file is correct: "reverse primer" matches to the R2). looking at the barcodes and the barcode file, it looks like it is the reverse complement
#start 9:06pm, end 12:00am
qiime demux emp-paired \
--m-barcodes-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv \
--m-barcodes-category BarcodeSequence \
--i-seqs emp-paired-end-sequences.qza \
--o-per-sample-sequences demux \
--p-rev-comp-mapping-barcodes

#summarize, start 7:53am, end 8:10am
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

qiime tools view demux.qzv

#export it so you can look inside at the files
#started 8:54, end 9:08
qiime tools export \
demux.qza \
--output-dir exported-demux

#The amplicon size should be 230 bp (according to the earth microbiome website). So there will be primers in many reads.


path <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/exported-demuxdada2" 
list.files(path)
fnFs <- sort(list.files(path, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))

#FWD <- "GCATCGATGAAGAACGCAGC" #what is in the forward read, the revesrse complement of hte reverse read
#REV <- "TTACTTCCTCTAAATGACCAAG" # ditto

#these are the actual primers
FWD <- "CTTGGTCATTTAGAGGAAGTAA"
REV <- "GCTGCGTTCTTCATCGATGC"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer) # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse
               (dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString)) # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

#Pre-filter the ambiguous bases (Ns) and trim the right side of the reads, I'm trimming here b/c trimming is based on the read length not 300, so after removing primers, I don't want to trim an additional 80bp.
#start 8:52am, end 8:58am
#filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE,trimRight = c(10,80))

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[296]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[296]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[296]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[296]]))

#yes, the reverse complement of the forward primer is on the reverse read, and vice versa. 

cutadapt <- "/Users/farrer/miniconda3/envs/qiime2-2018.2/bin/cutadapt"
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

#At first I just ran cutadapt looking for the reverse complement of rev on R1, but then after also looking for the rev on R2 it was finding and cutting some of those primers, so I looked for all orientations.
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
#R1.flags <- paste("-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
#R2.flags <- paste("-A", FWD.RC)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

# Run Cutadapt, start 9:12am, end 9:23
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads, -n 1 for one primer
                             #  "-u", 15, "-U", 70,# i tried this to have cutadapt trim but it yielded fewer reads than any of the other trials
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

p=plotQualityProfile(cutFs[200:201])
p + geom_vline(xintercept=290)

p=plotQualityProfile(cutRs[200:201])#not sure why this won't make a fig with 113, no reverse reads
p + geom_vline(xintercept=230)

#Filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#start 3:50pm, end 3:52
out6 <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)##with trimRight on the first use of filterAndTrim above c(10,80) I'm going with this one, the sequnces out are really similar to the original one and I feel better about removing the poor tail for the forward read
out6[200:220,]

#remove the sample where no reads passed filters
filtFs1<-filtFs[-which(filtFs=="/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/exported-demuxdada2/cutadapt/filtered/N.103.2015_296_L001_R1_001.fastq.gz")]
filtRs1<-filtRs[-which(filtRs=="/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/exported-demuxdada2/cutadapt/filtered/N.103.2015_296_L001_R2_001.fastq.gz")]

errF <- learnErrors(filtFs1, multithread = TRUE)
errR <- learnErrors(filtRs1, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)

##there is a little dip at the end of some of the error fits, especialy on errF, I didn't pursue this b/c it isn't too bad, but this is what I could do if I wanted to fix it. https://github.com/benjjneb/dada2/issues/584
#errF2 <- getErrors(errF, detailed=FALSE)
#errF2.qmax <- matrix(errF2[,ncol(errF2)], nrow=nrow(errF2), ncol=ncol(errF2))
#errF3 <- pmax(errF2, errF2.qmax)
#plotErrors(errF3, nominalQ = TRUE)

#Dereplicate identical reads
derepFs <- derepFastq(filtFs1, verbose = TRUE)
derepRs <- derepFastq(filtRs1, verbose = TRUE)

# Name the derep-class objects by the sample names
sample.names1<-sample.names[-which(sample.names=="N.103.2015")]

names(derepFs) <- sample.names1
names(derepRs) <- sample.names1

#Sample inference, start 4:04, end 4:08
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#merge paired reads, start 4:28, end 4:29
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

######construct OTU sequence table#####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(seqtab.nochim)))

dim(seqtab.nochim)
seqtab.nochim[1:2,1:2]

getN <- function(x) sum(getUniques(x))
out6<-out6[-which(rownames(out6)=="N.103.2015_296_L001_R1_001.fastq.gz"),]
track <- cbind(out6, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names1
head(track)


#####Assign taxonomy#####

unite.ref2 <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/sh_general_release_dynamic_s_02.02.2019.fasta" #

#start 10:23, end 11:40ish
taxa2 <- assignTaxonomy(seqtab.nochim, unite.ref2, multithread = TRUE, tryRC = TRUE) #bootstrap default is 50, which should be fine for sequences under 250bp
taxa.print2 <- taxa2 # Removing sequence rownames for display only
rownames(taxa.print2) <- NULL
tail(taxa.print2)
unique(taxa.print2[,1])

#start 9:21, end 11:20
taxa3 <- assignTaxonomy(seqtab.nochim, unite.ref2, multithread = TRUE, minBoot=70, tryRC = TRUE) #bootstrap default is 50, which should be fine for sequences under 250bp
taxa.print3 <- taxa3 # Removing sequence rownames for display only
rownames(taxa.print3) <- NULL
tail(taxa.print3)
unique(taxa.print3[,1])


###### Exporting things #####
write.csv(seqtab.nochim,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/exported-dada2/seqtab.nochim.csv")
#run this then save as a txt file in excel. for some reason the write.table is nt working with a dataframe that has row and column names

write.csv(taxa3,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/exported-dada2/taxa3.csv")





