
##### ITS #####

##### Demultiplexing #####
# (from MovingUphill3)

#Import the data into an artifact
#Navigate to: /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figures\&Stats/kingdata/QIIME2/ITS

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
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE,trimRight = c(10,80))#trimRight = c(15,70),trimRight = c(0,80)

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

#I need to rerun cutadapt again anyway, b/c for the bacteria I overwrote some of the files b/c I hadn't changed the folder, I deleted all the bacteria files but some had already overwritten the ITS files. 
#I might want to run this again looking for the forward primers as well, on the 16S I tried it with both since my first read had 1 forward primer, and cutadapt is finding a handful of forward primers in each sample

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
                           #  "-u", 15, "-U", 70,# i tried this but it yielded fewer reads than any of the other trials
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
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE) 
#try with a right trim on the reverse read, more reads retained
#using -u and -U terms on cutadapt resulted in fewer reads passing filters compared to any of the trials using no trimRight or with some trimRight. Not sure why. but upshot is don't use -u or -U on cutadapt. (I didn't save this run)
out2 <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), trimRight = c(0,50),
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE) 

out3 <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), trimRight = c(0,80),
                      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE) 

out4 <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  #with trimRight on the first use of filterAndTrim above c(15,70), does better than no trimRight, but not quite as good by ~300 reads compared to trimRight = c(0,80)
out5 <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)##with trimRight on the first use of filterAndTrim above c(0,80)
out6 <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)##with trimRight on the first use of filterAndTrim above c(10,80) I'm going with this one, the sequnces out are really similar to the original one and I feel better about removing the poor tail for the forward read


head(out4)

out[200:220,]
out2[200:220,]
out3[200:220,]
out4[200:220,]
out5[200:220,]
out6[200:220,] ##this is not in the environmtn here b/c my laptop was crashing when doing bioinformatics in dada2 with a large environment. it worked when doing it in its own new environemnt.

#remove the sample where no reads passed filters
filtFs1<-filtFs[-which(filtFs=="/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/exported-demuxdada2/cutadapt/filtered/N.103.2015_296_L001_R1_001.fastq.gz")]
filtRs1<-filtRs[-which(filtRs=="/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/exported-demuxdada2/cutadapt/filtered/N.103.2015_296_L001_R2_001.fastq.gz")]

errF <- learnErrors(filtFs1, multithread = TRUE)
errR <- learnErrors(filtRs1, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)

##there is a little dip at the end of some of the error fits, especialy on errF, I didn't pursue this b/c it isn't too bad, but this is what I could do if I wanted to fix it. https://github.com/benjjneb/dada2/issues/584
errF2 <- getErrors(errF, detailed=FALSE)
errF2.qmax <- matrix(errF2[,ncol(errF2)], nrow=nrow(errF2), ncol=ncol(errF2))
errF3 <- pmax(errF2, errF2.qmax)
plotErrors(errF3, nominalQ = TRUE)


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

temp<-nchar(getSequences(seqtab.nochim))
which(temp>500)
getSequences(seqtab.nochim)[1842]
sum(seqtab.nochim[,2578])
#many of these long sequences are matching things in blast, but others are only matching tiny portions (25 bp beginning and 25 bp end)

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
head(trackold)
trackold<-track

track[111:130,]
trackold[111:130,]

track[,"nonchim"]-trackold[,"nonchim"]


#####Assign taxonomy#####

unite.ref <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/sh_general_release_dynamic_s_01.12.2017.fasta" #
#start 9:52, end 11:00
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
#darn there was an error that it printed out, but I clicked on the file and I couldn't see the error anymore
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
tail(taxa.print)
which(taxa.print[,1]=="k__Plantae")
which(taxa.print[,1]=="k__Rhizaria")
which(taxa.print[,1]=="k__Protozoa")
which(is.na(taxa.print[,1])==T)
rownames(taxa[9619:9620,])
taxa[9619,]
#when I blasted some of these sequences online, they do match "uncultured fungal isolates" not plants or protists, they are long sequences, possibly errors or something but I'm fine with calling them fungi and using the updated unite fungal only database b/c there are more sequences in the updated database

unite.ref2 <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/sh_general_release_dynamic_s_02.02.2019.fasta" #

#start 10:23, end 11:40ish
taxa2 <- assignTaxonomy(seqtab.nochim, unite.ref2, multithread = TRUE, tryRC = TRUE) #bootstrap default is 50, which should be fine for sequences under 250bp
#No ERRORS!! Nice!
taxa.print2 <- taxa2 # Removing sequence rownames for display only
rownames(taxa.print2) <- NULL
tail(taxa.print2)
unique(taxa.print2[,1])

#I tried this twice on my laptop, it ran for over 24 hours each time and never finished, so I assume it might just be too large and taking up too much memory
#unite.ref3 <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/sh_general_release_dynamic_s_all_02.02.2019.fasta" #
#taxa3 <- assignTaxonomy(seqtab.nochim, unite.ref3, multithread = TRUE, tryRC = TRUE) 

#for the final taxonimic assignmnt, I used minBoot=70 to be consistent with qiime2 analysis of bacteria and euks
















##### 16S #####

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

#barcode is on the forward read, it is not the reverse complement
#I deleted the sample N.47.2015 b/c it was giving errors

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

#The amplicon size should be 390 bp (according to the earth microbiome website). So theoretically there shouldn't be primers in the reads, but the first one I checked did have primers in it. So I need to take them out.



#### using DADA2 for further processing #####
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.8")
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

pathb <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/bact/exported-demuxdada2" 
list.files(pathb)
fnFsb <- sort(list.files(pathb, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
fnRsb <- sort(list.files(pathb, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))

#f ATTAGAWACCCBNGTAGTCC #what is in the forward read, the revesrse complement of hte reverse read
#r TTACCGCGGCKGCTGRCAC # ditto

#these are the actual primers
FWDb <- "GTGYCAGCMGCCGCGGTAA"
REVb <- "GGACTACNVGGGTWTCTAAT"

FWD.orientsb <- allOrients(FWDb)
REV.orientsb <- allOrients(REVb)

fnFs.filtNb <- file.path(pathb, "filtN", basename(fnFsb)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtNb <- file.path(pathb, "filtN", basename(fnRsb))

#start 11:28pm, and 11:40pm
filterAndTrim(fnFsb, fnFs.filtNb, fnRsb, fnRs.filtNb, maxN = 0, multithread = TRUE)

rbind(FWD.ForwardReads = sapply(FWD.orientsb, primerHits, fn = fnFs.filtNb[[1]]),
      FWD.ReverseReads = sapply(FWD.orientsb, primerHits, fn = fnRs.filtNb[[1]]),
      REV.ForwardReads = sapply(REV.orientsb, primerHits, fn = fnFs.filtNb[[1]]),
      REV.ReverseReads = sapply(REV.orientsb, primerHits, fn = fnRs.filtNb[[1]]))

#yes, the reverse complement of hte forward primer is on the reverse read, and vice versa. 
#odd, one read has forward and reverse primers 

#cutadapt <- "/Users/farrer/miniconda3/envs/qiime2-2018.2/bin/cutadapt"
#system2(cutadapt, args = "--version") # Run shell commands from R
path.cutb <- file.path(pathb, "cutadapt")
if(!dir.exists(path.cutb)) dir.create(path.cutb)
fnFs.cutb <- file.path(path.cutb, basename(fnFsb))
fnRs.cutb <- file.path(path.cutb, basename(fnRsb))

FWD.RCb <- dada2:::rc(FWDb)
REV.RCb <- dada2:::rc(REVb)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flagsb <- paste("-g", FWDb, "-a", REV.RCb)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flagsb <- paste("-G", REVb, "-A", FWD.RCb)
# Run Cutadapt, start 12:39, end 12:59 (?)
for(i in seq_along(fnFsb)) {
  system2(cutadapt, args = c(R1.flagsb, R2.flagsb, "-n", 2, # -n 2 required to remove FWD and REV from reads, if you put -n 1 it means as soon as cutadapt finds one primer and removes it, it will stop and not look for any more.
                             "-o", fnFs.cutb[i], "-p", fnRs.cutb[i], # output files
                             fnFs.filtNb[i], fnRs.filtNb[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orientsb, primerHits, fn = fnFs.cutb[[1]]),
      FWD.ReverseReads = sapply(FWD.orientsb, primerHits, fn = fnRs.cutb[[1]]),
      REV.ForwardReads = sapply(REV.orientsb, primerHits, fn = fnFs.cutb[[1]]),
      REV.ReverseReads = sapply(REV.orientsb, primerHits, fn = fnRs.cutb[[1]]))

# Forward and reverse fastq filenames have the format:
cutFsb <- sort(list.files(path.cutb, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
cutRsb <- sort(list.files(path.cutb, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format:
#get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.namesb <- unname(sapply(cutFsb, get.sample.name))
head(sample.namesb)

p<-plotQualityProfile(cutFsb[189]) #162 look at S.153.2015_84_L001_R1_001.fastq.gz and unzip see why there are missing values, there are reads that have no bp in them and no quality scores, just blank with a heading
p + geom_vline(xintercept=235)+geom_vline(xintercept=249)
p<-plotQualityProfile(cutRsb[164])
p + geom_vline(xintercept=180) + geom_vline(xintercept=200)+ geom_vline(xintercept=233)

#Filter and trim
filtFsb <- file.path(path.cutb, "filtered", basename(cutFsb))
filtRsb <- file.path(path.cutb, "filtered", basename(cutRsb))

#start 12:17, end 3:30 (for the trimmed one, the untrimmed one took about 1.5 hrs). third one took 9 minutes (??!!)
outb <- filterAndTrim(cutFsb, filtFsb, cutRsb, filtRsb, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE) 
#try with a right trim on the forward and reverse read, more reads retained
out2b <- filterAndTrim(cutFsb, filtFsb, cutRsb, filtRsb, maxN = 0, maxEE = c(2, 2), trimRight = c(50,50),
                      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE) 
out3b <- filterAndTrim(cutFsb, filtFsb, cutRsb, filtRsb, maxN = 0, maxEE = c(2, 2), truncLen = c(249,200),
                       truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
out4b <- filterAndTrim(cutFsb, filtFsb, cutRsb, filtRsb, maxN = 0, maxEE = c(2, 2), truncLen = c(235,180),
                       truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE) #I think i'm fine with using this very stringent one - more reads ~1000 are resulting and it is chopping off more of the bad error regions

head(out4b)

outb[1:20,]
out2[1:20,]

#remove the sample where no reads passed filters
#filtFs1<-filtFs[-which(filtFs=="/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/exported-demuxdada2/cutadapt/filtered/N.103.2015_296_L001_R1_001.fastq.gz")]
#filtRs1<-filtRs[-which(filtRs=="/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/exported-demuxdada2/cutadapt/filtered/N.103.2015_296_L001_R2_001.fastq.gz")]

#start 10:58, end 11:01
errFb <- learnErrors(filtFsb, multithread = TRUE)
errRb <- learnErrors(filtRsb, multithread = TRUE)

plotErrors(errFb, nominalQ = TRUE)
###HEREE

#Dereplicate identical reads start 6:09 end 6:30
derepFsb <- derepFastq(filtFsb, verbose = TRUE)
derepRsb <- derepFastq(filtRsb, verbose = TRUE)

# Name the derep-class objects by the sample names
#sample.names1<-sample.names[-which(sample.names=="N.103.2015")]

names(derepFsb) <- sample.namesb
names(derepRsb) <- sample.namesb

#Sample inference, start 11:10, end 11:22
dadaFsb <- dada(derepFsb, err = errFb, multithread = TRUE)
dadaRsb <- dada(derepRsb, err = errRb, multithread = TRUE)

#merge paired reads, start 4:28, end 4:29
mergersb <- mergePairs(dadaFsb, derepFsb, dadaRsb, derepRsb, verbose=TRUE)

######construct OTU sequence table#####
seqtabb <- makeSequenceTable(mergersb)
dim(seqtabb)

#remove chimeras
seqtab.nochimb <- removeBimeraDenovo(seqtabb, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(seqtab.nochimb)))

dim(seqtab.nochimb)
seqtab.nochimb[1:2,1:2]

#getN <- function(x) sum(getUniques(x))
#out1<-out3[-which(rownames(out)=="N.103.2015_296_L001_R1_001.fastq.gz"),]
trackb4 <- cbind(out4b, sapply(dadaFsb, getN), sapply(dadaRsb, getN), sapply(mergersb, getN), rowSums(seqtab.nochimb))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(trackb4) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(trackb4) <- sample.namesb
head(trackb)
head(trackold)
trackb4
trackb3
trackb2
trackbold
#all outputs that I tried had no pattern between plant density or diveristy and bacterial richness, I chose to use the settings for trackb4 (truncate at 235 and 180 b/c it yielded the greatest number of reads)
