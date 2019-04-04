
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


#### using DADA2 for further processing #####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

path <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/ITS/exported-demuxdada2" 
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

#start 11:10pm, and 11:16pm
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)


primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#yes, the reverse complement of hte forward primer is on the reverse read, and vice versa. 

cutadapt <- "/Users/farrer/miniconda3/envs/qiime2-2018.2/bin/cutadapt"
system2(cutadapt, args = "--version") # Run shell commands from R


path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-A", FWD.RC)
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 1, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

i=1
system2(cutadapt, args = c(R1.flags, "-o", fnFs.cut[i], fnFs.filtN[i])) 

system2(less, args= c(Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/ITS/exported-demuxdada2/filtN/N.0.2015_259_L001_R1_001.fastq.gz))

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

















# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))




