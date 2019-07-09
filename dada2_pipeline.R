#install.packages("devtools")
library("devtools")
#devtools::install_github("benjjneb/dada2") # change the ref argument to get other versions
library(dada2)
install.packages("Rcpp")
source("https://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
#biocLite("phyloseq")
library(phyloseq)
library(ShortRead)

##define path and import read names
path <- "fastq"
fns <- list.files(path)
fns
fnFs <- sort(list.files(path, pattern="_L001_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)

##quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
#Assign the filenames for the filtered fastq.gz files.
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, #truncLen=c(240,160),
                     maxN=0, maxEE=c(1,1), truncQ=2, trimLeft=c(20,21),
                     rm.phix=TRUE,
                     compress=TRUE, matchIDs=TRUE,multithread=TRUE)
head(out)
tail(out)

#Learn errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Apply the core sample inference algorithm to the dereplicated data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]

#merge pairs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

#construct sequence table
seqtab <- makeSequenceTable(mergers)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(294,304)] #again, being fairly conservative wrt length
dim(seqtab2)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab2)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab2)

#track reads through pipeline
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "sample_QC_info.csv")
##asign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "fastq/filtered/GeoSymbio_ITS2_LocalDatabase_verForPhyloseq.fasta",  minBoot=10,multithread=TRUE,tryRC=TRUE,outputBootstraps=FALSE)
head(taxa)
dim(taxa)
unname(taxa)
saveRDS(seqtab.nochim, file="seqtab_nochim.rds")
saveRDS(taxa, file="taxa_blastCorrected.rds")
head(taxa)

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing output files
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

asv_tax <- taxa
head(asv_tax)
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

#First produce a blastdatabase with the OTUs
#00080513$ makeblastdb -in ASVs.fa -parse_seqids -dbtype nucl

# Then blast the OTUs against the database
#00080513$ blastn -db ASVs.fa -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 90 -perc_identity 84 -query ASVs.fa

##LULU
library(devtools)
#install_github("tobiasgf/lulu")  
library(lulu)
otutab <- read.csv("ASVs_counts.tsv",sep='\t',header=TRUE,as.is=TRUE, row.names = 1)
matchlist <- read.table("match_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
curated_result <- lulu(otutab, matchlist,minimum_relative_cooccurence=0.01)
curated_result$discarded_otus

write.csv(curated_result$curated_table, "curated_table.csv")
