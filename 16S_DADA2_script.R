# DADA2 pipeline https://benjjneb.github.io/dada2/tutorial_1_8.html
# amplicon = 16S

# Load libraries
library(dada2); packageVersion("dada2")

# Primer sequences for reference
FWD <- "GTGYCAGCMGCCGCGGTAA" # 19 bases
REV <- "GGACTACNVGGGTWTCTAAT" # 20 bases

# File parsing
path <- "/home/Data/" # CHANGE ME to the directory containing your demultiplexed fastqs
list.files(path)

# Filter and trim
fnFs <- sort(list.files(path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq"))

# Extract sample names
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# Examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[1:9]) 
plotQualityProfile(fnRs[1:9])

# Quality filtering
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filtering - optimise parameters by run
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19,20), truncLen=c(250,250), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

# Check reads out
write.csv(out, file="/home/output/filtered_out.csv")

# If some reads do not pass the filter - likely just blanks with low read numbers (check these)
# Filter to only samples that passed
exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

# Error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]

# Apply the core sequence-variant inference algorithm to the dereplicated data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # Check the table dimensions
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths

# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab) # Check proportion non-chimera

# Track number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

seqtab.nochim_t <- t(seqtab.nochim)
write.csv(seqtab.nochim_t, file="/home/output/seqtab_16S.csv")

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/home/training datasets/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "/home/output/ASVs_16S_assigntax.fa")

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "/home/output/ASVs_counts_16S_assigntax.tsv", sep="\t", quote=F, col.names=NA)

rownames(taxa) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(taxa, "/home/output/ASVs_taxonomy_16S_assigntax.tsv", sep="\t", quote=F, col.names=NA)

