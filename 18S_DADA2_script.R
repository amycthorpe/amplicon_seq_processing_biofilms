# DADA2 for big data pipeline https://benjjneb.github.io/dada2/bigdata_paired.html
# amplicon = 18S

# Load libraries
library(dada2); packageVersion("dada2")

# Primer sequences for reference
FWD <- "CGCGGTAATTCCAGCTCCA" # 19 bases
REV <- "TTGGYRAATGCTTTCGC" # 17 bases

# File parsing
pathF <- "/home/18S/data/forward" # CHANGE ME to the directory containing your demultiplexed forward fastqs
pathR <- "/home/18S/data/reverse" # CHANGE ME
list.files(pathF)
list.files(pathR)

filtpathF <- file.path(pathF, "filtered")
filtpathR <- file.path(pathR, "filtered")

fastqFs <- sort(list.files(pathF, pattern="_R1_001.fastq", full.names=F))
fastqRs <- sort(list.files(pathR, pattern="_R2_001.fastq", full.names=F))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Check quality
plotQualityProfile(fastqFs[20:30]) # add full.names=T to fastqFs and Rs lines for this to work, remember to change back to F
plotQualityProfile(fastqRs[20:30])

# Filtering - optimize parameters by run
out <- filterAndTrim(
  fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
  rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
  trimLeft=c(19,17), truncLen=c(300,300), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE, 
  compress=TRUE, verbose=TRUE, multithread=TRUE
)

write.csv(out, file="/home/18S/output/filtered_out.csv") # CHANGE ME

# File parsing
filtpathF <- "/home/18S/data/forward/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/home/18S/data/reverse/filtered" # CHANGE ME

filtFs <- list.files(filtpathF, pattern="_R1_001.fastq", full.names=TRUE)
filtRs <- list.files(filtpathR, pattern="_R2_001.fastq", full.names=TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names
names(filtRs) <- sample.names

set.seed(100)

# Infer sequence variants - learning error rates is computationally intensive, as it requires multiple iterations of the core algorithm - 1 million 100nt reads (or 100M total bases) is more than adequate
# Learn error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merge paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]], verbose=TRUE)
  ddF <- dada(derepF, err=errF, multithread=TRUE, verbose=TRUE)
  derepR <- derepFastq(filtRs[[sam]], verbose=TRUE)
  ddR <- dada(derepR, err=errR, multithread=TRUE, verbose=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose=TRUE)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # Check table dimensions
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths 

# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE, minFoldParentOverAbundance=4)
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab) # Check proportion non-chimera

seqtab.nochim_t <- t(seqtab.nochim)
write.csv(seqtab.nochim_t, file="/home/18S/output/seqtab_18S.csv") # CHANGE ME

# Track number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
rownames(track) <- sample.names
head(track)

write.csv(track, "/home/18S/output/track_reads_18S.csv") # CHANGE ME

# Merge with optimization samples and assign taxonomy
# Load seq tab
seqtab.optim <- read.table("/home/18S/data/optimisation/seqtab_18S.csv", header=T, row.names=1, sep=",", check.names=FALSE)

# Transpose so samples are rows
seqtab.optim <- t(seqtab.optim)

# Subset to only relevant samples - Zymo extraction kit
seqtab.optim <- seqtab.optim[grep("^Z", rownames(seqtab.optim)), ]

# Merge runs
seqtab.all <- mergeSequenceTables(seqtab.nochim, seqtab.optim)

# Remove chimeras after merge
seqtab.all.nochim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.all.nochim) # Check table dimensions

# Transpose and save merged seqtab
seqtab.all.nochim_t <- t(seqtab.all.nochim)
write.csv(seqtab.all.nochim_t, file="/home/18S/output/seqtab18S_merged.csv")

# Assign taxonomy to merged tables
taxa <- assignTaxonomy(seqtab.all.nochim, "/home/training datasets/pr2_version_5.0.0_SSU_dada2.fasta.gz", multithread=TRUE) # CHANGE ME to path of training dataset

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

asv_seqs <- colnames(seqtab.all.nochim)
asv_headers <- vector(dim(seqtab.all.nochim)[2], mode="character")

for (i in 1:dim(seqtab.all.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "/home/18S/output/ASVs_18S_assigntax_merged.fa")

asv_tab <- t(seqtab.all.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "/home/18S/output/ASVs_counts_18S_assigntax_merged.tsv", sep="\t", quote=F, col.names=NA)

rownames(taxa) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(taxa, "/home/18S/output/ASVs_taxonomy_18S_assigntax_merged.tsv", sep="\t", quote=F, col.names=NA)
