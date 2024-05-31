# DADA2 for big data pipeline https://benjjneb.github.io/dada2/bigdata_paired.html
# amplicon = ITS

# Load libraries
library(dada2)
packageVersion("dada2")

# Primer sequences for reference
FWD <- "GTGARTCATCGAATCTTTG" # 19 bases
REV <- "TCCTCCGCTTATTGATATGC" # 20 bases

# File paths
pathF <- "/home/ITS/data/forward" # CHANGE ME to the directory containing your demultiplexed forward fastqs
pathR <- "/home/ITS/data/reverse" # CHANGE ME
list.files(pathF)
list.files(pathR)

# File parsing
filtpathF <- file.path(pathF, "filtered")
filtpathR <- file.path(pathR, "filtered")

fastqFs <- sort(list.files(pathF, pattern="_R1_001.fastq", full.names = FALSE))
fastqRs <- sort(list.files(pathR, pattern="_R2_001.fastq", full.names = FALSE))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Check quality
plotQualityProfile(fastqFs[20:30]) # Add full.names=T to fastqFs and Rs lines for this to work, remember to change back to FALSE
plotQualityProfile(fastqRs[20:30])

# Filtering
out <- filterAndTrim(fwd = file.path(pathF, fastqFs), filt = file.path(filtpathF, fastqFs),
                     rev = file.path(pathR, fastqRs), filt.rev = file.path(filtpathR, fastqRs),
                     trimLeft = c(19,20), truncLen = c(280,270), maxEE = 2, truncQ = 2, maxN = 0, rm.phix = TRUE, # CHANGE ME
                     compress = TRUE, verbose = TRUE, multithread = TRUE)

write.csv(out, file = "/home/ITS/output/filtered_out.csv") # CHANGE ME

# File parsing
filtpathF <- "/home/ITS/data/forward/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/home/ITS/data/reverse/filtered" # CHANGE ME

filtFs <- list.files(filtpathF, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names
names(filtRs) <- sample.names

set.seed(100)

# Infer sequence variants and merge paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]], verbose = TRUE)
  ddF <- dada(derepF, err = errF, multithread = TRUE, verbose = TRUE)
  derepR <- derepFastq(filtRs[[sam]], verbose = TRUE)
  ddR <- dada(derepR, err = errR, multithread = TRUE, verbose = TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose = TRUE)
  mergers[[sam]] <- merger
}
rm(derepF)
rm(derepR)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # Check table dimensions
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths 

# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE, minFoldParentOverAbundance = 4)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # Check proportion non-chimera

seqtab.nochim_t <- t(seqtab.nochim)
write.csv(seqtab.nochim_t, file = "/home/ITS/output/seqtab_ITS.csv") # CHANGE ME

# Track number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
rownames(track) <- sample.names
head(track)

write.csv(track, "/home/ITS/output/track_reads_ITS.csv") # CHANGE ME

# Merge with optimization samples and assign taxonomy
# Load sequence table
seqtab.optim <- read.table("/home/ITS/data/optimisation/seqtab_ITS.csv", header = TRUE, row.names = 1, sep = ",", check.names = FALSE)

# Transpose so samples are rows
seqtab.optim <- t(seqtab.optim)

# Subset to only relevant samples - zymo extraction kit
seqtab.optim <- seqtab.optim[grep("^Z", rownames(seqtab.optim)), ]

# Merge runs
seqtab.all <- mergeSequenceTables(seqtab.nochim, seqtab.optim)

# Remove chimeras after merge
seqtab.all.nochim <- removeBimeraDenovo(seqtab.all, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.all.nochim) # Check table dimensions

# Transpose and save merged seqtab
seqtab.all.nochim_t <- t(seqtab.all.nochim)
write.csv(seqtab.all.nochim_t, file = "/home/ITS/output/seqtabITS_merged.csv")

# Assign taxonomy to merged tables
taxa <- assignTaxonomy(seqtab.all.nochim, "/home/training datasets/sh_general_release_dynamic_all_25.07.2023.fasta", multithread = TRUE) # CHANGE ME to path of training dataset

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

asv_seqs <- colnames(seqtab.all.nochim)
asv_headers <- vector(dim(seqtab.all.nochim)[2], mode = "character")

for (i in 1:dim(seqtab.all.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "/home/ITS/output/ASVs_ITS_assigntax_merged.fa")

asv_tab <- t(seqtab.all.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "/home/ITS/output/ASVs_counts_ITS_assigntax_merged.tsv", sep = "\t", quote = FALSE, col.names = NA)

rownames(taxa) <- gsub(pattern = ">", replacement = "", x = asv_headers)

write.table(taxa, "/home/ITS/output/ASVs_taxonomy_ITS_assigntax_merged.tsv", sep = "\t", quote = FALSE, col.names = NA)
