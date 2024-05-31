# Filtering and rarefaction script

# Define libraries to install
packages <- c("plyr", "dplyr", "vegan", "seqinr")

# Install and load libraries
for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# Set working directory
setwd("/home/16S/") # CHANGE ME

set.seed(123)

# Import data
tax_tab <- read.table("output/ASVs_taxonomy_16S_assigntax.tsv", header = TRUE, row.names = 1, check.names = FALSE, sep = "\t") # CHANGE ME
counts_tab <- read.table("output/ASVs_counts_16S_assigntax.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE) # CHANGE ME
meta <- read.table("seq_metadata_16S_2023.csv", header = TRUE, sep = ",", fill = TRUE, check.names = FALSE, row.names = 1) # CHANGE ME

# Correct taxonomy headings (18S only)
# colnames(tax_tab) <- c("Domain", "Supergroup", "Division", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ASVs per sample (check number of reads for any blanks or samples which did not work well)
ASVs_per_sample <- as.data.frame(colSums(counts_tab))
write.csv(ASVs_per_sample, "output/filtering_rarefaction/ASVs_per_sample.csv")

# Rarefaction
# Check row sums (samples need to be rows so transpose first)
counts_tab <- t(counts_tab)
rowsums <- as.data.frame(rowSums(counts_tab))
plot(sort(rowsums[,1]))

min(rowsums[,1])
mean(rowsums[,1])
hist(rowSums(counts_tab))

# Rarefaction curve
plot <- rarecurve(counts_tab, step = 1000, col = "blue", cex = 0.7, subsample = 1000, label = FALSE, xlab = "Number of Sequences",
                  ylab = "Observed Species Richness",
                  main = "Rarefaction Curve")

# Rarefy
counts_tab_rare <- rrarefy(counts_tab, 10000) # CHANGE ME depending on point of richness plateau according to rarefaction curve
counts_tab_rare <- as.data.frame(counts_tab_rare)

# Remove samples below rarefied depth (if any)
counts_tab_rare <- counts_tab_rare[rowSums(counts_tab_rare[]) == 10000,] # CHANGE ME

# Check rarefied rowsums - min and mean should now be the same
rowsums <- as.data.frame(rowSums(counts_tab_rare))
plot(sort(rowsums[,1]))

min(rowsums[,1])
mean(rowsums[,1])

# Remove 0 abundance ASVs from counts after rarefying (if samples below rarefaction depth were removed, these are ASVs only in those samples)
counts_tab_rare <- counts_tab_rare[, colSums(counts_tab_rare[]) > 0]

# Subset tax to new rarefied ASV list (only if you removed samples)
counts_tab_rare <- as.data.frame(t(counts_tab_rare))
tax_tab_rare <- tax_tab[rownames(tax_tab) %in% rownames(counts_tab_rare),]

# Subset meta to rarefied sample list (only if you removed samples)
counts_tab_rare <- as.data.frame(t(counts_tab_rare))
meta_rare <- meta[rownames(meta) %in% rownames(counts_tab_rare),]

# Save rarefied only tables
write.csv(counts_tab_rare, "output/filtering_rarefaction/EA_rarefied/EA_ASV_counts_16S_rare10000.csv") # CHANGE ME
write.csv(tax_tab_rare, "output/filtering_rarefaction/EA_rarefied/EA_ASV_tax_16S_rare10000.csv") # CHANGE ME
write.csv(meta_rare, "output/filtering_rarefaction/EA_rarefied/EA_meta_16S_rare10000.csv") # CHANGE ME

# Filter
# If any remaining blanks, filter out
# meta_rare_filt <- subset(meta_rare, Blank_type != "Blank")
# meta_rare_filt <- subset(meta_rare_filt, Sample_ID != "D6300")
# counts_tab_rare_filt <- counts_tab_rare[rownames(counts_tab_rare) %in% rownames(meta_rare_filt),]
# counts_tab_rare_filt <- counts_tab_rare_filt[, colSums(counts_tab_rare_filt[]) > 0]
# tax_tab_rare_filt <- tax_tab_rare[rownames(tax_tab_rare) %in% colnames(counts_tab_rare_filt),]

# Filter based on taxonomy - 16S
tax_tab_rare[is.na(tax_tab_rare)] <- "unassigned" # NA must be converted to unassigned first
tax_tab_rare_filt <- subset(tax_tab_rare, Kingdom == "Bacteria")
tax_tab_rare_filt <- subset(tax_tab_rare_filt, Order != "Chloroplast")
tax_tab_rare_filt <- subset(tax_tab_rare_filt, Phylum != "unassigned")

# Filter based on taxonomy - ITS
# tax_tab_rare[is.na(tax_tab_rare)] <- "unassigned"
# tax_tab_rare_filt <- subset(tax_tab_rare, Kingdom == "k__Fungi")
# tax_tab_rare_filt <- subset(tax_tab_rare_filt, Phylum != "unassigned")

# Filter based on taxonomy - 18S
# tax_tab_rare[is.na(tax_tab_rare)] <- "unassigned"
# tax_tab_rare_filt <- subset(tax_tab_rare, Domain == "Eukaryota")
# tax_tab_rare_filt <- subset(tax_tab_rare_filt, Phylum != "unassigned")

# Filter based on taxonomy - rbcL
# tax_tab_rare[is.na(tax_tab_rare)] <- "unassigned"
# tax_tab_rare_filt <- subset(tax_tab_rare, Phylum == "Bacillariophyta")

# Subset counts based on filtered tax (so both tax and counts tables have the same samples in after filtering)
counts_tab_rare <- t(counts_tab_rare)
counts_tab_rare_filt <- counts_tab_rare[rownames(counts_tab_rare) %in% rownames(tax_tab_rare_filt),]
counts_tab_rare_filt <- as.data.frame(t(counts_tab_rare_filt))

# Save filtered and rarefied tables
write.csv(counts_tab_rare_filt, "output/filtering_rarefaction/EA_rarefied/EA_rarefied_filtered/EA_ASV_counts_16S_rare10000_filt.csv") # CHANGE ME
write.csv(tax_tab_rare_filt, "output/filtering_rarefaction/EA_rarefied/EA_rarefied_filtered/EA_ASV_tax_16S_rare10000_filt.csv") # CHANGE ME
write.csv(meta_rare, "output/filtering_rarefaction/EA_rarefied/EA_rarefied_filtered/EA_meta_16S_rare10000_filt.csv") # CHANGE ME
