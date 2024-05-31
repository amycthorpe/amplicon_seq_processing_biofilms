# Analysis of biofilms for EA report
# Amplicon - 16S


#### Set up ####
# Define libraries to install
packages <- c("dplyr", "vegan", "ggplot2", "phyloseq", "RColorBrewer", "reshape2", "textshape")

# Install and load libraries
for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# Set path
path <- "/home/16S/output/filtering_rarefaction/EA/EA_rarefied/EA_rarefied_filtered/" # CHANGE AMPLICON FOLDER

set.seed(123)


#### Import data ####
# Data is already rarefied and filtered
counts_tab <- data.frame(data.table::fread(file.path(path, "EA_ASV_counts_16S_rare10000_filt.csv"), header = TRUE, check.names = FALSE, sep = ","),
                         row.names = 1) # CHANGE ACCORDING TO AMPLICON
tax_tab <- read.table(file.path(path, "EA_ASV_tax_16S_rare10000_filt.csv"), header = TRUE, row.names = 1, check.names = FALSE, sep = ",")
seq_meta <- read.table(file.path(path, "EA_meta_16S_rare10000_filt.csv"), header = TRUE, check.names = FALSE, sep = ",")

# Environmental metadata
env_meta <- read.table("/home/EA_biofilm_nextseqruns/metadata/all_wims_summary.csv", header = TRUE, check.names = FALSE, sep = ",")
land_cover <- read.table("/home/EA_biofilm_nextseqruns/metadata/RSN_CEH_Land_Cover_Map_2020_Site_Data.csv", header = TRUE, check.names = FALSE, sep = ",")
habitat <- read.table("/home/EA_biofilm_nextseqruns/metadata/RSN_UKNEA_Broad_Habitat_2020_Site_Data.csv", header = TRUE, check.names = FALSE, sep = ",")
geology <- read.table("/home/EA_biofilm_nextseqruns/metadata/RSN_BGS_Simplified_Geology_Site_Data.csv", header = TRUE, check.names = FALSE, sep = ",")


#### ASV filtering ####
# Remove ASVs with <3 reads
counts_tab <- as.data.frame(t(counts_tab))
counts_tab_sub <- counts_tab[rowSums(counts_tab) > 3, ]

# Remove ASVs that only appear <3 samples
counts_tab_sub2 <- counts_tab_sub[rowSums(counts_tab_sub > 0) >= 3, ]


#### Dealing with unassigned tax ####
# Function to add last identifiable taxonomy after 'unassigned' 
replace_unassigned <- function(row) {
  last_non_unassigned <- NA
  for (i in length(row):1) {
    if (row[i] != "unassigned") {
      last_non_unassigned <- row[i]
      break
    }
  }
  row <- ifelse(row == "unassigned", paste0("unassigned ", last_non_unassigned), row)
  return(row)
}

# Apply function
for (i in 1:nrow(tax_tab)) {
  tax_tab[i, ] <- replace_unassigned(tax_tab[i, ])
}


#### Metadata prep ####
# Format metadata - focus on 90 day mean
env_meta <- env_meta[, c("RSN_ID", "Sample_date", "Determinand", "90_day_mean")]

# Unmelt
env_meta_unmelt <- dcast(env_meta, RSN_ID + Sample_date ~ Determinand)

# Merge seq metadata (sample names) and env metadata by RSN_ID and sample date
meta_merged <- merge(seq_meta, env_meta_unmelt, by = c("RSN_ID", "Sample_date"), all.x = TRUE)

# Merge geodatabase metadata
meta_merged_geo <- merge(land_cover[, c(3, 4)], habitat[, c(3, 4)], by = "Text_RSN_ID", all = TRUE)
meta_merged_geo <- merge(meta_merged_geo, geology[, c(3, 4)], by = "Text_RSN_ID", all = TRUE)

colnames(meta_merged_geo)[1] <- "RSN_ID"
colnames(meta_merged_geo)[2] <- "land_cover"
colnames(meta_merged_geo)[3] <- "habitat"
colnames(meta_merged_geo)[4] <- "geology"

# Merge all metadata
meta_merged_all <- merge(meta_merged, meta_merged_geo, by = "RSN_ID", all.x = TRUE)

# Restore rownames
meta <- meta_merged_all[, -3] %>% `rownames<-`(., meta_merged_all[, 3])


#### Phyloseq object ####
TAX <- tax_table(as.matrix(tax_tab))
OTU <- otu_table(counts_tab_sub2, taxa_are_rows = TRUE)
samples <- sample_data(meta)
ps <- phyloseq(OTU, TAX, samples)

ps

# Relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Group at selected taxonomic level
tax_level <- "Phylum" # CHANGE TAX LEVEL
ps_glom <- tax_glom(ps_rel, tax_level)

# Melt
melt_rel <- psmelt(ps_glom)


#### Subset to top taxa ####
# Get top 15
top15 <- names(sort(taxa_sums(ps_glom), decreasing = TRUE))[1:15]

# Find names of top15
tax_ps <- tax_table(ps_glom)
top_names <- as.vector(tax_ps[top15, tax_level])

# Subset to top15
data_top15 <- transform_sample_counts(ps_glom, function(OTU) OTU / sum(OTU))
top_taxa <- prune_taxa(top15, data_top15)

# Create others grouping (x__ so it appears at the bottom of legends etc)
tax_df <- as.data.frame(tax_tab)
tax_df$Phylum[!(tax_df$Phylum %in% top_names)] <- "x__Others" # CHANGE TAX LEVEL
tax_m <- as.matrix(tax_df)

# New phyloseq object with top taxa and others
TAX <- tax_table(tax_m)
ASV <- otu_table(counts_tab_sub2, taxa_are_rows = TRUE)
samples <- sample_data(meta)
ps_top <- phyloseq(ASV, TAX, samples)

ps_top

# Relative abundance
ps_top_rel <- transform_sample_counts(ps_top, function(x) x / sum(x))

# Group at selected taxonomic level
ps_top_glom <- tax_glom(ps_top_rel, tax_level)

# Melt
melt_rel_top <- psmelt(ps_top_glom)


#### Community composition ####
# Make colour palette with enough colours
colourCount <- length(unique(melt_rel_top$Phylum)) # CHANGE TAX LEVEL
getPalette <- colorRampPalette(brewer.pal(16, "Spectral"))

# Set order of groups
land_cover_order <- c("Saltmarsh", "Supra-Littoral Sediment", # coastal margins
                      "Bog", "Fen, Marsh & Swamp", "Freshwater", # freshwaters
                      "Broadleaved Woodland", "Coniferous Woodland", # woodlands
                      "Heather", "Heather Grassland", # mountains, moorlands and heaths
                      "Calcareous Grassland", "Neutral Grassland", "Acid Grassland", # semi-natural grasslands
                      "Arable and Horticulture", "Improved Grassland", # enclosed farmland
                      "Urban", "Suburban") # urban

geology_order <- c("Calcareous","Chalk (Calcareous)","Peat","Salt","Siliceous")

habitat_order <- c("Coastal Margins","Freshwaters - Openwaters, Wetlands and Floodplains", "Woodlands", "Mountains, Moorlands and Heaths","Semi-natural Grasslands","Enclosed Farmland","Urban")

# Loop for each geodatabase variable to plot
datasets <- list(land_cover = land_cover, habitat = habitat, geology = geology)

for (name in names(datasets)) {
  # Merge by selected variable
  ps_merge <- merge_samples(ps_top, name)
  
  # Relative abundance
  ps_merge_rel <- transform_sample_counts(ps_merge, function(x) x / sum(x))
  
  # Group at selected taxonomic level
  ps_merge_glom <- tax_glom(ps_merge_rel, tax_level)
  
  # Melt
  melt_rel <- psmelt(ps_merge_glom)
  
  # Define order
  if (name == "land_cover") {
    plot_order <- land_cover_order
  } else if (name == "habitat") {
    plot_order <- habitat_order
  } else if (name == "geology") {
    plot_order <- geology_order
  }
  
  # Plot
  plot <- ggplot(melt_rel, aes(x = factor(Sample, levels = plot_order),
                               y = Abundance, fill = Phylum)) + # CHANGE TAX LEVEL
    geom_bar(colour = "black", stat = "identity", size = 0.1) +
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 10, colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8, colour = "black"),
          title = element_text(size = 8, colour = "black")) +
    labs(y = "Relative abundance") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = getPalette(colourCount))
  
  # Save plot
  if (name == "geology") {
    ggsave(paste0("output/community/phylum_", name, ".pdf"), plot, width = 5, height = 5, units = "in")
    ggsave(paste0("output/community/phylum_", name, ".png"), plot, width = 5, height = 5, units = "in", dpi = 600)
  } else if (name == "land_cover") {
    ggsave(paste0("output/community/phylum_", name, ".pdf"), plot, width = 5.5, height = 6, units = "in")
    ggsave(paste0("output/community/phylum_", name, ".png"), plot, width = 5.5, height = 6, units = "in", dpi = 600)
  } else if (name == "habitat") {
    ggsave(paste0("output/community/phylum_", name, ".pdf"), plot, width = 6, height = 8, units = "in")
    ggsave(paste0("output/community/phylum_", name, ".png"), plot, width = 6, height = 8, units = "in", dpi = 600)
  }
}


#### Alpha diversity ####
# Calculate diversity
diversity <- estimate_richness(ps)

# Fix rownames (restore hyphens)
rownames(diversity) <- gsub("\\.", "-", rownames(diversity))
write.csv(diversity, "output/diversity/diversity_table.csv")

# Combine diversity with metadata and restore row names
diversity_meta <- merge(meta, diversity, by = "row.names")
diversity_meta <- column_to_rownames(diversity_meta, "Row.names")

# Loop to plot Shannon diversity
for (name in names(datasets)) {
  
  # Calculate number of samples for each category
  sample_counts <- as.data.frame((table(diversity_meta[[name]])))
  write.csv(sample_counts, file = paste0("output/diversity/frequency_", name, ".csv"), row.names = FALSE)
  
  # Define order
  if (name == "land_cover") {
    plot_order <- land_cover_order
  } else if (name == "habitat") {
    plot_order <- habitat_order
  } else if (name == "geology") {
    plot_order <- geology_order
  }
  
  diversity_meta[[name]] <- factor(diversity_meta[[name]], levels = plot_order)
  
  # Plot
  plot <- ggplot(diversity_meta, aes_string(x = name, y = "Shannon")) +
    geom_boxplot() +
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.title.y = element_text(size = 10, colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8, colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 10, colour = "black", face = "bold"),
          title = element_text(size = 8, colour = "black")) +
    labs(y = "Shannon diversity") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 7), breaks = seq(0, 7, by = 1)) +
    scale_x_discrete(labels = plot_order)
  
  # Save plot
  if (name == "geology") {
    ggsave(paste0("output/diversity/shannon_", name, ".pdf"), plot, width = 160, height = 90, units = "mm")
    ggsave(paste0("output/diversity/shannon_", name, ".png"), plot, width = 160, height = 90, units = "mm", dpi = 600)
  } else if (name == "land_cover") {
    ggsave(paste0("output/diversity/shannon_", name, ".pdf"), plot, width = 160, height = 90, units = "mm")
    ggsave(paste0("output/diversity/shannon_", name, ".png"), plot, width = 160, height = 90, units = "mm", dpi = 600)
  } else if (name == "habitat") {
    ggsave(paste0("output/diversity/shannon_", name, ".pdf"), plot, width = 160, height = 90, units = "mm")
    ggsave(paste0("output/diversity/shannon_", name, ".png"), plot, width = 160, height = 90, units = "mm", dpi = 600)
  }
}

# Loop to plot richness
for (name in names(datasets)) {
  
  # Define order
  if (name == "land_cover") {
    plot_order <- land_cover_order
  } else if (name == "habitat") {
    plot_order <- habitat_order
  } else if (name == "geology") {
    plot_order <- geology_order
  }
  
  diversity_meta[[name]] <- factor(diversity_meta[[name]], levels = plot_order)
  
  # Plot
  plot <- ggplot(diversity_meta, aes_string(x = name, y = "Observed")) +
    geom_boxplot() +
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.title.y = element_text(size = 10, colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8, colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 10, colour = "black", face = "bold"),
          title = element_text(size = 8, colour = "black")) +
    labs(y = "Richness") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1400), breaks = seq(0, 1400, 200)) +
    scale_x_discrete(labels = plot_order)
  
  # Save plot
  if (name == "geology") {
    ggsave(paste0("output/diversity/richness_", name, ".pdf"), plot, width = 160, height = 90, units = "mm")
    ggsave(paste0("output/diversity/richness_", name, ".png"), plot, width = 160, height = 90, units = "mm", dpi = 600)
  } else if (name == "land_cover") {
    ggsave(paste0("output/diversity/richness_", name, ".pdf"), plot, width = 160, height = 90, units = "mm")
    ggsave(paste0("output/diversity/richness_", name, ".png"), plot, width = 160, height = 90, units = "mm", dpi = 600)
  } else if (name == "habitat") {
    ggsave(paste0("output/diversity/richness_", name, ".pdf"), width = 160, height = 90, units = "mm")
    ggsave(paste0("output/diversity/richness_", name, ".png"), width = 160, height = 90, units = "mm", dpi = 600)
  }
}


#### Diversity stats ####
library(dunn.test)
library(multcompView)

# Shannon
# Kruskal-Wallis test
KW_test <- kruskal.test(Shannon ~ land_cover, data = diversity_meta)
KW_test

dunn_test <- dunn.test(diversity_meta$Shannon, g = diversity_meta$land_cover, method = "bonferroni")
dunn_test

pairs <- data.frame(dunn_test[["comparisons"]])
zvalues <- data.frame(dunn_test[["Z"]])
pvalues <- data.frame(dunn_test[["P"]])
padjust <- data.frame(dunn_test[["P.adjusted"]])
dunn_results <- cbind(pairs, zvalues, pvalues, padjust)
write.csv(dunn_results, "output/diversity/shannon_dunn.csv")

# Compact letter display
# Dash and spaces interfere with function
dunn_test$comparisons <- gsub("Supra-Littoral Sediment", "Supra Littoral Sediment", dunn_test$comparisons)
Names <- gsub(" ", "", dunn_test$comparisons)

Diff <- dunn_test$P.adjusted < 0.05
names(Diff) <- Names
CLD <- multcompLetters(Diff)
shannon_cld <- data.frame(CLD$Letters)
colnames(shannon_cld)[which(names(shannon_cld) == "CLD.Letters")] <- "Shannon"

# Richness
# Kruskal-Wallis test
KW_test <- kruskal.test(Observed ~ land_cover, data = diversity_meta)
KW_test

dunn_test <- dunn.test(diversity_meta$Observed, g = diversity_meta$land_cover, method = "bonferroni")
dunn_test

pairs <- data.frame(dunn_test[["comparisons"]])
zvalues <- data.frame(dunn_test[["Z"]])
pvalues <- data.frame(dunn_test[["P"]])
padjust <- data.frame(dunn_test[["P.adjusted"]])
dunn_results <- cbind(pairs, zvalues, pvalues, padjust)
write.csv(dunn_results, "output/diversity/richness_dunn.csv")

# Compact letter display
# Dash and spaces interfere with function
dunn_test$comparisons <- gsub("Supra-Littoral Sediment", "Supra Littoral Sediment", dunn_test$comparisons)
Names <- gsub(" ", "", dunn_test$comparisons)

Diff <- dunn_test$P.adjusted < 0.05
names(Diff) <- Names
CLD <- multcompLetters(Diff)
richness_cld <- data.frame(CLD$Letters)
colnames(richness_cld)[which(names(richness_cld) == "CLD.Letters")] <- "Richness"
letters <- cbind(shannon_cld, richness_cld)
write.csv(CLD$Letters, "output/diversity/shannon_richness_cld.csv")


#### NMDS ####
# Extract ASV table and transpose so samples are rows
ASV_tab <- as.data.frame(ps_rel@otu_table)
ASV_t <- t(ASV_tab)

# DCA
decorana(ASV_t) # Check first axis length (<1 indicates low variability and PCA may be best, 1-2 indicates moderate variability and DCA may be best, >2 indicates strong variability and CCA may be best)

# NMDS
ord <- metaMDS(ASV_t, k = 4, maxit = 1000, sfgrmin = 1e-9, sratmax = 1) #,previous.best=ord)
ord

# Save(ord, file = "output/ordination/ord.rda")
# Load("output/ordination/ord.rda")

# Subset metadata
variables <- c("RSN_ID", "Sample_date", "Region", "Area",
               "Ammonia(N)",
               "C - Org Filt",
               "Nitrate-N",
               "Nitrite-N",
               "Orthophospht",
               "Oxygen Diss",
               "pH",
               "SiO2 Rv",
               "Temp Water",
               "Alky pH 4.5",
               "Chloride Ion",
               "land_cover", "habitat", "geology")

meta_sub <- meta[variables]

# Fit environmental variables
ord.fit <- envfit(ord ~ `Ammonia(N)` +
                    `C - Org Filt` +
                    `Nitrate-N` +
                    `Nitrite-N` +
                    `Orthophospht` +
                    `Oxygen Diss` +
                    `pH` +
                    `SiO2 Rv` +
                    `Temp Water` +
                    `Alky pH 4.5` +
                    `Chloride Ion`,
                  data = meta_sub, perm = 999, na.rm = T)
ord.fit

plot(ord, display = "sites")
plot(ord.fit)

plot(ord, display = "sites")
ordihull(ord, meta_sub$land_cover, display = "sites")

# Save dataframe of vector r and p values
df <- ord.fit$vectors
df_r <- as.data.frame(df[["r"]])
df_pvals <- as.data.frame(df[["pvals"]])
df_arrows <- as.data.frame(df[["arrows"]])

ord.fit_df <- cbind(df_arrows, df_r, df_pvals)

write.csv(ord.fit_df, "output/ordination/ord.fit_df.csv")

# Extract scores
nmds_scores <- as.data.frame(ord[["points"]])

# Combine with meta
nmds_scores_meta <- merge(nmds_scores, meta_sub, by = "row.names")
nmds_scores_meta <- column_to_rownames(nmds_scores_meta, "Row.names")

# Set up vector arrows
vectors <- as.data.frame(scores(ord.fit, "vectors")) * ordiArrowMul(ord.fit)
arrowhead = arrow(length = unit(0.02, "npc"))

# Loop to plot NMDS
for (name in names(datasets)) {
  
  # Make colour palette with enough colours
  colourCount = length(unique(nmds_scores_meta[[name]]))
  getPalette = colorRampPalette(brewer.pal(colourCount, "Spectral"))
  
  # Define order
  if (name == "land_cover") {
    plot_order <- land_cover_order
  } else if (name == "habitat") {
    plot_order <- habitat_order
  } else if (name == "geology") {
    plot_order <- geology_order
  }
  
  # Plot
  plot <- ggplot(nmds_scores_meta, aes(x = MDS1, y = MDS2)) +
    theme_bw() +
    geom_point(size = 0.1, alpha = 0.7, stroke = 1, aes_string(x = "MDS1", y = "MDS2", colour = name)) +
    theme(axis.title = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 10, colour = "black"),
          axis.text.y = element_text(size = 10, colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 10, colour = "black"),
          legend.title = element_blank()) +
    geom_vline(xintercept = c(0, 0), linetype = "dashed", colour = "black", size = 0.6) +
    geom_hline(yintercept = c(0, 0), linetype = "dashed", colour = "black", size = 0.6) +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), data = vectors, size = 0.4, arrow = arrowhead) +
    geom_text(data = vectors, aes(x = 1.1 * NMDS1, y = 1.1 * NMDS2), label = row.names(vectors), size = 3, color = "black") +
    labs(x = "NMDS1", y = "NMDS2") +
    scale_colour_manual(values = getPalette(colourCount), breaks = plot_order)
  
  # Save plot
  if (name == "geology") {
    ggsave(paste0("output/ordination/NMDS_", name, ".pdf"), plot, width = 7, height = 5, units = "in")
    ggsave(paste0("output/ordination/NMDS_", name, ".png"), plot, width = 7, height = 5, units = "in", dpi = 600)
  } else if (name == "land_cover") {
    ggsave(paste0("output/ordination/NMDS_", name, ".pdf"), plot, width = 160, height = 110, units = "mm")
    ggsave(paste0("output/ordination/NMDS_", name, ".png"), plot, width = 160, height = 110, units = "mm", dpi = 600)
  } else if (name == "habitat") {
    ggsave(paste0("output/ordination/NMDS_", name, ".pdf"), plot, width = 9, height = 5, units = "in")
    ggsave(paste0("output/ordination/NMDS_", name, ".png"), plot, width = 9, height = 5, units = "in", dpi = 600)
  }
}


#### ANOSIM ####
# Land cover
ASV_dist <- vegdist(ASV_t)
ano <- with(meta, anosim(ASV_dist, grouping = land_cover, permutations = 999, distance = "bray"))
summary(ano)

# Habitat
ano <- with(meta, anosim(ASV_dist, grouping = habitat, permutations = 999, distance = "bray"))
summary(ano)

# Geology
ano <- with(meta, anosim(ASV_dist, grouping = geology, permutations = 999, distance = "bray"))
summary(ano)


#### Bioenv ####
# Subset metadata to numeric variables of interest
variables <- c("Ammonia(N)",
               "C - Org Filt",
               "Nitrate-N",
               "Nitrite-N",
               "Orthophospht",
               "Oxygen Diss",
               "pH",
               "SiO2 Rv",
               "Temp Water",
               "Alky pH 4.5",
               "Chloride Ion")

meta_sub <- meta[variables]

# Remove rows with NA
meta_sub2 <- meta_sub[complete.cases(meta_sub), ]

# Subset to remaining samples
ASV_t2 <- ASV_t[rownames(ASV_t) %in% rownames(meta_sub2), ]

# Run bioenv with Pearson correlations
test <- bioenv(ASV_t2, meta_sub2, method = "pearson", index = "bray", trace = T, metric = "euclidean", parallel = 16)
test

# Run bioenv with Spearman correlations
test2 <- bioenv(ASV_t2, meta_sub2, method = "spearman", index = "bray", trace = T, metric = "euclidean", parallel = 16)
test2


#### Load additional libraries ####
# Some of these packages interfere with packages above and lead to an error in NMDS analysis
library(microeco)
library(file2meco)
library(labdsv)
library(TITAN2)
library(parallel)
library(tidyverse)
library(mikropml)
library(pheatmap)


#### Microeco ####
# Microeco object
meco <- phyloseq2meco(ps)
t1 <- trans_env$new(dataset = meco, add_data = meta_sub)

# Genera correlation heatmap
t1$cal_cor(use_data = "Genus", filter_thres = 0.005, cor_method = "pearson", p_adjust_method = "fdr", p_adjust_type = "Env")
plot <- t1$plot_cor(pheatmap = T)
plot

ggsave("output/community/microeco_genera_heatmap.png", plot, width = 160, height = 190, units = "mm", dpi = 600)
ggsave("output/community/microeco_genera_heatmap.pdf", plot, width = 160, height = 190, units = "mm")

# Save full correlation dataframe
t1$cal_cor(use_data = "Genus", filter_thres = 0.005, cor_method = "pearson", p_adjust_method = "fdr", p_adjust_type = "Env")
cors <- t1[["res_cor"]]
write.csv(cors, "output/community/microeco_genera_cors.csv")

# Diversity correlation heatmap (microeco diversity measures are identical to those calculated with vegan previously)
meco <- phyloseq2meco(ps)
meco$cal_alphadiv()
t2 <- trans_env$new(dataset = meco, add_data = meta_sub)
t2$cal_cor(add_abund_table = meco$alpha_diversity, add_data = meta_sub)
plot <- t2$plot_cor()
plot

ggsave("output/diversity/microeco_diversity_heatmap.png", plot, width = 160, height = 190, units = "mm", dpi = 600)
ggsave("output/diversity/microeco_diversity_heatmap.pdf", plot, width = 160, height = 190, units = "mm")

alpha_cors <- t2$res_cor
write.csv(alpha_cors, "output/diversity/microeco_diversity_cors.csv")

t3 <- trans_alpha$new(dataset = meco, group = "land_cover")
alpha_groups <- data.frame(t3$data_stat)
write.csv(alpha_groups, "output/diversity/microeco_diversity_means.csv")

# Autocorrelation
t1$cal_autocor()


#### Indval ####
# Subset metadata
variables <- c("RSN_ID", "Sample_date", "Region", "Area", "Ammonia(N)",
               "C - Org Filt", "Cond @ 25C", "N Oxidised", "Nitrate-N",
               "Nitrite-N", "O Diss %sat", "Orthophospht", "Oxygen Diss",
               "pH", "Phosphorus-P", "SiO2 Rv", "Temp Water", "Alky pH 4.5",
               "Chloride Ion", "NH3 un-ion", "Nitrogen - N", "land_cover",
               "habitat", "geology")

meta_sub <- meta[variables]

# Run
test <- indval(ASV_t, meta_sub$land_cover)


#### TITAN ####
# Phyloseq object
TAX <- tax_table(as.matrix(tax_tab))
OTU <- otu_table(counts_tab_sub2, taxa_are_rows = T)
samples <- sample_data(meta)
ps <- phyloseq(OTU, TAX, samples)

ps

# Relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Group at selected taxonomic level
tax_level <- "Genus" # CHANGE TAX LEVEL
ps_glom <- tax_glom(ps_rel, tax_level)

# Extract ASV table
ASV_tab <- as.data.frame(ps_glom@otu_table) # CHANGE BETWEEN ps_glom OR ps_rel
ASV_t <- t(ASV_tab)

# Subset metadata to numeric variables of interest
variables <- "Nitrate-N"
meta_sub <- meta[variables]

# Remove rows with NA
meta_sub2 <- as.data.frame(na.omit(meta_sub))

# Subset to remaining samples
ASV_t2 <- as.data.frame(ASV_t[rownames(ASV_t) %in% rownames(meta_sub2),])

# Remove ASVs with occurrence <3
ASV_t3 <- ASV_t2[, colSums(ASV_t2 > 0) >= 3]

# Test a subset
meta_sub2_subset <- head(meta_sub2, 50)
ASV_t3_subset <- ASV_t3[rownames(ASV_t3) %in% rownames(meta_sub2_subset),]
ASV_t3_subset <- ASV_t3_subset[, colSums(ASV_t3_subset > 0) >= 3]

start <- date()
test <- titan(meta_sub2_subset, ASV_t3_subset, minSplt = 5, numPerm = 500, boot = T, nBoot = 500, imax = F, ivTot = F, pur.cut = 0.95, rel.cut = 0.95, ncpus = 16)
end <- date()

# Run on the full sample list grouped at genera level (1173 samples x 1714 genera) - 24-28 hours
start <- date()
titan_genera_nitrate <- titan(meta_sub2, ASV_t3, minSplt = 5, numPerm = 500, boot = T, nBoot = 500, imax = F, ivTot = F, pur.cut = 0.95, rel.cut = 0.95, ncpus = 16)
end <- date()

# saveRDS(titan_genera_nitrate, file = "output/titan/titan_genera_nitrate.rds")
# titan_genera_nitrate <- readRDS("output/titan/titan_genera_nitrate.rds")

# Titan analysis
sppmax_df <- as.data.frame((titan_genera_nitrate[["sppmax"]]))

titan_genera_nitrate$sumz.cp

# Merge with tax tab
tax_glom <- as.data.frame(ps_glom@tax_table) # CHANGE BETWEEN ps_glom OR ps_rel
sppmax_tax <- merge(sppmax_df, tax_tab, by = "row.names")
write.csv(sppmax_tax, "output/titan/titan_genera_nitrate_tax.csv")

# Add taxonomy to titan object
asvs <- colnames(titan_genera_nitrate$taxa)
tax_glom_sub <- t(tax_glom)[, asvs]
titan_genera_nitrate$taxonomy <- tax_glom_sub

genus_row <- as.character(tax_glom_sub["Genus", asvs])
colnames(titan_genera_nitrate$taxa) <- genus_row
rownames(titan_genera_nitrate$sppmax) <- genus_row

# Plot
density <- plot_sumz_density(titan_genera_nitrate, ribbon = FALSE, points = FALSE,
                             axis.text.y = 7,
                             xlabel = "Nitrate-N (mg/L)",
                             xlim = c(-1, 35))
density

taxa <- plot_taxa_ridges(titan_genera_nitrate,
                         axis.text.y = 7,
                         xlabel = "Nitrate-N (mg/L)",
                         xlim = c(-1, 35))
taxa

ggsave("output/titan/titan_genera_nitrate_plot.png", taxa, width = 160, height = 190, units = "mm", dpi = 600)
ggsave("output/titan/titan_genera_nitrate_plot.pdf", taxa, width = 160, height = 190, units = "mm")


# Indicator values
# Only pure and reliable
sppmax_tax2 <- subset(sppmax_tax, filter != 0)

# Genera as rownames
# rownames(sppmax_tax2) <- NULL
# sppmax_tax2 <- tibble::column_to_rownames(sppmax_tax2, "Genus")

# Add minus sign to negative responders
sppmax_tax2$IndVal[sppmax_tax2$filter == 1] <- paste0("-", sppmax_tax2$IndVal[sppmax_tax2$filter == 1])
sppmax_tax2$IndVal <- as.numeric(sppmax_tax2$IndVal)

# Plot
sorted_df <- sppmax_tax2[order(sppmax_tax2$filter, 
                               ifelse(sppmax_tax2$filter == 1, sppmax_tax2$`50%`, -sppmax_tax2$`50%`)), ]

sorted_df$Genus <- factor(sorted_df$Genus, levels = rev(unique(sorted_df$Genus)))

colourCount = length(unique(sorted_df$Phylum))
getPalette = colorRampPalette(brewer.pal(colourCount, "Spectral"))

plot <- ggplot(sorted_df, aes(x = Genus, y = IndVal, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10, colour = "black"),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 8, colour = "black"),
    title = element_text(size = 10, colour = "black")
  ) +
  labs(y = "Indicator value") +
  scale_y_continuous(expand = c(0, 0), limits = c(-70, 70), breaks = seq(-80, 80, by = 20)) +
  coord_flip() +
  scale_fill_manual(values = getPalette(colourCount))
plot

ggsave("output/titan/titan_genera_nitrate_indval_plot.png", plot, width = 160, height = 190, units = "mm", dpi = 600)
ggsave("output/titan/titan_genera_nitrate_indval_plot.pdf", plot, width = 160, height = 190, units = "mm")


#### Mikropml ####
# Nitrates
final_meta <- meta %>% 
  select(`Nitrate-N`) %>% 
  na.omit() %>%
  tibble::rownames_to_column(., var = "sample")

ASV_tab <- as.data.frame(ps_glom@otu_table) # CHANGE BETWEEN ps_glom OR ps_rel
tax_genera <- as.data.frame(ps_glom@tax_table)

final_counts <- ASV_tab %>% # Relative abundance table
  t() %>% as.data.frame() %>%  # Transposing the table
  tibble::rownames_to_column(., var = "sample")

final_counts <- final_counts[rownames(final_counts) %in% rownames(final_meta),] # Only ASVs with metadata

df <- left_join(final_meta, final_counts, by = "sample") %>% 
  tibble::column_to_rownames(., var = "sample")

# Preprocess
preprocessed_df <- preprocess_data(dataset = df, outcome_colname = "Nitrate-N")
preprocessed_df

# Run
results_nitrate_genera <- run_ml(preprocessed_df$dat_transformed,
                                 "glmnet",
                                 outcome_colname = "Nitrate-N",
                                 seed = 42,
                                 find_feature_importance = TRUE
)

performance_nitrate_genera <- results_nitrate_genera$performance
feat_importance_nitrate_genera <- results_nitrate_genera$feature_importance
feat_importance_nitrate_genera_tax <- merge(feat_importance_nitrate_genera, tax_genera, by.x = "feat", by.y = "row.names", all.x = TRUE)

# saveRDS(results_nitrate_genera, file = "output/mikropml/results_nitrate_genera.rds")
# results_nitrate_genera <- readRDS("output/mikropml/results_nitrate_genera.rds")
# write.csv(performance_nitrate_genera, "output/mikropml/performance_nitrate_genera.csv")
# write.csv(feat_importance_nitrate_genera_tax, "output/mikropml/feat_importance_nitrate_genera.csv")

# Plot feature importance
feat_imp_1 <- results_nitrate_genera[["feature_importance"]]
feat_imp_1_tax <- merge(feat_imp_1, tax_genera, by.x = "feat", by.y = "row.names", all.x = TRUE)
perf_metric_name <- results_nitrate_genera[["trained_model"]]$metric
perf_actual <- results_nitrate_genera[["performance"]] %>% pull(perf_metric_name)

feat_imp_1_tax_sub <- feat_imp_1_tax %>%
  filter(perf_metric_diff > 0.001)

colourCount <- length(unique(feat_imp_1_tax_sub$Phylum))
getPalette <- colorRampPalette(brewer.pal(14, "Spectral"))

plot <- feat_imp_1_tax_sub %>%
  mutate(significance = case_when(
    pvalue < 0.001 ~ "***",
    pvalue < 0.01 ~ "**",
    pvalue < 0.05 ~ "*",
    TRUE ~ ""),
    Genus = paste0(significance, Genus)) %>%
  mutate(Genus = as.factor(Genus) %>% forcats::fct_reorder(perf_metric_diff)) %>%
  ggplot(aes(x = perf_metric, xmin = lower, xmax = upper, y = Genus, color = Phylum))+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        title = element_text(size = 10, colour = "black")) +
  geom_pointrange()+
  geom_vline(xintercept = perf_actual, linetype = "dashed")+
  labs(x = "Permutation performance", y = "Genera")+
  scale_color_manual(values = getPalette(colourCount))
plot

ggsave("output/mikropml/feat_importance_plot_nitrate.png", plot, width = 160, height = 190, units = "mm", dpi = 600)
ggsave("output/mikropml/feat_importance_plot_nitrate.pdf", plot, width = 160, height = 190, units = "mm")

