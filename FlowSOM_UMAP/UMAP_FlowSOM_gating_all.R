#### Script for gating using UMAP plot and/or FlowSOM ####


# Loading libraries and setting a seed to obtain the same result consistently
library(ggplot2)#, lib.loc="/home/people/s134891/R/x86_64-pc-linux-gnu-library/3.5") comment for computerome
library(FlowSOM)
library(umap)
library(flowCore)
library(ConsensusClusterPlus)
library(ggridges)
library(matrixStats)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(lme4)
library(multcomp)
library(openxlsx)

# Read helper functions
source('/home/molimm2/pedersen/CyTOF/Scripts/helper_functions_CYTOF_analysis.R')

#### Variables to set ####
data_dir <- '/t1-data/user/woodcock/CyTOF_data/Christina_processing/'   # Directory containing fcs files to process
out_dir <- '/home/molimm2/pedersen/CyTOF/Data/flowsom_umap/'            # Directory for saving files and plots, make sure it exists

# How many cells from each sample to use for UMAP
n_down = 5000

# Number of metaclusters to use
nmc <- 12

# Reading metadata
meta_data <- read.xlsx('/home/molimm2/pedersen/CyTOF/Data/metadata_batch1-6_samples.xlsx')

# Markers to include (should be formatted as the fcs descrption field but without the channel - e.g. for 89Y_CD45 write CD45)
CD14_CD19_markers <- c('CD45', 'CD14_CD19', 'CD57', 'CD45RA', 'CD4', 'CD127', 'TCRgd', 'Va7.2', 'CD45RO', 'CD56', 'CD161', 'CD3', 'HLADR', 'CD8a', 'CD16')
CD3_markers <- c('CD45', 'CD57', 'CCR6', 'Bcl2', 'CD45RA', 'GranzymeB', 'CD4', 'CD127', 'CXCR5', 'CXCR3', 'FoxP3', 'Ki67', 'CD45RO', 'CCR7', 'CD161', 'CD38', 'CD25', 'HLADR', 'CCR4', 'CD8a')
CD8_markers <- c('CD57', 'Bcl2', 'CD45RA', 'GranzymeB', 'CD39', 'CD127', 'ICOS', 'CXCR5', 'CD95', 'CD103', 'MicrobialTet', 'PD1', 'CXCR3', 'CD27', 'FoxP3', 'CTLA4', 'TumourTet', 'Integrin_b7', 'CD49a', 'Ki67', 'CD45RO', 'CCR7', 'CD161', 'CD38', 'CD49d', 'CD25', 'CCR10', 'HLADR', 'CCR4', 'CLA')
CD4_markers <- c('CD57', 'CCR6', 'Bcl2', 'CD45RA', 'GranzymeB', 'CD39', 'CD127', 'ICOS', 'CXCR5', 'CD95', 'CD103', 'PD1', 'CXCR3', 'CD27', 'FoxP3', 'CTLA4', 'Integrin_b7', 'CD49a', 'Ki67', 'CD45RO', 'CCR7', 'CD161', 'CD38', 'CD49d', 'CD25', 'CCR10', 'HLADR', 'CCR4', 'CLA')
all_markers <- c('CD45', 'CD14_CD19', 'CD57', 'CCR6', 'Bcl2', 'CD45RA', 'GranzymeB', 'CD4', 'CD39', 'CD127', 'ICOS', 'CXCR5', 'CD95', 'CD103', 'TCRgd', 'Va7.2', 'MicrobialTet', 'PD1', 'CXCR3', 'CD27', 'FoxP3', 'CTLA4', 'TumourTet', 'Integrin_b7', 'CD49a', 'Ki67', 'CD45RO', 'CD56', 'CCR7', 'CD161', 'CD38', 'CD3', 'CD49d', 'CD25', 'CCR10', 'HLADR', 'CCR4', 'CLA', 'CD8a', 'CD16')

use_markers <- CD14_CD19_markers

# Specifying files to use
setwd(data_dir)
files <- list.files('.', pattern='\\.fcs', recursive=T); files <- files[grep('batch[1-6]', files)]
files <- files[grep('Unassigned', files, invert = T)]

# Setting names for plots and outputs
data_name <- 'batch1-6'; pretty_name <- 'batches 1-6'


# Define cluster colors (here there are 30 colors)
color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")




#### Data processing ####

# Read all the data files 
fcs_raw <- read.flowSet(files, transformation = FALSE, truncate_max_range = FALSE)

# Get nicer names for markers
panel_fcs <- pData(parameters(fcs_raw[[1]])); panel_fcs$desc <- gsub("\\d+[A-Za-z]+_", "", panel_fcs$desc)

# Transform data using asinh
fcs <- fsApply(fcs_raw, function(x, cofactor=5){
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  expr <- asinh(expr[, use_markers] / cofactor)
  exprs(x) <- expr
  x
})

# Extract expression data
combined_expr <- fsApply(fcs, exprs)


# Normalizing data (between 0 and 1) for better heatmaps
combined_expr01 <- rescale_data(combined_expr, 0, 1)



# Prepare data for FlowSOM
fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)

# Build SOM
# set.seed(123)
# som <- BuildSOM(fsom, colsToUse = use_markers, silent = F)
# save(som, file = file.path(out_dir, paste0('FlowSOM_SOM_', data_name, '.RData')))

# Load pre-made SOM
load(file = file.path(out_dir, paste0('FlowSOM_SOM_', data_name, 'files.RData')))

## Get the cell clustering into 100 SOM codes
cell_clustering_som <- som$map$mapping[,1]
codes <- som$map$codes

# Get directory for plots from meta-clustering
plot_outdir <- file.path(out_dir, paste0("consensus_plots_", data_name)); if (!dir.exists(file.path(out_dir, paste0("consensus_plots_", data_name)))) {dir.create(file.path(out_dir, paste0("consensus_plots_", data_name)))}

## Meta-clustering into meta-clusters with ConsensusClusterPlus
# mc <- ConsensusClusterPlus(t(codes), maxK = 20, reps = 100,
#                            pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png",
#                            clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
#                            distance = "euclidean", seed = 123)
# 
# save(mc, file = file.path(out_dir, paste0('FlowSOM_MC_', data_name, '.RData')))
# Load meta-clusters
load(file = file.path(out_dir, paste0('FlowSOM_MC_', data_name, 'files.RData')))


# ## Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[cell_clustering_som]





# # Heatmap for FlowSOM clusters
# png(paste0(out_dir, '/heatmap_', nmc, '.png'), width = 6, height = 6, res = 300, units = 'in')
# plot_clustering_heatmap_wrapper(expr = combined_expr[, use_markers],
#                                expr01 = combined_expr01[, use_markers],
#                                cell_clustering = cell_clustering1, color_clusters = color_clusters)
# dev.off()


# # Distribution plots for markers in all clusters (slow)
# png(paste0(out_dir, '/cluster_distr_', nmc, '.png'), width = 6, height = 6, res = 300, units = 'in')
# plot_clustering_distr_wrapper(expr = combined_expr[, use_markers],
#                               cell_clustering = cell_clustering1)
# dev.off()



# Get sample names
sample_ids <- rep(gsub('.fcs', '', basename(names(fsom$metaData))), fsApply(fcs_raw, nrow))

# Downsampling for UMAP
# Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids)


#### Counting cells in groups of clusters from each file - manual definition ####
# list_of_clusters <- lapply(inds, function(x) { table(cell_clustering1[x]) })
# 
# for (s in unique(sample_ids)) {
#   all <- sum(list_of_clusters[[s]])
#   cd1419 <- sum(list_of_clusters[[s]][c(1,2,3,6,7,8,9,10,11)])
#   cd3 <- sum(list_of_clusters[[s]][c(1,2,6,7,8,9,10,11)])
#   cd4 <- sum(list_of_clusters[[s]][c(7,11)])
#   cd8 <- sum(list_of_clusters[[s]][c(1,2,6)])
# 
#   print(paste(s, all, cd1419, cd3, cd4, cd8, collapse = '\t'))
# }





#### UMAP ####

## How many cells to downsample per-sample
umap_ncells <- pmin(table(sample_ids), n_down)

## Get subsampled indices
set.seed(123)
umap_inds <- lapply(names(inds), function(i){
  sample(inds[[i]], umap_ncells[i], replace = FALSE)
})

umap_inds <- unlist(umap_inds)
umap_expr <- combined_expr[umap_inds, use_markers]

# UMAP of data
#set.seed(123)
#umap <- umap::umap(umap_expr, n_neighbors = 15, min_dist = 0.2, metric = 'euclidean', method = 'naive')
#save(umap, file = paste0(out_dir, '/UMAP_', n_down, '_', data_name, 'samples.RData'))
load(paste0(out_dir, '/UMAP_', n_down, '_', data_name, 'samples.RData'))


# Ggplots of UMAP
# Make data frame to use
dr <- data.frame(UMAP1 = umap$layout[, 1], UMAP2 = umap$layout[, 2],
                 umap_expr)
dr$sample_id <- sample_ids[umap_inds]
dr$Cluster <- factor(cell_clustering1[umap_inds], levels = 1:nmc)



# Plot UMAP (no colors)
# p <- ggplot(dr, aes(x = UMAP1, y = UMAP2)) +
#   geom_point(size = 0.1, alpha = 0.2) +
#   ggtitle(paste0('UMAP of ', n_down, ' cells from each sample in ', pretty_name)) +
#   xlab('UMAP1') + ylab('UMAP2') +
#   theme_bw() +
#   theme(legend.title = element_text(size = 15),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 14))
#ggsave(plot = p, paste0(out_dir, '/umap_', n_down, '_', data_name, 'samples.png'), width = 8, height = 7)


# Plot the UMAP (no colors) per sample
#p2 <- p + facet_wrap(~ sample_id)
#ggsave(plot = p2,  paste0(out_dir, '/umap_', n_down, '_persample_', data_name, 'samples.png'), width = 24, height = 21)



#### Coloring by marker expr. ####
# Plots colored on marker expression for each marker
if (!dir.exists(paste0(out_dir, '/marker_expr_umap_', n_down))) {dir.create(paste0(out_dir, '/marker_expr_umap_', n_down))}

for (i in 1:length(use_markers)) {
  ggplot(dr, aes(x = UMAP1, y = UMAP2, colour = eval(parse(text=use_markers[i])))) +
    geom_point(size = 0.1, alpha = 0.4) +
    scale_color_gradientn(use_markers[i],
                          colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)) +
    ggtitle(paste0('UMAP for ', pretty_name, ' (', n_down, ' cells per file)')) +
    labs(colour=use_markers[i]) +
    xlab('UMAP1') + ylab('UMAP2') +
    theme_bw() +
    theme(legend.title = element_text(size = 15),
          legend.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.key.height = unit(3, "line"))
  #ggsave(paste0(out_dir, '/marker_expr_umap_', n_down, '/umap_', n_down, '_', use_markers[i], '_', data_name, 'samples.png'))
}



# Plot of UMAP colored by FlowSOM clustering
p <- ggplot(dr, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 0.1, alpha = 0.2) +
  ggtitle(paste0('UMAP of ', n_down, ' cells from each sample in ', pretty_name)) +
  xlab('UMAP1') + ylab('UMAP2') +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))
ggsave(plot = p, paste0(out_dir, '/umap_', n_down, '_', nmc, '_', data_name, 'samples.png'), width = 8, height = 7)

# UMAP plot per sample (colored by clustering)
# p2 <- p + facet_wrap(~ sample_id)
# ggsave(plot = p2,  paste0(out_dir, '/umap_', n_down, '_', nmc, '_persample_', data_name, 'samples.png'), width = 24, height = 21)





#### Downstream processing of samples ####

#### Differential analysis on 12 clusters ####


# Getting cluster counts and percentages
counts_table <- table(cell_clustering1, sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100

counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)


# Plot color bars for fractions
ggdf <- melt(data.frame(cluster = rownames(props), props),
             id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")
ggdf$sample_id <- gsub('^X', '', ggdf$sample_id); ggdf$sample_id <- gsub('\\.', '-',  ggdf$sample_id)
ggdf$cluster <- factor(ggdf$cluster, levels = 1:nmc)

# Add info from metadata
mm <- match(ggdf$sample_id, gsub('.fcs', '', meta_data$SampleID))
ggdf$CB <- factor(meta_data$CB[mm])
ggdf$PatientID <- factor(meta_data$PatientID[mm])


p <- ggplot(ggdf, aes(x = sample_id, y = proportion, fill = cluster)) +
            geom_bar(stat = "identity") +
            facet_wrap(~CB, scales = "free_x") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            scale_fill_manual(values = color_clusters)

#ggsave(plot = p,  paste0(out_dir, '/bar_graph_', nmc, '_patient_', data_name, 'samples.png'), width = 24, height = 21)


# Focus only on 81 samples PD1/Combo
sub_counts_table <- counts_table[,colnames(counts_table) %in% gsub('.fcs', '', meta_data$SampleID[meta_data$Group=="Initital"])]
sub_props_table <- t(t(sub_counts_table) / colSums(sub_counts_table)) * 100

sub_counts <- as.data.frame.matrix(sub_counts_table)
sub_props <- as.data.frame.matrix(sub_props_table)

ggdf <- melt(data.frame(cluster = rownames(sub_props), sub_props),
             id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")
ggdf$cluster <- factor(ggdf$cluster, levels = 1:nmc)


mm <- match(ggdf$sample_id, gsub('.fcs', '', meta_data$SampleID))
ggdf$CB <- factor(meta_data$CB[mm])
ggdf$PatientID <- factor(meta_data$PatientID[mm])




formula_glmer_binomial2 <- y/total ~ condition + (1|patient_id) + (1|sample_id)









#### Cluster merging ####
cluster_labels <- c('CD8+ T-cells', 'CD8+ T-cells', 'Other', 'Other', 'Other', 'CD8+ T-cells', 'CD4+ T-cells', 'CIK cells',
                    'MAIT CD8+ T-cells',  'Other',  'CD4+ T-cells', 'Other')
cluster_merging <- cbind.data.frame(seq(1,nmc), factor(cluster_labels)); colnames(cluster_merging) <- c('org_cl', 'new_cl')
# 
# # Get new clusters for data
# cell_clustering2 <- cluster_merging1$new_cl[match(cell_clustering1, cluster_merging1$org_cl)]
# code_clustering2 <- cluster_merging1$new_cl[match(code_clustering1, cluster_merging1$org_cl)]
# 
# 
# 
# 
# # Differential analysis - clinical benefit vs. no benefit/not known
# model.matrix( ~ CB, data = meta_data)
# contrast_names <- c("CBvsNone")
# K <- matrix(c(0,1), nrow = 1, byrow = TRUE, dimnames = list(contrast_names))
# 
# 
# # Cell population abudance
# 
# FDR_cutoff <- 0.05



#### New clustering with FlowSOM starting from CD8+ only cells ####
cl_to_keep <- cluster_merging$org_cl[cluster_merging$new_cl=='CD8+ T-cells']

# Extract only CD8+ T-cells into a new FlowSet
fcs_CD8 <- fsApply(fcs_raw, function(x, cofactor=5){
                   colnames(x) <- panel_fcs$desc
                   expr <- exprs(x)
                   expr <- asinh(expr[, all_markers] / cofactor)
                   cells <- cell_clustering1[inds[gsub('\\.fcs', '', x@description$GUID)][[1]]] 
                   expr <- expr[cells %in% cl_to_keep,]
                   exprs(x) <- expr
                   x
})

use_markers <- CD8_markers


# Prepare data for FlowSOM
fsom_CD8 <- ReadInput(fcs_CD8, transform = FALSE, scale = FALSE)

# Build SOM
# set.seed(123)
# som <- BuildSOM(fsom_CD8, colsToUse = use_markers, silent = F)
# save(som, file = file.path(out_dir, paste0('FlowSOM_SOM_', data_name, '_CD8.RData')))

# Load pre-made SOM
load(file = file.path(out_dir, paste0('FlowSOM_SOM_', data_name, '_CD8.RData')))

## Get the cell clustering into 100 SOM codes
cell_clustering_som <- som$map$mapping[,1]
codes <- som$map$codes

# Get directory for plots from meta-clustering
plot_outdir <- file.path(out_dir, paste0("consensus_plots_CD8_", data_name)); if (!dir.exists(file.path(out_dir, paste0("consensus_plots_CD8_", data_name)))) {dir.create(file.path(out_dir, paste0("consensus_plots_CD8_", data_name)))}

## Meta-clustering into meta-clusters with ConsensusClusterPlus
# mc <- ConsensusClusterPlus(t(codes), maxK = 20, reps = 100,
#                            pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png",
#                            clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
#                            distance = "euclidean", seed = 123)
# 
# save(mc, file = file.path(out_dir, paste0('FlowSOM_MC_', data_name, '_CD8.RData')))
# Load meta-clusters
load(file = file.path(out_dir, paste0('FlowSOM_MC_', data_name, '_CD8.RData')))


# ## Get cluster ids for each cell
nmc = 10
code_clustering2 <- mc[[nmc]]$consensusClass
cell_clustering2 <- code_clustering2[cell_clustering_som]







#### Investigation of the location of tetramer positive cells ####
fcs_all <- fsApply(fcs_raw, function(x, cofactor=5){
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  expr <- asinh(expr[,all_markers] / cofactor)
  exprs(x) <- expr
  x
})
combined_expr_all <- fsApply(fcs_all, exprs)

# Get clusters of tumor_tet_positive cells
tumor_tet_pos <- CD8_expr[,'TumourTet'] > asinh(50/5)
expr_tumortetpos <- CD8_expr[tumor_tet_pos,]
table(cell_clustering2[tumor_tet_pos])

CD8_expr01 <- rescale_data(CD8_expr, 0, 1)

# # Heatmap for FlowSOM clusters
png(paste0(out_dir, '/heatmap_', nmc, '_CD8.png'), width = 6, height = 6, res = 300, units = 'in')
plot_clustering_heatmap_wrapper(expr = CD8_expr[, use_markers],
                               expr01 = CD8_expr01[, use_markers],
                               cell_clustering = cell_clustering2, color_clusters = color_clusters)
dev.off()


