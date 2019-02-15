#### Functions from Nowicka et al. (2017) for nice plots in CYTOF analysis ####

# Plotting the distribution of marker expression in each cluster 
plot_clustering_distr_wrapper <- function(expr, cell_clustering){
  # Calculate the median expression
  cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  # Calculate cluster frequencies
  freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering,
                            labels = paste0(levels(cell_clustering), "  (", freq_clust, "%)"))
  ### Data organized per cluster
  ggd <- melt(data.frame(cluster = cell_clustering, expr),
              id.vars = "cluster", value.name = "expression",
              variable.name = "antigen")
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  ggd$reference <- "no"
  ### The reference data
  ggd_bg <- ggd
  ggd_bg$cluster <- "reference"
  ggd_bg$reference <- "yes"
  
  ggd_plot <- rbind(ggd, ggd_bg)
  ggd_plot$cluster <- factor(ggd_plot$cluster,
                             levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))
  
  ggplot() +
    geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster,
                                             color = reference, fill = reference), alpha = 0.3) +
    facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +
    theme_ridges() +
    theme(axis.text = element_text(size = 7),  
          strip.text = element_text(size = 7), legend.position = "none")
  
}


# Clustering heatmap function
plot_clustering_heatmap_wrapper <- function(expr, expr01,
                                            cell_clustering, color_clusters, cluster_merging = NULL){
  
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop ,
                       "%)")
  
  # Annotation for the original clusters
  annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  # Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Cluster_merging <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Cluster_merging <- color_clusters2
  }
  
  pheatmap(expr_heat, color = color_heat, cluster_cols = FALSE,
           cluster_rows = cluster_rows, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 8, fontsize_number = 6,  legend_breaks = legend_breaks,
           annotation_row = annotation_row, annotation_colors = annotation_colors)
  
}


# Function for rescaling data to between min and max
rescale_data <- function(expr, min, max) {
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < min] <- min
  expr01[expr01 > max] <- max
  
  return(expr01)
}









