# Clustering and plotting CyTOF data
This folder contains scripts that are adapted from Nowicka et al. (2017): https://f1000research.com/articles/6-748.

Overall, you can use these scripts to perform FlowSOM-based clustering (with meta-clustering) for a set of FCS files. This can be followed by creation of heatmaps and density plots for marker expression in each cluster, which can help in the labeling of the clusters.
You can then use the script for creating a UMAP of the data and creating relevant visualizations of the clustering in dimesnionality-reduced space. This includes

* UMAPs without coloring - both for individual samples and for the whole set
* UMAPS colored by marker expression
* UMAPs colored by FlowSOM cluster


After analyzing the data, you may subset to look only at certain clusters (perhaps even run FlowSOM again with a different set of markers).

Future: It will also be possible to perform differential analysis with this script to investigate whether abundance/expression changes with different conditions.
