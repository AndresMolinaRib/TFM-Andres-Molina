# Code for different parts of the pipeline

Here is brief description of each script. The order of execution is the one where the scripts are presented. 

## MS-DAP
Code for extracting data from .raw files and written in "dea_final_prueba.xlxs". They are w. It uses the library MS-DAP to obtain legible data from DIA-NN.

## PCA
It separates each of the comparisons, obtains boxplot of normalized data and explore PCA plots.

## Extract stats
Select the metrics of interest from each of the comparisons.

## Venn
Exclude proteins for each of the groups, select common proteins for each of the statistics, calculate the combined metrics and extract the exclusive proteins for each group.

## Enrichment analysis
Do the enrichment analysis and plot the graphics for functions and pathways of down and upregultated proteins in each comparison.

## Heatmap
Generates heatmaps and volcano plots for each condition. It is executed in last term because it is needed to explore the differences and validate the combined statistics method.
