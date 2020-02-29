###################################################################################
## Name: scRNA_processing.R
## Goal: Process single-cell mRNA dataset to select genes to use as predictors for model
##       1) Read in scRNA data as seurat object, filter out cells with likely bad quality
##       2) Normalize data and log tranform
##       3) Select most variable 1000 & 2000 genes
##       4) Subset normalized bulk mRNA dataset to the most variable genes to be used for model building
##       5) Sense check quality of single cell data by MA and UMAP plots
## Author: Claire Su
## Date Created: 02/14/2020
## Date Last Modified: 02/28/2020
## Notes:
###################################################################################
## Load Libraries
library(readr)
library(tidyverse)
library(Seurat)
###################################################################################
## Read in scRNA data using seurat functions
#  Data is in 10X data format
sc_counts <- Read10X(data.dir = "/home/isu/miRNA_project/scRNA_data/SUM149MCF7/filtered_gene_bc_matrices/hg19/")
class(sc_counts)
dim(sc_counts) # 32738  1531
# summary of total expression per single cell
summary(colSums(ctrl_counts))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1096    2746    4246    4894    6406   22815 
hist(colSums(ctrl_counts),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")
###################################################################################
## Create Seurat Object 
#  Filter to genes detected in three or more cells
sc_seurat <- CreateSeuratObject(counts = sc_counts, project = "hg19", min.cells = 3)
class(sc_seurat)
###################################################################################
## Calculate QC metrics and filter out bad samples
# calculate mitochondrail QC measure
sc_seurat[["percent.mt"]] <- PercentageFeatureSet(sc_seurat, pattern = "^MT-")
head(sc_seurat@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(sc_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter out: cells that have unique feature counts over 2,500 or less than 200 
# & cell with > 5% mitochondria count as they are likely low quality
sc_seurat_sub <- subset(sc_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
dim(sc_seurat_sub)#[1] 13843  1265
###################################################################################
## Get intersection of sc and TCGA Genes
# Get TCGA gene names from mrna_selected dataset
mrna_selected<-read.table(paste0("/home/isu/miRNA_project/normalized_data/mrna_selected.csv"), fill = TRUE, header=TRUE,sep=',',stringsAsFactors = FALSE,quote="")
tcga_gene_list <- mrna_selected$gene_id %>% map(function(x) str_split(x,"\\|")[[1]][1]) %>% purrr::flatten() %>% as.data.frame() %>% t() %>% as.data.frame()

# Also create a dataframe to map TCGA gene names to scRNA gene names for later scRNA prediction use
colnames(tcga_gene_list) = "gene"
tcga_gene_id <- data.frame(tcga_gene_id=mrna_selected$gene_id)
tcga_gene_list <- cbind(tcga_gene_list, tcga_gene_id)
tcga_gene_list$gene <- as.character(tcga_gene_list$gene)
tcga_gene_list$tcga_gene_id <- as.character(tcga_gene_list$tcga_gene_id)
saveRDS(tcga_gene_list, file="/home/isu/miRNA_project/intmd_output/tcga_sc_gene_mapping_df.RDS")

# Create a list of genes in both dataset
sc_rna_gene_list <- rownames(x = sc_seurat_sub) %>% as.data.frame()
colnames(sc_rna_gene_list) = "gene"
intersected_gene_list <- intersect(tcga_gene_list$gene, sc_rna_gene_list$gene)
sc_seurat_intersected <- sc_seurat_sub %>% subset(features=intersected_gene_list)
dim(sc_seurat_intersected) # 11228  1265
###################################################################################
## Normalize data by a global-scaling normalization method “LogNormalize” 
sc_seurat_norm_scaled <- NormalizeData(sc_seurat_intersected, normalization.method = "LogNormalize", scale.factor = 1)
###########################################################
## Identify most highly variable gene 
# Directly modeling the mean-variance relationship to control for zero-inflation
# Selecting both 1000 and 2000 most variable genes
sc_selected_2000 <- FindVariableFeatures(sc_seurat_norm_scaled, selection.method = "vst", nfeatures = 2000)
sc_selected_1000 <- FindVariableFeatures(sc_seurat_norm_scaled, selection.method = "vst", nfeatures = 1000)
###########################################################
## Save the list of most variables genes to be used for subset
most_var_2000 <-head(VariableFeatures(sc_selected_2000), 2000)
most_var_1000 <-head(VariableFeatures(sc_selected_1000, 1000)
saveRDS(most_var_2000,file="/home/isu/miRNA_project/intmd_output/most_var_2000.rds")
saveRDS(most_var_1000,file="/home/isu/miRNA_project/intmd_output/most_var_1000.rds")
###########################################################
## Subset mRNA dataset by the most variable genes
# Filter the normalized and matched mRNA dataset to genes that are in the most variable 1000 & 2000 gene sets
# Save the subsetted mRNA datasets
mrna_selected$gene_name <- tcga_gene_list$gene
mrna_selected_new_2000 <- mrna_selected %>% filter(gene_name %in% most_var_2000) #2000 9860
mrna_selected_new_1000 <- mrna_selected %>% filter(gene_name %in% most_var_1000) #1000 9860
row.names(mrna_selected_new_2000) <- mrna_selected_new_2000[,1]
row.names(mrna_selected_new_1000) <- mrna_selected_new_1000[,1]
mrna_selected_new_2000 <- mrna_selected_new_2000[,-c(1,9860)]
mrna_selected_new_1000 <- mrna_selected_new_1000[,-c(1,9860)]
saveRDS(mrna_selected_new_2000, file="/home/isu/miRNA_project/intmd_output/mrna_selected_new_2000.RDS")
saveRDS(mrna_selected_new_1000, file="/home/isu/miRNA_project/intmd_output/mrna_selected_new_1000.RDS")
###########################################################
## Diagnostics/Sense-check Plots
# Scale data to produce dimension reduction plots
sc_selected_2000_scaled <- ScaleData(sc_selected_2000, features = all.genes)
sc_selected_1000_scaled <- ScaleData(sc_selected_1000, features = all.genes)

# Reduce data dimension by PCA
sc_selected_2000_reduced <- RunPCA(sc_selected_2000_scaled, features = VariableFeatures(object = sc_selected_2000_scaled))
sc_selected_1000_reduced <- RunPCA(sc_selected_1000_scaled, features = VariableFeatures(object = sc_selected_1000_scaled))

# Create MA Plots to check make sure that the most variable genes are not zero-inflated
ma_plot_2000 <- VariableFeaturePlot(sc_selected_2000) 
ma_plot_1000 <- VariableFeaturePlot(sc_selected_1000)

# UMAP -- since the expression clustered into two distinct clusters
#         which makes sense since there are two distinct cell lines in the sample
sc_seurat_2000_UMAP <- RunUMAP(sc_selected_2000_reduced, dims = 1:3)
sc_seurat_1000_UMAP <- RunUMAP(sc_selected_1000_reduced, dims = 1:3)
umap_plot_2000 <- DimPlot(sc_seurat_2000_UMAP, reduction = "umap") ## clearly separated into two clusters
umap_plot_1000 <- DimPlot(sc_seurat_1000_UMAP, reduction = "umap") 


