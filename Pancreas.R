# Required libraies #

library(Seurat)

library(biomaRt)

library(ggplot2)

# Setting working directory #

setwd("C:\\Users\\HCD\\OneDrive\\Masaüstü\\R_Projects\\sc_rna_seq\\Pancreas")

# import RDS file for pancreas

pancreas.combined <- readRDS("pancreas.rds")


pancreas.combined <- FindVariableFeatures(pancreas.combined, selection.method = "vst", nfeatures = 2000)


pancreas.combined <- ScaleData(pancreas.combined, verbose = FALSE)

pancreas.combined <- RunPCA(pancreas.combined, npcs = 30, verbose = FALSE)

pancreas.combined<-FindNeighbors(pancreas.combined, reduction = "pca", dims = 1:23)

pancreas.combined<-FindClusters(pancreas.combined, resolution = 0.1)

pancreas.combined<-RunUMAP(pancreas.combined, reduction = "pca", dims = 1:23)

# Converting ensemble number into gene_name #

{
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")
  ensembl_info <- getBM(attributes = attributes, 
                        filters = "ensembl_gene_id", 
                        values = rownames(pancreas.combined), 
                        mart = ensembl)
  
  current_row_names <- rownames(GetAssayData(object = pancreas.combined))
  
  new_row_names <- ensembl_info$external_gene_name[match(current_row_names, ensembl_info$ensembl_gene_id)]
  rownames(pancreas.combined@assays$RNA@data) <- new_row_names
}


# Naming Idents(Clusters) #

pancreas.combined <- RenameIdents(pancreas.combined, 
                                `0` = "Type A enteroendocrine cells", 
                                `1` = "Pancreatic ductal cells", 
                                `2` = "Acinar cell", 
                                `3` = "Type B pancreatic cell types", 
                                `4` = "Type A enteroendocrine cells", 
                                `5` = "Native cells", 
                                `6` = "Mesenchymal cells", 
                                `7` = "type D enteroendocrine cell") 
                              

DimPlot(pancreas.combined, reduction = "umap",label = TRUE)

# Observation of the gene expression based on cluters #

FeaturePlot(object = pancreas.combined, 
            features = c("desired genes"),
            sort.cell = TRUE,
            min.cutoff = 'q2', 
            label = TRUE,
            repel = TRUE)

# Other observation way of gene expression with DotPlot #

DotPlot(pancreas.combined, features = c("desired genes"),
        cols = c("white","blue"), 
        dot.scale = 5)+
  RotatedAxis()+
  coord_flip()+
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10)) + 
  xlab('Gene') +  
  ylab('Cluster')
