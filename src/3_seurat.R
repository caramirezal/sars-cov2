## Processing single cell data from PanglaoDB repository using Seurat
## Bulk download can be done from here:
## https://panglaodb.se/bulk.html

#######################################################################################
##                                                                                   ##
##                      Initial settings                                             ##
##                                                                                   ##
#######################################################################################

## Dependencies
library(Seurat)
library(miceadds)   ## for load.Rdata
library(Matrix)

## Setting data input

##  path to github repository directory
path2project <- '~/sc/sars-cov2/'     ## change in cluster
## directory to save seurat objects
panglaodb_seurat_dir <- '~/sc/sars-cov2/analyisis/panglaodb_seurat'

#######################################################################################
##                                                                                   ##
##                     Auxiliary functions                                           ##
##                                                                                   ##
#######################################################################################

## Seurat standard workflow
st_workflow <- function(
        seurat_object,
        n_features = 3000,
        n_pca_dims = 15
){
        cat('Normalizing and finding variable features\n')
        seurat.p <- NormalizeData(seurat_object) %>%
                FindVariableFeatures(selection.method = 'vst',
                                     nfeatures = n_features)
        cat('Scaling and projection\n')
        seurat.p <- ScaleData(seurat.p, 
                              verbose = FALSE) %>% 
                RunPCA(npcs = n_pca_dims, 
                       verbose = FALSE) %>%
                RunUMAP(reduction = "pca", 
                        dims = 1:n_pca_dims) 
        return(seurat.p)
}

## define seurat object
define_seurat <- function(file.name, sample_name) {
        load.Rdata(file.name, objname = 'sparseMat')
        colnames(sparseMat) <- paste0(sample_name, '-', 
                                      colnames(sparseMat))
        sample.seu <- CreateSeuratObject(
                counts = sparseMat, 
                project = sample_name, 
                assay = 'RNA', 
                min.cells = 30, 
                min.features = 300 
        )
        return(sample.seu)
}

#######################################################################################
## loading metadata
## Metadata can be obtained by cloning the pangladb repository:
## https://github.com/oscar-franzen/PanglaoDB
## Run 1_download_data.sh in order to get the metadata.txt file in the correct dir
metadata <- read.table(paste0(path2project, 
                              'data/panglaodb/PanglaoDB/data/metadata.txt'), 
                       sep = ',', colClasses = 'character')
colnames(metadata) <- c('SRA_accession', 'SRS_accession',
                        'Tissue_origin', 'scRNA-seq_protocol',
                        'Species', 'Sequencing_instrument',
                        'Number_of_expressed_genes', 'Median_number_of_expressed_genes',
                        'Number_of_clusters_in_sample', 'is_tumor',
                        'is_primary', 'is_cell_line')

#####################################################################################
##                                                                                 ##
##                 Transform GE matrix to seurat object                            ##
##                                                                                 ##
#####################################################################################

## **1. Records with no SRS ID were removed (30)**
## **2. Only count matrices were kept (45 records given in RPKM were left)**
## matrices to Seurat 
for (i in 1:nrow(metadata)){
        file.name <- paste0('var/www/html/SRA/SRA.final/',
                            metadata$SRA_accession[i],
                            '_',
                            metadata$SRS_accession[i], 
                            '.sparse.RData')
        file.name.full <- paste0(path2project,
                                 'temp/',         ## change this in cluster
                                 file.name)
        if (file.exists(file.name.full)) {
                sample.seu <- define_seurat(
                        file.name = file.name.full,
                        sample_name = paste0(metadata$SRA_accession[i],
                                             '_',
                                             metadata$SRS_accession[i])
                )
                sample.seu$'SRA_accession' <- metadata$SRA_accession[i]
                sample.seu$'SRS_accession' <- metadata$SRS_accession[i]
                sample.seu$'Tissue_origin' <- metadata$Tissue_origin[i]
                sample.seu$'scRNA_seq_protocol' <- metadata$`scRNA-seq_protocol`[i]
                sample.seu$'Sequencing_instrument' <- metadata$Sequencing_instrument[i]
                sample.seu$'Number_of_expressed_genes' <- metadata$Number_of_expressed_genes[i]
                
                cat('Finishing ', file.name, 'sparse to seurat conversion\n')
                seurat.file <- paste0(panglaodb_seurat_dir, '/',
                                      metadata$SRA_accession[i], 
                                      '_', metadata$SRS_accession[i], 
                                      '_seurat.rds')
                cat('Saving file to ', seurat.file, '\n')
                write_rds(x = sample.seu, path = seurat.file)
        }
}

seurat.files <- list.files(panglaodb_seurat_dir)
seurat.files.full <- paste0(panglaodb_seurat_dir, '/',seurat.files)
seurat.list <- lapply(seurat.files.full, read_rds)


## selecting tissue cell types
#tissues <- c('lung', 'intestine')
#species <- grepl('Homo sapiens', metadata$Species)
#tissues.s <- grepl(paste(tissues, collapse = '|'), tolower(metadata$Tissue_origin))
#metadata.s <- metadata[tissues.s & species, ]
#dim(metadata.s)
#head(metadata.s)

#####################################################################################
##                                                                                 ##
##              Coembedding samples                                                ##
##                                                                                 ##
####################################################################################
for (i in 1:length(seurat.list)) {
        seurat.list[[i]] <- NormalizeData(seurat.list[[i]], 
                                          verbose = FALSE)
        seurat.list[[i]] <- FindVariableFeatures(seurat.list[[i]], 
                                                 selection.method = "vst", 
                                                 nfeatures = 5000, verbose = FALSE)
}

## taking the intersection of variable features in all samples
var_feat <- lapply(seurat.list, VariableFeatures)
var_fea_int <- Reduce(intersect, var_feat)

merge.seu <- merge(x= seurat.list[[1]], y = seurat.list[1:length(seurat.list)])
merge.seu <- ScaleData(merge.seu, features = var_fea_int, do.scale=TRUE)
merge.seu <- RunPCA(merge.seu, features = var_fea_int, verbose = FALSE)
merge.seu <- RunUMAP(merge.seu, dims = 1:30, metric = 'cosine')
merge.seu <- FindNeighbors(merge.seu, reduction = 'umap', dims = 1:2)
merge.seu <- FindClusters(merge.seu, resolution = 0.05)
DimPlot(merge.seu, group.by = 'seurat_clusters', reduction = 'umap')
DimPlot(merge.seu, group.by = 'orig.ident', reduction = 'umap')
DimPlot(merge.seu, group.by = 'SRA_accession', reduction = 'umap')
DimPlot(merge.seu, group.by = 'Tissue_origin', reduction = 'umap')
FeaturePlot(merge.seu, features = 'ACE2-ENSG00000130234.10', reduction = 'umap')

