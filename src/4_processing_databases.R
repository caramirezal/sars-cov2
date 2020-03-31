## Processing databases
library(Seurat)

## seurat standard preprocessing
seu_std_processing <- function(seurat_object) {
        seurat_object <- NormalizeData(seurat_object) 
        seurat_object <- ScaleData(seurat_object)
        seurat_object <- FindVariableFeatures(seurat_object)
        #seurat_object <- RunPCA(seurat_object, npcs = 30, verbose = F)
        #seurat_object <- RunUMAP(seurat_object, 
        #                   reduction = 'pca', 
        #                   dims = 1:30, 
        #                   metric = 'cosine')
        #seurat_object <- FindNeighbors(seurat_object, 
        #                               dims = 1:30, 
        #                               reduction = 'pca')
        #seurat_object <- FindClusters(seurat_object, 
        #                              resolution = 0.1)
        return(seurat_object)
}

path2project <- '~/sc/sars-cov2/'
## Processing kim 2015 data
kim <- read.table(
        file = paste0(path2project, 
                      'data/kim2015/GSE69405_PROCESSED_GENE_TPM_ALL.txt'),
        colClasses = 'character',
        header = TRUE
)
gene_names <- kim$gene_name
kim.mtx <- apply(kim[, 4:ncol(kim)], 2, as.numeric)
rownames(kim.mtx) <- gene_names 
rm(kim, gene_names)
kim_seu <- CreateSeuratObject(
        counts = kim.mtx, 
        project = 'Lung carcinoma', 
        min.cells = 1, 
        min.features = 300, 
        assay = 'RNA'
)
kim_seu$'cell_type' <- 'Lung carcinoma'
kim_seu <- seu_std_processing(kim_seu)

## Reading Zilionis 2019 data
zilionis <- readMM(paste0(path2project, 
                          '/data/zillionis2019/GSE127465_human_counts_normalized_54773x41861.mtx'))
row_names <- read.table(paste0(path2project, 
                               'data/zillionis2019/human_cell_metadata_54773x25.tsv'), 
                        header = TRUE, 
                        quote = '', 
                        colClasses = 'character', 
                        sep = '\t'
)
row_names <- row_names$Barcode
rownames(zilionis) <- row_names
gene_names <- read_table(paste0(path2project, 
                                '/data/zillionis2019/GSE127465_gene_names_human_41861.tsv'), 
                                col_names = FALSE)
colnames(zilionis) <- gene_names$X1
zilionis <- t(zilionis)
## Subsampling ------- REMOVE on cluster
gene.sample <- sample(1:nrow(zilionis), 500)
samp.sample <- sample(1:ncol(zilionis), 500)
zilionis_seu <- CreateSeuratObject(
        counts = zilionis[gene.sample, samp.sample], 
        project = 'LTI Myeloid cells', 
        assay = 'RNA', 
        min.cells = 1, 
        min.features = 1
)
zilionis_seu$'cell_type' <- 'LTI Myeloid Cells'

