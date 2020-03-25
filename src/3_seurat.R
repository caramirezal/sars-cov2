## dependencies
library(Seurat)
library(miceadds)   ## for load.Rdata
library(Matrix)


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

path2project <- '~/sc/sars-cov2/'
list.files(path2project)

metadata <- read.table(paste0(path2project, 
                              'data/panglaodb/PanglaoDB/data/metadata.txt'), 
                       sep = ',')
colnames(metadata) <- c('SRA_accession', 'SRS_accession',
                        'Tissue_origin', 'scRNA-seq_protocol',
                        'Species', 'Sequencing_instrument',
                        'Number_of_expressed_genes', 'Median_number_of_expressed_genes',
                        'Number_of_clusters_in_sample', 'is_tumor',
                        'is_primary', 'is_cell_line')

## selecting tissue cell types
tissues <- c('lung', 'intestine')
species <- grepl('Homo sapiens', metadata$Species)
tissues.s <- grepl(paste(tissues, collapse = '|'), tolower(metadata$Tissue_origin))
metadata.s <- metadata[tissues.s & species, ]
dim(metadata.s)
head(metadata.s)

file.name <- paste0('var/www/html/SRA/SRA.final/',
                    metadata.s$SRA_accession[[1]],
                    '_',
                    metadata.s$SRS_accession[[1]], 
                    '.sparse.RData')

tar.path <- '~/Downloads/panglaodb_bulk_x2EfwF.tar'
dir.out <- '~/sc/sars-cov2/temp/'
extract.comm <- paste('tar  -C ', dir.out,
                       '-xvf', tar.path,
                       file.name, sep = ' ')
#system(extract.comm)

load.Rdata(paste0(dir.out, file.name), objname = 'sparseMat')
sample.seu <- CreateSeuratObject(
        counts = sparseMat, 
        project = 'sars-cov2', 
        assay = 'RNA', 
        min.cells = 30, 
        min.features = 300 
)

## extract sparse matrices from panglaodb by providing sra and srs ids
extract_panglaodb <- function(
        sra, 
        srs, 
        tar.path = '~/Downloads/panglaodb_bulk_x2EfwF.tar',
        dir.out = '~/sc/sars-cov2/temp/') {
        if ( file.exists(tar.path) & 
              dir.exists(dir.out) ) {
                file.name <- paste0('var/www/html/SRA/SRA.final/',
                                    sra,
                                    '_',
                                    srs, 
                                    '.sparse.RData')
                
                cat('Executing: ', extract.comm, '\n')
                extract.comm <- paste('tar  -C ', dir.out,
                                      '-xvf', tar.path,
                                      file.name, sep = ' ')
                system(extract.comm)   
        }
}

## define seurat object
define_seurat <- function(file.name, sample_name) {
        load.Rdata(file.name, objname = 'sparseMat')
        sample.seu <- CreateSeuratObject(
                counts = sparseMat, 
                project = sample_name, 
                assay = 'RNA', 
                min.cells = 30, 
                min.features = 300 
        )
        return(sample.seu)
}


#for (i in 1:nrow(metadata.s)) {
#        extract_panglaodb(
#                sra = metadata.s$SRA_accession[i],
#                srs = metadata.s$SRS_accession[i]
#        )
#}

path2sparse <- paste0(path2project, 'temp//var/www/html/SRA/SRA.final/')
files <- list.files(path2sparse)
files.path <- paste0(path2sparse, files)

seurat.list <- list()
for (i in 1:length(files.path)) {
        seurat.list <- c(seurat.list, 
                         define_seurat(files.path[i],
                                     as.character(metadata.s$Tissue_origin[i]))
                         )
}
#seurat.list <- sapply(files.path, function(f) define_seurat(f))
names(seurat.list) <- gsub('~/sc/sars-cov2/temp//var/www/html/SRA/SRA.final/',
                           '', names(seurat.list))
names(seurat.list) <- gsub('\\.sparse\\.RData',
                           '', names(seurat.list))

for (i in 1:length(seurat.list)) {
        seurat.list[[i]] <- NormalizeData(seurat.list[[i]], 
                                          verbose = FALSE)
        seurat.list[[i]] <- FindVariableFeatures(seurat.list[[i]], 
                                                 selection.method = "vst", 
                                                 nfeatures = 5000, verbose = FALSE)
}

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
FeaturePlot(merge.seu, features = 'ACE2-ENSG00000130234.10', reduction = 'umap')
