library(tidyverse)
library(SCENIC)

##################################################################
## data preprocessing
maike <- read.csv('data/maike2020/nina_thimme_raw_counts.csv', 
                   header = TRUE)
r_names <- maike$X
maike <- select(maike, -X)
maike_mtx <- apply(maike, 2, as.numeric)
rownames(maike_mtx) <- gsub('__chr.*', '', r_names) 
rm(maike, r_names)

##################################################################
## Dowloading human reference
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")

if ( ! dir.exists('ref')) {
        dir.create('ref/')
}
if ( ! file.exists('ref/hg19-500bp-upstream-7species.mc9nr.feather')){
        for(featherURL in dbFiles)
        {
                download.file(featherURL, 
                              destfile= paste0('ref/', 
                                                   basename(featherURL))
                              ) # saved in current dir
        }  
}



## Setting SCENIC
org="hgnc" # or hgnc, or dmel
dbDir="." # RcisTarget databases location
myDatasetTitle="testingSCENIC" # choose a name for your analysis
dbs <- list.files('ref', full.names = TRUE)
scenicOptions <- initializeScenic(org=org, 
                                  dbDir=dbDir, 
                                  dbs=dbs, 
                                  datasetTitle=myDatasetTitle, 
                                  nCores=4)


### Co-expression network
genes.unique <- unique(rownames(maike_mtx))
maike_mtx.u <- maike_mtx[genes.unique,]
genesKept <- geneFiltering(maike_mtx.u, scenicOptions)
exprMat_filtered <- maike_mtx.u[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions) 
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
