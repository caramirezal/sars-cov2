## Downloading data
path2project <- '/media/ag-cherrmann/cramirez/sars-cov2/'

kim_dir <- paste0(path2project, 'data/kim2015/')
if ( ! dir.exists(kim_dir) ) {
        dir.create(kim_dir)
        url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE69405&format=file&file=GSE69405%5FPROCESSED%5FGENE%5FTPM%5FALL%2Etxt%2Egz'
        download.file(url, 
                      destfile = paste0(kim_dir, 'GSE69405_PROCESSED_GENE_TPM_ALL.txt.gz') )
}

stewart_dir <- paste0(path2project, 'data/stewart2020/')
if ( ! dir.exists(stewart_dir) ) {
        dir.create(stewart_dir)
        url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE138474&format=file'
        download.file(url, 
                      destfile = paste0(stewart_dir, 'GSE138474_RAW.tar'))
        untar(paste0(stewart_dir, 'GSE138474_RAW.tar'), 
              exdir = stewart_dir)
}

## Downloading landscape data
#landscape_data_dir <- paste0(path2project, 'data/landscape_data/')
#if ( ! dir.exists(landscape_data_dir) ) {
#        dir.create(landscape_data_dir)
#        url <- 'https://ndownloader.figshare.com/articles/7235471?private_link=2892d1cf91752e08a963'
#        download.file(url, destfile = paste0(landscape_data_dir, 'landscape_data'))
#}


