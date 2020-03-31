## Downloading data
path2project <- '~/sc/sars-cov2'

if ( ! dir.exists('data/kim2015/')) {
        dir.create('data/kim2015/')
        url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE69405&format=file&file=GSE69405%5FPROCESSED%5FGENE%5FTPM%5FALL%2Etxt%2Egz'
        download.file(url, destfile = 'data/kim2015/GSE69405_PROCESSED_GENE_TPM_ALL.txt.gz')
}

if ( ! dir.exists('data/stewart2020/')) {
        dir.create('data/stewart2020/')
        url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE138474&format=file'
        download.file(url, destfile = 'data/stewart2020/GSE138474_RAW.tar')
        untar('data/stewart2020/GSE138474_RAW.tar', exdir = 'data/stewart2020/')
        untar('data/stewart2020/GSM4104144_SC16.LB17028.matrix.tar.gz', 
              exdir = 'data/stewart2020/')
}
