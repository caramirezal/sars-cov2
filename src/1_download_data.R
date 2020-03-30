## Downloading data
library(R.utils)

## Kim2015
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70968
## Lung adenocarcinoma
if ( ! dir.exists('data/kim2015/')) {
        dir.create('data/kim2015/')
        url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE69405&format=file&file=GSE69405%5FPROCESSED%5FGENE%5FTPM%5FALL%2Etxt%2Egz'
        download.file(url, destfile = 'data/kim2015/GSE69405_PROCESSED_GENE_TPM_ALL.txt.gz')
}

## Bulk RNA-Seq of cell lines infected with SARS-COV2
if ( ! dir.exists('data/tenOever2020/') ) {
        dir.create('data/tenOEver202')
        url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE147507&format=file&file=GSE147507%5FRawReadCounts%2Etsv%2Egz'
        download.file(url, 'data/tenOEver202/GSE147507_RawReadCounts.tsv.gz')
        gunzip('data/tenOEver202/GSE147507_RawReadCounts.tsv.gz',
               destname = 'data/tenOEver202/GSE147507_RawReadCounts.tsv',
               remove = FALSE)
}

