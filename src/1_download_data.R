## Downloading data

if ( ! dir.exists('data/kim2015/')) {
        dir.create('data/kim2015/')
        url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE69405&format=file&file=GSE69405%5FPROCESSED%5FGENE%5FTPM%5FALL%2Etxt%2Egz'
        download.file(url, destfile = 'data/kim2015/GSE69405_PROCESSED_GENE_TPM_ALL.txt.gz')
}

