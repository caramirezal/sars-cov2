## Pathway enrichement analysis of DEGs in Calu-3 and Caco-2 cells

## Dependencies
library(dplyr)
library(rio)
library(VennDiagram)
library(enrichR)
library(ggplot2)
library(viridis)
library(VennDetail)
library(UpSetR)

path2project <- '/Users/carlosramirez/sc/sars-cov2/'
setwd(path2project)

set.seed(333)

## Reading DEGs 

## Calu sars-cov2 signature
calu.inf.sign.4hrs <- readRDS('analysis/calu.inf.sign.4hrs.rds') 
head(calu.inf.sign.4hrs)
dim(calu.inf.sign.4hrs)

calu.inf.sign.8hrs <- readRDS('analysis/calu.inf.sign.8hrs.rds')
dim(calu.inf.sign.8hrs)

calu.inf.sign.12hrs <- readRDS('analysis/calu.inf.sign.12hrs.rds')
dim(calu.inf.sign.12hrs)

## H1299 viral signature
h1299.inf.sign.12hrs <- readRDS('analysis/H1299.inf.sign.12hrs.rds')
dim(h1299.inf.sign.4hrs)

h1299.inf.sign.24hrs <- readRDS('analysis/H1299.inf.sign.24hrs.rds')
dim(h1299.inf.sign.24hrs)

h1299.inf.sign.36hrs <- readRDS('analysis/H1299.inf.sign.36hrs.rds')
dim(h1299.inf.sign.36hrs)

## Calu vs h1299 DEGs
calu.vs.h1299.degs <- readRDS('analysis/deg_calu_h1299.rds')
head(calu.vs.h1299.degs)
dim(calu.vs.h1299.degs)

h1299.sign <- readRDS('analysis/deg_wyler/sc_deg_analysis/h1299.infection.signature.rds')
head(h1299.sign)
dim(h1299.sign)

## Number or unique genes in Calu vs H1299 signature
unique(calu.vs.h1299.degs$gene) %>% length()
#[1] 824


#####################################################################
## 
calu.list <- list(calu_vs_h1299 = degs.list$calu.vs.h1299.degs$gene, 
                  calu_inf_sign_12hrs = degs.list$calu.inf.sign.12hrs$gene,
                  h1299_inf_sign = h1299.sign$gene)


filtered_signature <- upset(fromList(calu.list), 
                               sets = c('calu_vs_h1299', 
                                        'calu_inf_sign_12hrs',
                                        'h1299_inf_sign'),
                               keep.order = TRUE #, 
                               #      intersections = list(list('calu_vs_h1299'),
                               #                           list('calu_inf_sign_4hrs'),
                               #                           list('calu_inf_sign_8hrs'),
                               #                           list('calu_inf_sign_12hrs'),
                               #                           list('calu_vs_h1299', 'calu_inf_sign_4hrs'),
                               #                           list('calu_vs_h1299', 'calu_inf_sign_8hrs'),
                               #                           list('calu_vs_h1299', 'calu_inf_sign_12hrs'))
)
filtered_signature
######################################################################
## Upset plot

degs.list <- list(
        calu.inf.sign.4hrs=calu.inf.sign.4hrs,
        calu.inf.sign.8hrs=calu.inf.sign.8hrs,
        calu.inf.sign.12hrs=calu.inf.sign.12hrs,
        h1299.inf.sign.12hrs=h1299.inf.sign.12hrs,
        h1299.inf.sign.24hrs=h1299.inf.sign.24hrs,
        h1299.inf.sign.36hrs=h1299.inf.sign.36hrs,
        calu.vs.h1299.degs=calu.vs.h1299.degs
)

## Function to perform and plot gene set enrichment
genEnrich <- function(geneList){
        up.genes <- geneList %>%
                arrange(desc(avg_logFC)) %>%
                filter(!grepl('SCoV2', gene)) %>%
                head(100) %>%
                select(gene) %>%
                unlist() %>%
                unname()
        enrRF <- enrichr(up.genes, 
                         databases = c("GO_Biological_Process_2018"))
        p <- plotEnrich(enrRF$GO_Biological_Process_2018, 
                   showTerms = 20, 
                   numChar = 50, 
                   y = "Count", 
                   orderBy = "Count") 
        return(p)
}

gea <- lapply(degs.list, genEnrich)
gea <- lapply(1:length(gea), function(i) 
                gea[i][[1]] + 
                      ggtitle(names(gea)[i]) +
                      scale_fill_viridis()
              )
names(gea) <- names(degs.list)

pdf('figures/enrichment_filtered_degs.pdf',
    width = 16, height = 22)
gridExtra::grid.arrange(
        gea$calu.inf.sign.4hrs,
        gea$h1299.inf.sign.12hrs,
        gea$calu.inf.sign.8hrs,
        gea$h1299.inf.sign.24hrs,
        gea$calu.inf.sign.12hrs,
        gea$h1299.inf.sign.36hrs,
        ncol=2
)
dev.off()

###################################################################
## Upset plot

## Intersection of calu vs H1299 and 1299 infection signatures
calu.list <- list(calu_vs_h1299 = degs.list$calu.vs.h1299.degs$gene, 
             #calu_inf_sign_4hrs = degs.list$calu.inf.sign.4hrs$gene,
             #calu_inf_sign_8hrs = degs.list$calu.inf.sign.8hrs$gene,
             calu_inf_sign_12hrs = degs.list$calu.inf.sign.12hrs$gene)


removed_calu_inf_sign <- upset(fromList(calu.list), 
      sets = c('calu_vs_h1299', 
               'calu_inf_sign_4hrs',
               'calu_inf_sign_8hrs', 
               'calu_inf_sign_12hrs'),
      keep.order = TRUE #, 
#      intersections = list(list('calu_vs_h1299'),
#                           list('calu_inf_sign_4hrs'),
#                           list('calu_inf_sign_8hrs'),
#                           list('calu_inf_sign_12hrs'),
#                           list('calu_vs_h1299', 'calu_inf_sign_4hrs'),
#                           list('calu_vs_h1299', 'calu_inf_sign_8hrs'),
#                           list('calu_vs_h1299', 'calu_inf_sign_12hrs'))
)
removed_calu_inf_sign

## Intersection of calu vs H1299 and calu infection signatures
h1299.list <- list(calu_vs_h1299 = degs.list$calu.vs.h1299.degs$gene, 
                  h1299_inf_sign_12hrs = degs.list$h1299.inf.sign.12hrs$gene,
                  h1299_inf_sign_24hrs = degs.list$h1299.inf.sign.24hrs$gene,
                  h1299_inf_sign_36hrs = degs.list$h1299.inf.sign.36hrs$gene)


removed_h1299_inf_sign <- upset(fromList(h1299.list), 
      sets = c('calu_vs_h1299', 'h1299_inf_sign_12hrs',
               'h1299_inf_sign_24hrs', 'h1299_inf_sign_36hrs'),
      keep.order = TRUE #, 
#      intersections = list(list('calu_vs_h1299'),
#                           list('h1299_inf_sign_12hrs'),
#                           list('h1299_inf_sign_24hrs'),
#                           list('h1299_inf_sign_36hrs'),
#                           list('calu_vs_h1299', 'h1299_inf_sign_12hrs'),
#                           list('calu_vs_h1299', 'h1299_inf_sign_24hrs'),
#                           list('calu_vs_h1299', 'h1299_inf_sign_36hrs'),
#                           list('calu_vs_h1299', 'h1299_inf_sign_12hrs', 'h1299_inf_sign_36hrs'),
#                           list('calu_vs_h1299', 'h1299_inf_sign_24hrs', 'h1299_inf_sign_36hrs'),
#                           list('calu_vs_h1299', 'h1299_inf_sign_12hrs', 'h1299_inf_sign_24hrs'))
)

pdf('figures/removed_infection_signatures_upset_plot.pdf')
removed_calu_inf_sign
removed_h1299_inf_sign
dev.off()

###################################################################
## Removed genes from signature

calu_removed_genes <- c(
        intersect(calu.list$calu_vs_h1299, calu.list$calu_inf_sign_4hrs),
        intersect(calu.list$calu_vs_h1299, calu.list$calu_inf_sign_8hrs),
        intersect(calu.list$calu_vs_h1299, calu.list$calu_inf_sign_12hrs)
) %>% unique()
writeLines(calu_removed_genes, 
           'analysis/calu_removed_genes.txt')

h1299_removed_genes <- c(
        intersect(h1299.list$calu_vs_h1299, h1299.list$h1299_inf_sign_12hrs),
        intersect(h1299.list$calu_vs_h1299, h1299.list$h1299_inf_sign_24hrs),
        intersect(h1299.list$calu_vs_h1299, h1299.list$h1299_inf_sign_36hrs)
) %>% unique()
writeLines(h1299_removed_genes, 
           'analysis/h1299_removed_genes.txt')

