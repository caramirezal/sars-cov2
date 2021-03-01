## Evaluation of TMPRSS2 independent effect of SPINT2 over SARS-CoV-2
## viral load

## Dependencies
library(Seurat)
library(ppcor)
library(pheatmap)
library(ggpubr)

path2project <-  '/Users/carlosramirez/sc/sars-cov2/'
setwd(path2project)

## Auxiliary function to get viral load as the cumulative sum
## of viral reads
get_sign_scores <- function(seurat_object,
                            gene_set, ...){
        gene.subset <- FetchData(seurat_object, 
                                 gene_set, ...)
        signature.scores <- rowSums(gene.subset)
        return(signature.scores)       
}

## Loading SARS-CoV-2 infected Calu-3 cells scRNA-Seq data 
calu <- readRDS('data/200408.Seurat_Calu_CoV_1000_Merged.rds')
calu.scov2.inf <- subset(calu, strain == 'SARSCoV2' & 
                                 infect == 'infected' &
                                 orig.ident %in% c("Calu3-S2-12h-A", 
                                                   "Calu3-S2-12h-B")
)
rm(calu)

## sars-cov genes 
scov2_genes <- grep('CoV2', rownames(calu.scov2.inf), value = T)

## normalization
calu.scov2.inf <- NormalizeData(calu.scov2.inf) %>%
                        ScaleData()

## sars-cov2 expression
calu.scov2.inf$'sarscov2' <- get_sign_scores(calu.scov2.inf,
                                             gene_set = scov2_genes,
                                             slot = 'scale.data')

## Extracting TMPRSS2, SPINT2, ACE2 and viral load expression
gene.exp.df <- calu.scov2.inf %>%
                  FetchData(vars = c('TMPRSS2', 
                           'SPINT2', 'sarscov2')) 
pcors <- pcor(gene.exp.df, method = 'spearman')
pairs(gene.exp.df)
pcors$p.value
pcors$estimate        


## Stratification of TMPRSS2 approach
pdf('figures/spint2_sarscov2_on_tmprss2_stratification.pdf')
gene.exp.df %>%
        ggplot(aes(x=TMPRSS2, y=SPINT2,
                   colour=sarscov2)) +
        geom_point() + 
        theme_light() +
        geom_vline(xintercept = 0.1,
                   linetype = 'dashed')


gene.exp.df <- mutate(gene.exp.df, 
                      tmprss2_up=ifelse(TMPRSS2> 0.1, 
                                        'TMPRSS2 high', 
                                        'TMPRSS2 low'))

## Plotting SARS-CoV-2 vs SPINT2 
gene.exp.df %>% ggscatter(
        iris, x = "sarscov2", y = "SPINT2", 
        palette = "jco",
        add = "reg.line",
        color='steelblue',
        alpha=0.5
) + facet_wrap(~tmprss2_up, 
               scales = 'free') +
        stat_cor(label.y = 4.4) +
        stat_regline_equation(label.y = 4.2)
dev.off()

