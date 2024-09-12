## function
load_env <- function() {
    suppressMessages({
        library(clusterProfiler)
        library(enrichplot)
        library(ReactomePA)
        library(KEGG.db)
        library(org.Hs.eg.db)
        library(org.Mm.eg.db)
        library(dplyr)
    })
}

read_genes <- function(infile, logfoldchanges_threshold, pvals_threshold) {

    difdat <- read.csv(infile, header = TRUE, row.names = 1, sep = '\t')
    if (colnames(difdat)[1] == 'scores') {
        up_loc <- grep("logfoldchanges", colnames(difdat))
        up <- difdat %>%
            subset(logfoldchanges > logfoldchanges_threshold & pvals < pvals_threshold) %>%
            rownames()
    }

    up_geneList <- difdat[up, up_loc]
    names(up_geneList) <- up

    return(
        list(
            up = up,
            upgeneList = up_geneList
        )
    )
}

get_term <- function(obj, topshow, db) {
    if (db == 'GO') {
        obj <- obj %>%
            dplyr::arrange(get(cutoff)) %>%
            dplyr::group_by(ONTOLOGY) %>%
            dplyr::slice(1:topshow)
    }else {
        obj <- obj %>%
            dplyr::arrange(get(cutoff)) %>%
            dplyr::slice(1:topshow)
    }

    uterms <- obj@result$Description

    uterms
}

get_heatplot_df <- function(obj, uterms, geneList) {
    obj <- filter(obj, obj@result$Description %in% uterms)
    heatplot_obj <- heatplot(obj, showCategory = uterms, foldChange = geneList)
    heatplot_df <- heatplot_obj$data
    heatplot_df <- heatplot_df[!duplicated(heatplot_df), ]

    heatplot_df
}

## go
enrich_go <- function(glt) {
    go <- enrichGO(gene = glt, OrgDb = orgdb, keyType = 'ENTREZID', ont = 'ALL', pAdjustMethod = 'BH', readable = T, pvalueCutoff = 1, qvalueCutoff = 1)

    go_topn <- go %>%
        dplyr::filter(get(cutoff) < 0.05, Count >= count) %>%
        dplyr::arrange(get(cutoff)) %>%
        dplyr::group_by(ONTOLOGY) %>%
        dplyr::slice(1:topshow)

    uterms <- get_term(go, topshow, 'GO')
    obj <- filter(go, go@result$Description %in% uterms)
    heatplot_data <- get_heatplot_df(obj, uterms, geneList)

    return(
        list(
            go_topn@result,
            heatplot_data
        )
    )
}


## kegg
enrich_kg <- function(glt) {
   kg <- enrichKEGG(glt, organism = orgkg, keyType = 'kegg', pvalueCutoff = 1, pAdjustMethod = 'BH', use_internal_data = T, qvalueCutoff = 1)
   kg <- setReadable(kg, OrgDb=orgdb, keyType='ENTREZID')

   kg_topn <- kg %>%
        dplyr::filter(get(cutoff) < 0.05, Count >= count) %>%
        dplyr::arrange(get(cutoff)) %>%
        dplyr::top_n(-topshow, wt = get(cutoff))

   uterms <- get_term(kg, topshow, 'KEGG')
   obj <- filter(kg, kg@result$Description %in% uterms)
   heatplot_data <- get_heatplot_df(obj, uterms, geneList)

   return(
        list(
            kg_topn@result,
            heatplot_data
        )
   )
}


## reactome
enrich_rt <- function(glt) {
    rt <- enrichPathway(gene = glt, organism = react, pAdjustMethod = 'BH', pvalueCutoff = 1, qvalueCutoff = 1, readable = T)

    rt_topn <- rt %>%
        dplyr::filter(get(cutoff) < 0.05, Count >= count) %>%
        dplyr::arrange(get(cutoff)) %>%
        dplyr::top_n(-topshow, wt = get(cutoff))

    uterms <- get_term(rt, topshow, 'Reactome')
    obj <- filter(rt, rt@result$Description %in% uterms)
    heatplot_data <- get_heatplot_df(obj, uterms, geneList)

    return(
        list(
            rt_topn@result,
            heatplot_data
        )
    )
}

## wiki
enrich_wk <- function(glt) {
    wpgmt <- paste0('/opt/WikiPathways/', wiki_gmt, '.gmt')
    wp2gene <- read.gmt(wpgmt)
    wp2gene <- wp2gene %>% tidyr::separate('term', c("PathName", "version", "wpid", "org"), "%")
    wpid2gene <- wp2gene %>% dplyr::select(wpid, gene)
    wpid2name <- wp2gene %>% dplyr::select(wpid, PathName)
    wk <- enricher(gene = glt, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, pAdjustMethod = 'BH', pvalueCutoff = 1, qvalueCutoff = 1)
    wk <- setReadable(wk, OrgDb = orgdb, keyType = 'ENTREZID')

    wk_topn <- wk %>%
        dplyr::filter(get(cutoff) < 0.05, Count >= count) %>%
        dplyr::arrange(get(cutoff)) %>%
        dplyr::top_n(-topshow, wt = get(cutoff))

    uterms <- get_term(wk, topshow, 'Wiki')
    obj <- filter(wk, wk@result$Description %in% uterms)
    heatplot_data <- get_heatplot_df(obj, uterms, geneList)

    return(
        list(
            wk_topn@result,
            heatplot_data
        )
    )
}


## run
library(jsonlite)
args <- commandArgs(T)
inputs <- fromJSON(args[1])

species <- inputs$enrich_parameters$species
DE_file <- "DEGs_filter_result.tsv"
cutoff <- 'p.adjust'
topshow <- 10000
count <- 3

load_env()

genome <- data.frame(
    species_type = c('Homo sapiens', 'Mus musculus'),
    orgdb_package = c('org.Hs.eg.db', 'org.Mm.eg.db'),
    organisms = c('hsa', 'mmu'),
    organisms_rt = c('human', 'mouse'),
    organisms_gmt = c('Homo_sapiens', 'Mus_musculus')
)
orgdb <- genome$orgdb_package[genome$species_type == species]
orgkg <- genome$organisms[genome$species_type == species]
react <- genome$organisms_rt[genome$species_type == species]
wiki_gmt <- genome$organisms_gmt[genome$species_type == species]

logfoldchanges_threshold <- inputs$filter_parameters$logfoldchanges
pvals_threshold <- inputs$filter_parameters$pvals
gset <- list(name = read_genes(DE_file, logfoldchanges_threshold, pvals_threshold))

glist <- gset[['name']][['up']]
geneList <- gset[['name']][['upgeneList']]

glt_df <- bitr(glist, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = orgdb)
glt <- glt_df$ENTREZID

results_go <- enrich_go(glt)
write.table(results_go[1], file = "./enrich_go_result.tsv", col.names = T, row.names = F, quote = F, sep = '\t')

results_kg <- enrich_kg(glt)
write.table(results_kg[1], file = "./enrich_kg_result.tsv", col.names = T, row.names = F, quote = F, sep = '\t')

results_rt <- enrich_rt(glt)
write.table(results_rt[1], file = "./enrich_rt_result.tsv", col.names = T, row.names = F, quote = F, sep = '\t')

results_wk <- enrich_wk(glt)
write.table(results_wk[1], file = "./enrich_wk_result.tsv", col.names = T, row.names = F, quote = F, sep = '\t')


topshow <- 10
results_go <- enrich_go(glt)
write.table(results_go[2], file = "./heatplot_go_result.tsv", col.names = T, row.names = F, quote = F, sep = '\t')

results_kg <- enrich_kg(glt)
write.table(results_kg[2], file = "./heatplot_kg_result.tsv", col.names = T, row.names = F, quote = F, sep = '\t')

results_rt <- enrich_rt(glt)
write.table(results_rt[2], file = "./heatplot_rt_result.tsv", col.names = T, row.names = F, quote = F, sep = '\t')

results_wk <- enrich_wk(glt)
write.table(results_wk[2], file = "./heatplot_wk_result.tsv", col.names = T, row.names = F, quote = F, sep = '\t')
