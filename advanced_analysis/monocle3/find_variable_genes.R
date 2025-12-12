source('/opt/scripts/monocle3_function.R')

## inputs
library(jsonlite)
args <- commandArgs(T)
inputs <- fromJSON(args[1])

monocle3file <- "monocle3.rds"
nodelist <- paste(inputs$find_variable_genes_parameters$nodelist, collapse = ",")
partitions <- paste(inputs$find_variable_genes_parameters$partitions, collapse = ",")
resolution <- inputs$find_variable_genes_parameters$resolution
q_value <- inputs$
    find_variable_genes_parameters$
    q_value_threshold

outdir <- './'

## run
result_list <- diffana(
    monocle3file = monocle3file,
    nodelist = nodelist,
    partitions = partitions,
    resolution = resolution,
    q_value = q_value
)
write.table(result_list$gene_result, file = 'gene_result.xls', sep = "\t", quote = F, row.names = F)
