source('/opt/scripts/monocle3_function.R')

## inputs
library(jsonlite)
args <- commandArgs(T)
inputs <- fromJSON(args[1])

monocle3file <- "monocle3.rds"
genelist <- paste(inputs$genes_in_pseudotime_parameters$genelist, collapse = ",")

## run
genes_in_pseudotime(
    monocle3file = monocle3file,
    genelist = genelist
)
