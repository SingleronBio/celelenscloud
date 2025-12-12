source('/opt/scripts/monocle3_function.R')

## inputs
library(jsonlite)
args <- commandArgs(T)
inputs <- fromJSON(args[1])

# h5ad_file <- inputs$path
h5ad_file <- "./temp_file/result.h5ad"
annot_key <- inputs$runmonocle3_parameters$subset
subcelltype <- paste(inputs$runmonocle3_parameters$celltypes, collapse = ",")
subsample <- paste(inputs$runmonocle3_parameters$sampleid, collapse = ",")
reduction <- inputs$runmonocle3_parameters$reduction
if (is.null(reduction)) {
     reduction = FALSE
}
remove_batch <- inputs$runmonocle3_parameters$remove_batch
if (is.null(remove_batch)) {
     remove_batch = 'ignore'
}
knei <- inputs$runmonocle3_parameters$knei
if (is.null(knei)) {
     knei = 20
}

outdir <- './'

## run
result_list <- RunMonocle3(
    h5ad_file = h5ad_file,
    annot_key = annot_key,
    subcelltype = subcelltype,
    subsample = subsample,
    outdir = outdir,
    reduction = reduction,
    batch = remove_batch,
    knei = knei
)

write.table(result_list$nodelist_umap, file = 'nodelist_umap.xls', sep = "\t", quote = F, row.names = F)

