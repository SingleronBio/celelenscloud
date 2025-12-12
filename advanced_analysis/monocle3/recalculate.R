source('/opt/scripts/monocle3_function.R')

## inputs
library(jsonlite)
args <- commandArgs(T)
inputs <- fromJSON(args[1])

monocle3file <- "monocle3.rds"
nodelist <- paste(inputs$recalculate_parameters$nodelist, collapse = ",")
knei <- inputs$recalculate_parameters$knei
if (is.null(knei)) {
     knei = 20
}

outdir <- './'

## run
Recal(
    monocle3file = monocle3file,
    nodelist = nodelist,
    knei = knei
)
