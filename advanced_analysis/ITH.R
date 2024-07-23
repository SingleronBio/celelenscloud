ITH_dist <- function(rds, idents, type, ncell.min=20) {
    suppressPackageStartupMessages({
        library(Seurat)
        library(dplyr)
    })

    df = tibble(ITH_score=numeric(), type=character())
    Idents(rds) <- rds@meta.data[type] 
    
    for (ident in idents) {
        rds.ident <- subset(rds, ident=ident)
        ncell.ident <- length(colnames(rds.ident))

        if (ncell.ident < ncell.min) {
            print(paste(ident, "has", ncell.ident, "cells. Discarded."))
            next
        }

    tryCatch({

        rds.ident <- FindVariableFeatures(rds.ident, verbose=FALSE)
        rds.ident <- ScaleData(rds.ident, verbose=FALSE)
        rds.ident <- RunPCA(object=rds.ident, npcs=20, pc.genes=rds.ident@var.genes, verbose=FALSE)
        data <- (rds.ident[["pca"]]@cell.embeddings)
        DIST <- dist(data)
        dat <- as.vector(DIST)
        df.ident <- tibble(ITH_score=dat, type=ident)
        df <- rbind(df, df.ident)
        print(paste(ident, "done"))

    }, error=function(e) { print(e) })
  }
  return (df)
}
