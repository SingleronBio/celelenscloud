suppressMessages({
    library(argparser)
    library(GSEABase)
    library(GSVA)
    library(limma)
    library(pheatmap)
    library(svglite)
    library(RColorBrewer)
    library(jsonlite)
})


CreateDir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir, recursive = T)
    }
}


GetHeatmapColor <- function(hm_color, n=100) {
    col <- NULL
    if (hm_color == 'RedWhiteBlue') {
        col <- colorRampPalette(c("blue", "white", "red"))(n)
    }
    if (hm_color == 'RedWhiteBlue2') {
        col <- colorRampPalette(c("#0099CC", "white", "#CC0033"))(n)
    }
    if (hm_color == 'RedYellowBlue') {
        col <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(n)
    }
    return(col)
}


GetAnnotColors <- function(num_groups) {
    palette_name <- if (num_groups <= 8) "Set1" else if (num_groups <= 12) "Paired" else "Set3"
    colors <- brewer.pal(min(num_groups, length(brewer.pal.info[palette_name, "maxcolors"])), palette_name)
    return(colors)
}


ParseGmt <- function(gmt_file) {
    c3gsc2 <- getGmt(gmt_file,
                     collectionType = BroadCollection(category="c3"),
                     geneIdType = SymbolIdentifier())
    return(c3gsc2)
}


OrderGSVA <- function(gsvaOut) {
    fit <- as.data.frame(gsvaOut)
    fit$sd <- apply(fit, 1, sd)
    fit <- fit[order(-fit$sd), ]
    fit <- fit[, -which(colnames(fit) == "sd")]
    return(fit)
}


RunGSVA <- function(matrix, c3gsc2, outdir, prefix) {
    matrix <- as.matrix(matrix)
    dim <- list(rownames(matrix), colnames(matrix))
    matrix <- matrix(as.numeric(matrix), nrow = nrow(matrix), dimnames = dim)
    matrix <- avereps(matrix)
    matrix <- normalizeBetweenArrays(matrix)

    gsvaParm <- gsvaParam(matrix,
                          c3gsc2,
                          kcdf = "Gaussian")
    gsvaOut <- gsva(gsvaParm)
    gsvaOut <- OrderGSVA(gsvaOut = gsvaOut)
    write.table(data.frame(id = rownames(gsvaOut), gsvaOut, check.names = F),
                file=paste0(outdir, "/", prefix,".gsvaOut.csv"),
                sep=",", quote=F, row.names=F)
    return(gsvaOut)
}


DrawHeatMap <- function(gsvaOut, top, annotation_col, color = "RedWhiteBlue",
                        cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
                        show_rownames = TRUE, show_colnames = TRUE, outdir, prefix) {
    gsvaOut <- gsvaOut[1:top, ]

    hm_colors <- GetHeatmapColor(color)

    if (!is.null(annotation_col)) {
        unique_groups <- unique(annotation_col$Group)
        num_groups <- length(unique_groups)
        annot_colors_list <- GetAnnotColors(num_groups)
        annot_colors <- list(
            Group = setNames(annot_colors_list[1:num_groups], unique_groups)
        )
    } else {
        annot_colors <- NULL
    }

    svg_file <- file.path(outdir, paste0(prefix, ".gsvaOut.svg"))
    svglite(file=svg_file)
    pheatmap(gsvaOut,
             cluster_rows = cluster_rows,
             cluster_cols = cluster_cols,
             scale = scale,
             annotation_col = annotation_col,
             annotation_colors = annot_colors,
             show_rownames = show_rownames,
             show_colnames = show_colnames,
             color = hm_colors,
             border_color = "white"
    )
    dev.off()
}
