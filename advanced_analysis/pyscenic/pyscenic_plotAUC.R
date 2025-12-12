args <- commandArgs(T)

suppressPackageStartupMessages({
library(argparser)
library(AUCell)
library(SCENIC)
library(reshape2)
library(Cairo)
library(grid)
library(pheatmap)
library(Seurat)
library(viridis)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(jsonlite)
library(svglite)
})

argv <- arg_parser('')
argv <- add_argument(argv, '--auc', help='the cell auc csv file')
argv <- add_argument(argv, '--tf', help='argvions, featureplot of tfs and tf regulon activity, eq: ATF3 ')
argv <- add_argument(argv, '--top', help='topN regulon file for plotting heatmap of each celltype')
argv <- add_argument(argv, '--h5', help='the h5 file')
argv <- add_argument(argv, '--metadata', help='the obs metadata file')
argv <- add_argument(argv, '--annot_order', help='the annot order file')
argv <- add_argument(argv, '--ident', help='cluster or sample', default='cluster')
argv <- add_argument(argv, '--split', help='sample or group', default='F')
argv <- add_argument(argv, '--outdir', help='', default='.')
argv <- add_argument(argv, '--prefix', help='', default='multi')

argv <- add_argument(argv,"--clusterrow", help="T|F", default='T')
argv <- add_argument(argv,"--clustercol", help="T|F", default='T')
argv <- parse_args(argv)

createDir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir)
    }
}

splitPRO <- function(PRO, split, ident){
    order <- levels(PRO)
    neworder <- unlist(lapply(order, function(x)paste(x, levels(PRO@meta.data[, split]), sep='_')))
    df <- data.frame(Idents(PRO), PRO@meta.data[, split])
    colnames(df) <- c(ident, split)
    df[,3] <- paste(df[,1], df[,2], sep='_')
    PRO@meta.data[[paste0(ident,'_', split)]] <- df[,3]
    Idents(PRO) <- df[,3]
    levels(PRO) <- neworder
    print(levels(PRO))
    return(PRO)
}


getPlotParams <- function(PRO, tag){
    if(tag == 'ptSize'){
        clust<-summary(PRO@active.ident)
        cluster_cell<-as.data.frame(clust)
        cell_number<-sum(cluster_cell$clust)
        pt_use<-0.6
        pt_use <- ifelse(cell_number > 6500, 0.1, ifelse(cell_number > 5500, 0.15, ifelse(cell_number > 4000, 0.2, ifelse(cell_number > 2500, 0.3, ifelse(cell_number > 1000, 0.4, 0.6)))))
        #print(paste('cell number: ', cell_number,';ptSize: ', pt_use, sep=''))
        return(pt_use)
    }
}

plotFunc <- function(pic, file){
    CairoPDF(file)
    print(pic)
    dev.off()

    file <- gsub('pdf$', 'png', file)
    CairoPNG(file)
    print(pic)
    dev.off()
}

zScore <- function(df){
    for(i in 1:ncol(df)){
        df[,i] <- as.numeric(df[,i])
        vector <- (df[,i]-mean(df[,i]))/sd(df[,i])
        if(i==1){
            z <- data.frame(vector)
        }else{
            z <- cbind(z, data.frame(vector))
        }
    }
    colnames(z) <- colnames(df)
    rownames(z) <- rownames(df)
    return(z)
}

heatmapCell <- function(argv, celltype, color_protocol){
    auc <- t(read.table(argv$auc, head=T, row.names=1, sep='\t', quote="", check.names=F))
    regulon <- read.table(argv$top, head=F, quote="")[,1]
    auc <- auc[,regulon]
    auc <- auc[rownames(celltype),]

    df <- as.matrix(zScore(auc))
    cluster_rows <- ifelse(argv$clusterrow!='F', T, F)
    cluster_columns <- ifelse(argv$clustercol!='F', T, F)

    annoColor <- color_protocol[1:length(unique(celltype[,1]))]
    names(annoColor) <- levels(celltype[,1])
    annoColor <- list(celltype=annoColor)
    annoCelltype <- HeatmapAnnotation(celltype=celltype[,1], col=annoColor, which='row', show_annotation_name=F, annotation_legend_param=list(title=NULL))
    # color <- c('#000080', 'white', '#B22222')
    # color <- brewer.pal(n = 9, name = "YlGnBu")
    color <- colorRamp2(c(-4, 0, 4), c('#000080', 'white', '#B22222'))

    p <- Heatmap(
        as.matrix(df),
        col = color,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        show_row_names = F,
        show_column_names = T,
        show_heatmap_legend = T,
        left_annotation = annoCelltype,
        heatmap_legend_param=list(title='Z-score'),
        use_raster=F)

    mar_bottom <- max(nchar(regulon))/1.8 + 4
    lengths <- max(nchar(as.character(unique(celltype[,1]))))
    width <- ifelse(length(regulon)<5, 3, length(regulon)*0.28)
    width <- max(c(3, length(regulon)*0.28))
    width <- width + 0.1 * lengths
    file <- paste0(file.path(argv$outdir, argv$prefix), "_heatmap_top.pdf")
    pdf(file, width=width, bg='white')
    par(mar=c(mar_bottom,1,1,1))
    draw(p)
    dev.off()
    file <- paste0(file.path(argv$outdir, argv$prefix), "_heatmap_top.png")
    png(file, width=width*80)
    par(mar=c(mar_bottom,1,1,1))
    draw(p)
    dev.off()
    file <- paste0(file.path(argv$outdir, argv$prefix), "_heatmap_top.svg")
    svglite(file, width=width, bg='white')
    par(mar=c(mar_bottom,1,1,1))
    draw(p)
    dev.off()
    df1 <- as.matrix(df)
    df1 <- ifelse(df1 > 4, 4, ifelse(df1 < -4, -4, df1))
    df1 <- cbind(rownames(df1), df1)
    write.table(df1, paste0(file.path(argv$outdir, argv$prefix), "_heatmap_top.xls"), sep="\t", row.names=F)
    annot <- data.frame(celltype[,1])
    row.names(annot) <- row.names(df)
    annot <- cbind(rownames(annot), annot)
    colnames(annot) <- c('barcode', 'annot')
    write.table(annot, paste0(file.path(argv$outdir, argv$prefix), "_heatmap_top_annotation.xls"), sep="\t", row.names=F)
}

plotAverageHeatmap <- function(df, argv){
    p <- pheatmap(df,
        cluster_rows = F,
        cluster_cols = F,
        scale='row',
        border_color='gray90')
        row <- nrow(df)
        col <- ncol(df)
        height <- ifelse(row <40, 7, (row-40)*7/40+7)
        width <- ifelse(col<5, 4, (col-5)/3+4)
    file <- paste0(argv$outdir,'/', argv$prefix,'.average_expression.pdf')
    pdf(file, width=width, height=height,bg='white')
    print(p)
    dev.off()
    png <- gsub('pdf$', 'png', file)
    system(paste('convert -density 450', file, png, sep=' '))
    file <- paste0(argv$outdir,'/', argv$prefix,'.average_expression.svg')
    svglite(file, width=width, height=height,bg='white')
    print(p)
    dev.off()
    df <- cbind(rownames(df), df)
    write.table(df, paste0(argv$outdir,'/', argv$prefix,'.average_expression.xls'), sep='\t', row.names=F)
}

main <- function(){
    createDir(argv$outdir)

    RAW <- Read10X_h5(argv$h5)
    PRO <- CreateSeuratObject(RAW)
    metadata <- read.table(argv$metadata, sep='\t', header=T, row.names=1, comment.char = "", check.names = FALSE)
    PRO@meta.data <- metadata

    json_data <- fromJSON(argv$annot_order)
    column_order <- names(json_data)
    for (col in column_order) {
         PRO@meta.data[[col]] = factor(PRO@meta.data[[col]], levels=json_data[[col]])
    }
    ident <- argv$ident
    if(argv$split!='F'){
        PRO <- splitPRO(PRO, argv$split, ident)
        ident <- paste0(ident,'_', argv$split)
    }else{
        Idents(PRO) <- PRO@meta.data[, ident]
        ident <- ident
    }
    print(levels(PRO))
    averageExp <- AverageExpression(PRO, verbose=F)$RNA

    auc <- read.table(argv$auc, head=T, row.names=1, sep='\t', quote="", check.names=F)
    PRO <- subset(PRO, cells=colnames(auc))

    df <- data.frame(celltype=sort(Idents(PRO)))
    rownames(df) <- names(sort(Idents(PRO)))
    ident_color_name <- paste0(ident, "_colors")
    if (!is.null(PRO[[ident_color_name]])) {
        celltype_colors <- data.frame(
             celltype = PRO[[ident]],
             cluster_colors = PRO[[ident_color_name]]
        )
        colnames(celltype_colors) <- c('celltype', 'cluster_colors')
        unique_celltype_colors <- celltype_colors %>% distinct(celltype, cluster_colors)
        unique_celltype_colors$celltype <- factor(unique_celltype_colors$celltype,
                                          levels = levels(PRO))
        sorted_colors <- unique_celltype_colors[order(unique_celltype_colors$celltype), ]
        color_protocol <- sorted_colors$cluster_colors
    }else{
        color_protocol <- c(
        "#0067AA", "#FF7F00", "#00A23F", "#FF1F1D", "#A763AC", "#B45B5D", "#FF8AB6", "#B6B800", "#01C1CC",
        "#85D5F8", "#FFC981", "#C8571B", "#C6CCC3", "#727272", "#EFC800", "#8A5626", "#502E91", "#59A4CE",
        "#344B2B", "#FBE29D", "#FDD6E6", "#849C8C", "#F07C6F", "#000101", "#FF4500", "#4F94CD", "#FF8C00",
        "#ADFF2F", "#A020F0", "#2F4F4F", "#FFD700", "#006400", "#FF1493", "#8B0000", "#4682B4", "#FFDAB9",
        "#708090", "#836FFF", "#CDC673", "#CD9B1D", "#FF6EB4", "#CDB5CD", "#008B8B", "#43CD80", "#483D8B",
        "#66CD00", "#CDC673", "#CDAD00", "#CD9B9B", "#FF8247", "#8B7355", "#8B3A62", "#68228B", "#CDB7B5",
        "#CD853F", "#6B8E23", "#696969", "#7B68EE", "#9F79EE", "#B0C4DE", "#7A378B", "#66CDAA", "#EEE8AA",
        "#00FF00", "#EEA2AD", "#A0522D", "#000080", "#E9967A", "#00CDCD", "#8B4500", "#DDA0DD", "#EE9572",
        "#EEE9E9", "#8B1A1A", "#8B8378", "#EE9A49", "#EECFA1", "#8B4726", "#8B8878", "#EEB4B4", "#C1CDCD",
        "#8B7500", "#0000FF", "#EEEED1", "#4F94CD", "#6E8B3D", "#B0E2FF", "#76EE00", "#A2B5CD", "#548B54",
        "#BBFFFF", "#B4EEB4", "#00C5CD", "#008B8B", "#7FFFD4", "#8EE5EE", "#43CD80", "#68838B", "#00FF00",
        "#B9D3EE", "#9ACD32", "#00688B", "#FFEC8B", "#1C86EE", "#CDCD00", "#473C8B", "#FFB90F", "#EED5D2",
        "#CD5555", "#CDC9A5", "#FFE7BA", "#FFDAB9", "#CD661D", "#CDC5BF", "#FF8C69", "#8A2BE2", "#CD8500",
        "#B03060", "#FF6347", "#FF7F50", "#CD0000", "#F4A460", "#FFB5C5", "#DAA520", "#CD6889", "#32CD32",
        "#FF00FF", "#2E8B57", "#CD96CD", "#48D1CC", "#9B30FF", "#1E90FF", "#CDB5CD", "#191970", "#E8E8E8"
        )
    }
    heatmapCell(argv, df, color_protocol)

    regulonAUC <- importAUCfromText(argv$auc)
    print(head(regulonAUC))
    regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

    if(!is.na(argv$tf)){
        selected_tf <- unlist(strsplit(argv$tf, ','))
        selected_regulon <- c()
        filter_tf <- c()
        for(t in selected_tf){
            r <- colnames(regulonAUC)[grepl(paste0('^', t, '\\('), colnames(regulonAUC))]
            if(length(r)==0){
                print(paste0("TF not exits in AUC file: ", t))
                next
            }
            print(r)
            filter_tf <- c(filter_tf, t)
            regulonData <- getAUC(regulonAUC)[, r]
            names(regulonData) <- rownames(regulonAUC)
            PRO <- AddMetaData(object = PRO, metadata = regulonData, col.name = paste0(t, '_regulon'))
            selected_regulon <- c(selected_regulon, paste0(t, '_regulon'))
        }
        if(length(selected_regulon)>0){
            pt_use <- getPlotParams(PRO, 'ptSize')
            for(feature in c(filter_tf, selected_regulon)){
                p <- FeaturePlot(PRO,features=feature,cols = c("lightgrey","red"),min.cutoff="q1",max.cutoff="q99",pt.size = pt_use,reduction="umap")
                p <- p + scale_color_viridis()
                file <- file.path(argv$outdir, paste0(argv$prefix, '_', feature, '_featureplot.pdf'))
                plotFunc(p, file)
            }
        }
    }

    if(!is.na(argv$top)){
        selected_regulon <- read.table(argv$top, head=F)[, 1]
        print(selected_regulon)

        celltype <- data.frame(cell=colnames(PRO), celltype=Idents(PRO))
        regulonActivity_byCellType <- sapply(split(celltype[, 1], celltype[, 2]),
                                        function(cells){
                                            if(length(cells)==1){
                                                getAUC(regulonAUC)[cells, ]
                                            }else{
                                                colMeans(getAUC(regulonAUC)[cells, ])
                                            }
                                        })

        center=T
        if(ncol(regulonActivity_byCellType)==2){
            center=F
        }
        regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = center, scale=T))
        df <- regulonActivity_byCellType_Scaled[selected_regulon, ]
        row <- nrow(df)
        col <- ncol(df)
        height <- ifelse(row <40, 7, (row-40)*7/40+7)
        width <- ifelse(col<5, 4, (col-5)/3+4)

        file <- file.path(argv$outdir, paste0(argv$prefix, "_regulonActivity_CellTypeHeatmap.pdf"))
        pdf(file, width=width, height=height, bg='white')
        pheatmap::pheatmap(regulonActivity_byCellType_Scaled[selected_regulon, ],
            show_rownames = T,
            show_colnames = T,
            cluster_rows = T,
            cluster_cols = F,
            angle_col="90",
            color=colorRampPalette(c('#000080', 'white', '#B22222'))(100),
            treeheight_row=10, treeheight_col=10,
            border_color='gray90')
        dev.off()
        pngfile <- sub('pdf$', 'png', file)
        system(paste0('convert ', file, ' ', pngfile))

        file <- file.path(argv$outdir, paste0(argv$prefix, "_regulonActivity_CellTypeHeatmap.svg"))
        svglite(file, width=width, height=height, bg='white')
        pheatmap::pheatmap(regulonActivity_byCellType_Scaled[selected_regulon, ],
            show_rownames = T,
            show_colnames = T,
            cluster_rows = T,
            cluster_cols = F,
            angle_col="90",
            color=colorRampPalette(c('#000080', 'white', '#B22222'))(100),
            treeheight_row=10, treeheight_col=10,
            border_color='gray90')
        dev.off()

	ptgene<-pheatmap::pheatmap(regulonActivity_byCellType_Scaled[selected_regulon, ],
            show_rownames = T,
            show_colnames = T,
            cluster_rows = T,
            cluster_cols = F,
            color=colorRampPalette(c('#000080', 'white', '#B22222'))(100),
            treeheight_row=10, treeheight_col=10,
            border_color='gray90')
        gorder<-regulonActivity_byCellType_Scaled[selected_regulon, ][ptgene$tree_row$order,]
        x <- cbind(rownames(regulonActivity_byCellType), regulonActivity_byCellType)
        file <- file.path(argv$outdir, paste0(argv$prefix, "_regulonActivity_CellType.xls"))
        write.table(x, file=file, sep="\t", row.names=F)

        scaled_x <- regulonActivity_byCellType_Scaled[selected_regulon, ]
        scaled_x <- cbind(rownames(scaled_x), scaled_x)
        file_top <- file.path(argv$outdir, paste0(argv$prefix, "_regulonActivity_byCellType_Scaled_top.xls"))
        write.table(scaled_x, file=file_top, sep="\t", row.names=F)

        gene <- as.character(gsub('\\(.*', '', rownames(gorder)))
	print(gene)
	df <- AverageExpression(PRO,verbose=F,features =gene)$RNA
	print(head(df))
        plotAverageHeatmap(df,argv)
    }
}
main()
