##library(Seurat)
##library(ComplexHeatmap)
##library(circlize)
##library(dplyr)
##library(grid)
##library(gtable)
##library(pheatmap)


##scale
scale_before_heat <- function(matrix, scale_by, center = TRUE, scale = TRUE){
    if(scale_by == "columns"){
        matrix <- scale(matrix, center = center, scale = scale)
    }
    if(scale_by == "rows"){
        matrix <- t(scale(t(matrix), center = center, scale = scale))
    }
    return(matrix)
}

##figure转换单位，如果是pdf()用unit = ”inch“，如果是png()用unit = ”pt“
calc_ht_size <- function(ht, unit = "inch") {
    pdf(NULL)
    ht = draw(ht)
    w = ComplexHeatmap:::width(ht)
    w = convertX(w, unit, valueOnly = TRUE)
    h = ComplexHeatmap:::height(ht)
    h = convertY(h, unit, valueOnly = TRUE)
    dev.off()
    print(c(w,h))
    return(c(w, h))
}

##ComplexHeatmap输出保存结果获取整个热图的长宽并转化成保存时需要的单位
save_heatmap <- function(ht, pdf = TRUE, png = TRUE, prefix, outdir){
    if(!dir.exists(outdir)){
        dir.create(outdir, recursive=T)
    }
    if(pdf){
        size_pdf <- calc_ht_size(ht)
        file1 <- paste0(outdir,'/', prefix, '.heatmap.pdf')
        #rownames_width1 <- convertX(ht@row_names_param[["max_width"]], "inches", valueOnly = TRUE)
        #colnames_height1 <- convertY(ht@column_names_param[["max_height"]], "inches", valueOnly = TRUE)
        #max1 <- max(rownames_width1, colnames_height1)
        #print(max1)
        rownames <- ht@row_names_param[["labels"]]
        rownames_width1 <- convertX(max_text_width(rownames, gp=gpar(fontsize = ht@row_names_param[["gp"]][["fontsize"]])), "inches", valueOnly = TRUE)
        colnames <- ht@column_names_param[["labels"]]
        colnames_height1 <- convertX(max_text_width(colnames, gp=gpar(fontsize = ht@column_names_param[["gp"]][["fontsize"]])), "inches", valueOnly = TRUE)
        max1 <- max(rownames_width1, colnames_height1)
        print(max1)
        pdf(file1, width = size_pdf[1] + max1, height = size_pdf[2] + max1 , bg="white")
        #pdf(file1, width = size_pdf[1], height = size_pdf[2], bg="white")
        draw(ht, merge_legend = TRUE)
        dev.off()
    }
    if(png){
        size_png <- calc_ht_size(ht, unit = "cm")
        file2 <- paste0(outdir,'/', prefix, '.heatmap.png')
        #rownames_width2 <- convertX(ht@row_names_param[["max_width"]], "cm", valueOnly = TRUE)
        #colnames_height2 <- convertY(ht@column_names_param[["max_height"]], "cm", valueOnly = TRUE)
        #max2 <- max(rownames_width2, colnames_height2)
        #print(max2)
        rownames <- ht@row_names_param[["labels"]]
        rownames_width2 <- convertX(max_text_width(rownames, gp=gpar(fontsize = ht@row_names_param[["gp"]][["fontsize"]])), "cm", valueOnly = TRUE)
        colnames <- ht@column_names_param[["labels"]]
        colnames_height2 <- convertX(max_text_width(colnames, gp=gpar(fontsize = ht@column_names_param[["gp"]][["fontsize"]])), "cm", valueOnly = TRUE)
        max2 <- max(rownames_width2, colnames_height2)
        print(max2)
        png(file2, width = size_png[1]+max2, height = size_png[2]+max2, bg="white", units = 'cm', res = 300)
        #png(file2, width = size_png[1]+2.5, height = size_png[2]+2.5, bg="white", units = 'cm', res = 300)
        draw(ht, merge_legend = TRUE)
        #draw(ht)
        dev.off()
    }
}

##pheatmap保存图片
save_pheatmap <- function(ht, pdf = TRUE, png = TRUE, prefix, outdir){
    if(!dir.exists(outdir)){
        dir.create(outdir, recursive=T)
    }
    if(pdf){
        width <- convertX(unit(sum(ht[["gtable"]][["widths"]])), "inches", valueOnly = TRUE) + 0.5
        height <- convertY(unit(sum(ht[["gtable"]][["heights"]])), "inches", valueOnly = TRUE) + 0.5
        file1 <- paste0(outdir,'/', prefix, '.heatmap.pdf')
        pdf(file1, width = width ,height = height, bg="white")
        print(ht)
        dev.off()
    }
    if(png){
        width <- convertX(unit(sum(ht[["gtable"]][["widths"]])), "inches", valueOnly = TRUE) + 0.5
        height <- convertY(unit(sum(ht[["gtable"]][["heights"]])), "inches", valueOnly = TRUE) + 0.5
        file2 <- paste0(outdir,'/', prefix, '.heatmap.png')
        png(file2, width = width ,height = height, bg="white", units = 'in', res = 300)
        print(ht)
        dev.off()
    }
}

## matrix为热图主体的矩阵数据；
## heat_col是热图主体颜色的色板（colorRampPalette函数形式）；
## cluster_rows和cluster_columns是否对行列进行聚类；
## row_split和column_split是否对行列进行分割，如果没有注释则不分割；如果有注释，则根据注释信息进行分割；如果没有注释也想要分割，需要提供按照行列的顺序提供一个类别向量，长度需要与行列相同；
## rowAnnotation和colAnnotation分别为两个dataframe，两个df的行名分别对应矩阵的行和列（顺序一致），列为对应的注释信息，列名为注释名称亦为图例标题；
## rowAnnotation_col与colAnnotation_col分别为行、列注释提供色板，并且需要赋names
## show_rowAnnotation_name是否展示行注释名（注意这里不是注释）

draw_ComplexHeatmap <- function(matrix, scale_or_not = TRUE, scale_by = "columns", center = TRUE, scale = TRUE, heat_col = NULL,
                         cellwidth = 18, cellheight = 12, heatmap_width = 20, heatmap_height = 20,
                         cluster_rows = FALSE, cluster_columns = FALSE, show_parent_dend_line = TRUE,
                         cluster_row_slices = TRUE, cluster_column_slices = TRUE,
                         row_dend_side = "left", row_dend_width = 10, show_row_dend = TRUE,
                         column_dend_side = "top", column_dend_height = 10, show_column_dend = TRUE,
                         row_split = NULL, column_split = NULL, row_split_by_rowAnnotation = 1, col_split_by_colAnnotation = 1,
                         row_gap = 2, column_gap = 2,
                         row_km = 1, row_km_repeats = 100, column_km = 1, column_km_repeats = 100,
                         na_col = "grey", border = TRUE, border_outside_col = "black", border_inside_col = "white",
                         row_title = NULL, row_title_side = "left", row_title_fontsize = 13.2, row_title_rot = 90, row_title_fontface = "bold",
                         column_title = NULL, column_title_side = "top", column_title_fontsize = 13.2, column_title_rot = 0, column_title_fontface = "bold",
                         show_row_names = TRUE, row_names_side = "right", row_names_fontsize = 10, row_names_fontface = "bold", row_names_rot = 0,
                         show_column_names = TRUE, column_names_side = "bottom", column_names_fontsize = 10, column_names_fontface = "bold", column_names_rot = 90,
                         heat_legend_title_fontsize = 10, heat_legend_title_fontface = "bold", heat_legend_title = "",
                         heat_legend_direction = "vertical", 
                         heat_legend_labels_fontsize = 10, heat_legend_labels_fontface = "bold",
                         rowAnnotation = NULL, colAnnotation = NULL, rowAnnotation_col = NULL, colAnnotation_col = NULL,
                         show_rowAnnotation_name = FALSE, rowAnnotation_name_side = "top", rowAnnotation_name_rot = 0, rowAnnotation_name_fontsize = 10, rowAnnotation_name_fontface = "bold",
                         rowAnnotation_inside_color = NA, rowAnnotation_outside_border = TRUE, 
                         show_rowAnnotation_legend = TRUE, rowAnnotation_legend_title_fontsize = 10, rowAnnotation_legend_title_fontface = "bold", rowAnnotation_legend_direction = "vertical", 
                         rowAnnotation_legend_labels_fontsize = 10, rowAnnotation_legend_labels_font = 2, rowAnnotation_legend_order = NULL,
                         show_colAnnotation_name = FALSE, colAnnotation_name_side = "right", colAnnotation_name_rot = 0, colAnnotation_name_fontsize = 10, colAnnotation_name_fontface = "bold",
                         colAnnotation_inside_color = NA, colAnnotation_outside_border = TRUE,
                         show_colAnnotation_legend = TRUE, colAnnotation_legend_title_fontsize = 10, colAnnotation_legend_title_fontface = "bold", colAnnotation_legend_direction = "vertical", 
                         colAnnotation_legend_labels_fontsize = 10, colAnnotation_legend_labels_font = 2, colAnnotation_legend_order = NULL,
                         display_numbers = FALSE, number_color = "black", number_format = "%.2f", number_fontsize = 6,
                         display_sig = FALSE, sig = NULL
                         ){
                            ann_colors = c("#E76F51","#FFAFCC","#264653","#0077B6","#DDBEA9","#81B29A","#00B4D8","#DC2F02","#FCA311","#57CC99")
                            clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
                            color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
                            clustcol <- c(ann_colors, color_protocol, clustcol)
                            if(scale_or_not){
                                matrix <- scale_before_heat(matrix, scale_by = scale_by, center = center, scale = scale)
                            }
                            if(is.null(heat_col)){
                                heatmap_col <- c("#0099CC", "white", "#CC0033")
                                if(nrow(matrix) * ncol(matrix) == 1000){
                                    col_fun <- colorRampPalette(heatmap_col)(1000 + 1)
                                }else{
                                    col_fun <- colorRampPalette(heatmap_col)(1000)
                                }
                            }else{
                                col_fun <- heat_col
                            }
                            use_number <- 0
                            if(is.null(rowAnnotation_col) & !is.null(rowAnnotation)){
                                rowAnnotation_col <- list()
                                for(i in 1:length(colnames(rowAnnotation))){
                                    row_col_use <- clustcol[(use_number + 1):(length(unique(rowAnnotation[,i])) + use_number)]
                                    names(row_col_use) <- unique(rowAnnotation[,i])
                                    rowAnnotation_col[[i]] <- row_col_use
                                    use_number <- use_number + length(unique(rowAnnotation[,i]))
                                }
                                names(rowAnnotation_col) <- colnames(rowAnnotation)
                            }
                            if(is.null(colAnnotation_col) & !is.null(colAnnotation)){
                                colAnnotation_col <- list()
                                for(j in 1:length(colnames(colAnnotation))){
                                    col_col_use <- clustcol[(use_number + 1):(length(unique(colAnnotation[,j])) + use_number)]
                                    names(col_col_use) <- unique(colAnnotation[,j])
                                    colAnnotation_col[[j]] <- col_col_use
                                    use_number <- use_number + length(unique(colAnnotation[,j]))
                                }
                                names(colAnnotation_col) <- colnames(colAnnotation)
                            }
                            if(is.null(rowAnnotation_legend_order) & !is.null(rowAnnotation)){
                                rowAnnotation_legend_order <- list()
                                for(x in 1:length(colnames(rowAnnotation))){
                                    rowAnnotation_legend_order[[x]] <- unique(rowAnnotation[,x])
                                    }
                            }
                            if(is.null(colAnnotation_legend_order) & !is.null(colAnnotation)){
                                colAnnotation_legend_order <- list()
                                for(y in 1:length(colnames(colAnnotation))){
                                   colAnnotation_legend_order[[y]] <- unique(colAnnotation[,y])
                                   } 
                            }
                            if(!is.null(rowAnnotation)){
                                matrix <- matrix[match(rownames(rowAnnotation), rownames(matrix)),]
                                length_rowAnnotation <- c()
                                for(i in 1:ncol(rowAnnotation)){
                                    length_rowAnnotation <- c(length_rowAnnotation, length(unique(rowAnnotation[,i])))
                                    max_length_rowAnnotation <- max(length_rowAnnotation)
                                    if(nrow(matrix) > max_length_rowAnnotation){ncol_rowAnnotation_legend <- 1}else{ncol_rowAnnotation_legend <- ceiling(max_length_rowAnnotation/nrow(matrix))}
                                }
                                row_annotation <- rowAnnotation(df = rowAnnotation, col = rowAnnotation_col, gp = gpar(col = rowAnnotation_inside_color), border = rowAnnotation_outside_border,
                                                                annotation_name_gp = gpar(fontsize = rowAnnotation_name_fontsize, fontface = rowAnnotation_name_fontface),
                                                                show_annotation_name = show_rowAnnotation_name, annotation_name_side = rowAnnotation_name_side, annotation_name_rot = rowAnnotation_name_rot,
                                                                show_legend = show_rowAnnotation_legend,
                                                                annotation_legend_param = list(title_gp = gpar(fontsize = rowAnnotation_legend_title_fontsize, fontface = rowAnnotation_legend_title_fontface), ncol = ncol_rowAnnotation_legend,
                                                                                              direction = rowAnnotation_legend_direction, labels_gp = gpar(fontsize = rowAnnotation_legend_labels_fontsize)))
                                for(z in 1:length(colnames(rowAnnotation))){
                                    row_annotation@anno_list[[colnames(rowAnnotation)[z]]]@legend_param[["labels_gp"]][["font"]] <- as.integer(rowAnnotation_legend_labels_font)
                                    row_annotation@anno_list[[colnames(rowAnnotation)[z]]]@color_mapping@levels <- rowAnnotation_legend_order[[z]]
                                }  
                            }else{
                                row_annotation <- NULL
                            }
                            if(!is.null(colAnnotation)){
                                matrix <- matrix[,match(rownames(colAnnotation), colnames(matrix))]
								length_colAnnotation <- c()
								print(matrix)
								for(i in 1:ncol(colAnnotation)){
                                    length_colAnnotation <- c(length_colAnnotation, length(unique(colAnnotation[,i])))
                                    max_length_colAnnotation <- max(length_colAnnotation)
									print(nrow(matrix))
									print(max_length_colAnnotation)
                                    if(nrow(matrix) > max_length_colAnnotation){ncol_colAnnotation_legend <- 1}else{ncol_colAnnotation_legend <- ceiling(max_length_colAnnotation/nrow(matrix))}
                                }
                                col_annotation <- HeatmapAnnotation(df = colAnnotation, col = colAnnotation_col, gp = gpar(col = colAnnotation_inside_color), border = colAnnotation_outside_border,
                                                                annotation_name_gp = gpar(fontsize = colAnnotation_name_fontsize, fontface = colAnnotation_name_fontface),
                                                                show_annotation_name = show_colAnnotation_name, annotation_name_side = colAnnotation_name_side, annotation_name_rot = colAnnotation_name_rot,
                                                                show_legend = show_colAnnotation_legend,
                                                                annotation_legend_param = list(title_gp = gpar(fontsize = colAnnotation_legend_title_fontsize, fontface = colAnnotation_legend_title_fontface), ncol = ncol_colAnnotation_legend,
                                                                                              direction = colAnnotation_legend_direction, labels_gp = gpar(fontsize = colAnnotation_legend_labels_fontsize)))
                                for(a in 1:length(colnames(colAnnotation))){
                                    col_annotation@anno_list[[colnames(colAnnotation)[a]]]@legend_param[["labels_gp"]][["font"]] <- as.integer(colAnnotation_legend_labels_font)
                                    col_annotation@anno_list[[colnames(colAnnotation)[a]]]@color_mapping@levels <- colAnnotation_legend_order[[a]]
                                }
                            }else{
                                col_annotation <- NULL
                            }
                            if(is.null(row_split)){rowsplit <- NULL}
                            if(length(row_split) == 1){
                                if(row_split == "km"){
                                    rowsplit <- NULL
                                }else if(row_split == "rowAnnotation"){
                                    rowsplit <- factor(rowAnnotation[,row_split_by_rowAnnotation], levels = unique(rowAnnotation[,row_split_by_rowAnnotation]))
                                }else{rowsplit <- row_split}
                            }else if(length(row_split) > 1){rowsplit <- row_split}
                            if(is.null(column_split)){colsplit <- NULL}
                            if(length(column_split) == 1){
                                if(column_split == "km"){
                                    colsplit <- NULL
                                }else if(column_split == "colAnnotation"){
                                    colsplit <- factor(colAnnotation[,col_split_by_colAnnotation], levels = unique(colAnnotation[,col_split_by_colAnnotation]))
                                }else{colsplit <- column_split}
                            }else if(length(column_split) > 1){
                                colsplit <- column_split
                            }
                            if(display_numbers){
                                layer_fun = function(j, i, x, y, width, height, fill){
                                    grid.text(sprintf(number_format, pindex(matrix, i, j)), x, y,
                                        gp = gpar(col = number_color, fontsize = number_fontsize))
                                }
                            }else if(display_sig){
                                layer_fun = function(j, i, x, y, width, height, fill){
                                    grid.text(pindex(sig, i, j), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
                                }
                            }else{
                                layer_fun <- NULL
                            }
                            #if(display_sig){
                            #    cell_fun = function(j, i, x, y, width, height, fill){
                            #       grid.text(sig[i,j], x, y, gp = gpar(fontsize = 8, fontface = "bold"))
                            #   }
                            #}else{
                            #    cell_fun = NULL
                            #}
                            if(is.null(cellwidth)){
                                width <- NULL
                                border_inside_col <- NA
                                heatmap_width_use <- unit(heatmap_width,"cm")
                            }else{
                                width <- ncol(matrix)*unit(cellwidth, "pt")
                                heatmap_width_use <- unit(1, "npc")
                            }
                            if(is.null(cellheight)){
                                height <- NULL
                                border_inside_col <- NA
                                heatmap_height_use <- unit(heatmap_height,"cm")
                            }else{
                                height <- nrow(matrix)*unit(cellheight, "pt")
                                heatmap_height_use <- unit(1, "npc")
                            }
                            ht <- ComplexHeatmap::Heatmap(matrix, col = col_fun, na_col = na_col, border = border, border_gp = gpar(col = border_outside_col), rect_gp = gpar(col = border_inside_col),
                                                          row_title = row_title, row_title_side = row_title_side, row_title_gp = gpar(fontsize = row_title_fontsize, fontface = row_title_fontface), row_title_rot = row_title_rot,
                                                          column_title = column_title, column_title_side = column_title_side, column_title_gp = gpar(fontsize = column_title_fontsize, fontface = column_title_fontface), column_title_rot = column_title_rot,
                                                          show_row_names = show_row_names, row_names_side = row_names_side, row_names_gp = gpar(fontsize = row_names_fontsize, fontface = row_names_fontface), row_names_rot = row_names_rot,
                                                          show_column_names = show_column_names, column_names_side = column_names_side, column_names_gp = gpar(fontsize = column_names_fontsize, fontface = column_names_fontface), column_names_rot = column_names_rot,
                                                          cluster_rows = cluster_rows, cluster_columns = cluster_columns, show_parent_dend_line = show_parent_dend_line,
                                                          row_split = rowsplit, column_split = colsplit, row_gap = unit(row_gap, 'mm'), column_gap = unit(row_gap, 'mm'),
                                                          row_km = row_km, row_km_repeats = row_km_repeats, column_km = column_km, column_km_repeats = column_km_repeats,
                                                          width = width, height = height, heatmap_width = heatmap_width_use, heatmap_height = heatmap_height_use,
                                                          left_annotation = row_annotation, top_annotation = col_annotation,
                                                          cluster_row_slices = cluster_row_slices, cluster_column_slices = cluster_column_slices,
                                                          row_dend_side = row_dend_side, row_dend_width = unit(row_dend_width, "mm"), show_row_dend = show_row_dend,
                                                          column_dend_side = column_dend_side, column_dend_height = unit(column_dend_height, "mm"), show_column_dend = show_column_dend,
                                                          row_names_max_width = max_text_width(rownames(matrix), gp = gpar(fontsize = row_names_fontsize, fontface = row_names_fontface)) * 1.3,
                                                          column_names_max_height = max_text_width(colnames(matrix),gp = gpar(fontsize = column_names_fontsize, fontface = column_names_fontface)) * 1.3,
                                                          heatmap_legend_param = list(title_gp = gpar(fontsize = heat_legend_title_fontsize, fontface = heat_legend_title_fontface),
                                                                                      title = heat_legend_title,
                                                                                      direction = heat_legend_direction,
                                                                                      labels_gp = gpar(fontsize = heat_legend_labels_fontsize, fontface = heat_legend_labels_fontface)),
                                                          layer_fun = layer_fun
                                                          )
                            if(show_row_names){ht@row_names_param[["anno"]]@width <- ht@row_names_param[["max_width"]]}
                            if(show_column_names){ht@column_names_param[["anno"]]@height <- ht@column_names_param[["max_height"]]}
                            return(ht)
                         }
                        
draw_pheatmap <- function(matrix, scale = "none", heat_col = NULL, kmeans_k = NA, breaks = NA, border_color = "white", cellwidth = 18, cellheight = 12,
                          cluster_rows = FALSE, cluster_cols = FALSE, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete",
                          cutree_rows = NA, cutree_cols = NA, treeheight_row = 50, treeheight_col = 50, legend = TRUE, legend_breaks = NA, legend_labels = NA,
                          annotation_row = NULL, annotation_col = NULL, annotation_colors = NULL, annotation_legend = TRUE, annotation_names_row = FALSE, annotation_names_col = TRUE,
                          show_rownames = TRUE, show_colnames = TRUE, main = NA, fontsize = 10, fontsize_row = 10, fontsize_col = 10, 
                          angle_col = "90", display_numbers = FALSE, number_format = "%.2f", number_color = "black", fontsize_number = 6, gaps_row = NULL, gaps_col = NULL, row_split_by_rowAnnotation = 1, col_split_by_colAnnotation = 1,
                          labels_row = NULL, labels_col = NULL, na_col = "grey"){
                            ann_colors = c("#E76F51","#FFAFCC","#264653","#0077B6","#DDBEA9","#81B29A","#00B4D8","#DC2F02","#FCA311","#57CC99")
                            clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
                            color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
                            clustcol <- c(ann_colors, color_protocol, clustcol)
                            if(is.null(heat_col)){
                                heatmap_col <- c("#0099CC", "white", "#CC0033")
                                if(nrow(matrix) * ncol(matrix) == 1000){
                                    col_fun <- colorRampPalette(heatmap_col)(1000 + 1)
                                }else{
                                    col_fun <- colorRampPalette(heatmap_col)(1000)
                                }
                            }else{
                                col_fun <- heat_col
                            }
                            use_number <- 0
                            annotation_colors_list <- list()
                            if(is.null(annotation_colors) & !is.null(annotation_row)){
                                matrix <- matrix[match(rownames(annotation_row), rownames(matrix)),]
                                for(i in 1:length(colnames(annotation_row))){
                                    row_col_use <- clustcol[(use_number + 1):(length(unique(annotation_row[,i])) + use_number)]
                                    names(row_col_use) <- unique(annotation_row[,i])
                                    annotation_colors_list[[i]] <- row_col_use
                                    use_number <- use_number + length(unique(annotation_row[,i]))
                                }
                                names(annotation_colors_list) <- colnames(annotation_row)
                            }
                            if(is.null(annotation_colors) & !is.null(annotation_col)){
                                matrix <- matrix[,match(rownames(annotation_col), colnames(matrix))]
                                for(j in 1:length(colnames(annotation_col))){
                                    col_col_use <- clustcol[(use_number + 1):(length(unique(annotation_col[,j])) + use_number)]
                                    names(col_col_use) <- unique(annotation_col[,j])
                                    annotation_colors_list[[j+length(colnames(annotation_row))]] <- col_col_use
                                    use_number <- use_number + length(unique(annotation_col[,j]))
                                }
                                names(annotation_colors_list) <- c(colnames(annotation_row), colnames(annotation_col))
                            }
                            if(is.null(annotation_row) & is.null(annotation_col)){annotation_colors_list <- NA}
                            if(!is.null(annotation_colors)){annotation_colors_list <- annotation_colors}
                            if(length(gaps_row) == 1){
                                if(gaps_row == "annotation"){
                                    fre_row <- c()
                                    for(a in unique(annotation_row[, row_split_by_rowAnnotation])){
                                        fre_row <- c(fre_row, as.numeric(table(annotation_row[, row_split_by_rowAnnotation])[a])) 
                                    }
                                    position_row<-c()
                                    for(x in 1:(length(fre_row)-1)){
                                        position_row[x]<-sum(fre_row[1:x])
                                    }
                                    gaps_row <- position_row
                                }else{gaps_row <- gaps_row}
                            }else{gaps_row <- gaps_row}
                            if(length(gaps_col) == 1){
                                if(gaps_col == "annotation"){
                                    fre_col <- c()
                                    for(b in unique(annotation_col[, col_split_by_colAnnotation])){
                                        fre_col <- c(fre_col, as.numeric(table(annotation_col[, col_split_by_colAnnotation])[b]))
                                    }
                                    position_col<-c()
                                    for(y in 1:(length(fre_col)-1)){
                                        position_col[y]<-sum(fre_col[1:y])
                                    }
                                    gaps_col <- position_col
                                }else{gaps_col <- gaps_col}
                            }else{gaps_col <- gaps_col}
                            if(is.null(annotation_row)){annotation_row <- NA}
                            if(is.null(annotation_col)){annotation_col <- NA}
                            
                            ht <- pheatmap(mat = matrix, scale = scale, color = col_fun, kmeans_k = kmeans_k, breaks = breaks, border_color = border_color, cellwidth = cellwidth, cellheight = cellheight,
                                           cluster_rows = cluster_rows, cluster_cols = cluster_cols, clustering_distance_rows = clustering_distance_rows, clustering_distance_cols = clustering_distance_cols,
                                           clustering_method = clustering_method, cutree_rows = cutree_rows, cutree_cols = cutree_cols, treeheight_row = treeheight_row, treeheight_col = treeheight_col,
                                           legend = legend, legend_breaks = legend_breaks, legend_labels = legend_labels,
                                           annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors_list, annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col,
                                           show_rownames = show_rownames, show_colnames = show_colnames, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
                                           angle_col = angle_col, display_numbers = display_numbers, number_format = number_format, number_color = number_color, fontsize_number = fontsize_number,
                                           gaps_row = gaps_row, gaps_col = gaps_col, labels_row = labels_row, labels_col = labels_col, na_col = na_col)
                            return(ht)
                          }
