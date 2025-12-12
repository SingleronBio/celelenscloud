################################subfunction################################################
load_env <- function() {
    suppressMessages({
        library(ggplot2)
        library(reshape2)
        library(Cairo)
        library(Seurat)
        library(dplyr)
        library(grid)
        library(cowplot)
        library(monocle3)
        library(viridis)
        library(stringr)
        library(jsonlite)
        library(data.table)
        library(rhdf5)
        library(RColorBrewer)
        library(svglite)
        library(Matrix)
        library(harmony)
    })

    set.seed(1122334455)

    #### 调色板 ####
    col1 <- colorRampPalette(c("#7F0000", "red", "red", "#FF7F00", "#FF7F00", "yellow", "yellow", "cyan", "#007FFF", "blue", "#00007F"))
    corrcol <- colorRampPalette(c("red", "orange", "blue", "white", "white"))
    clustcol <- c("#0067AA", "#FF7F00", "#00A23F", "#FF1F1D", "#A763AC", "#B45B5D", "#FF8AB6", "#B6B800", "#01C1CC", "#85D5F8", "#FFC981", "#C8571B", "#C6CCC3", "#727272", "#EFC800", "#8A5626", "#502E91", "#59A4CE", "#344B2B", "#FBE29D", "#FDD6E6", "#849C8C", "#F07C6F", "#000101", "OrangeRed", "SlateBlue3", "DarkOrange", "GreenYellow", "Purple", "DarkSlateGray", "Gold", "DarkGreen", "DeepPink2", "Red4", "#4682B4", "#FFDAB9", "#708090", "#836FFF", "#CDC673", "#CD9B1D", "#FF6EB4", "#CDB5CD", "#008B8B", "#43CD80", "#483D8B", "#66CD00", "#CDC673", "#CDAD00", "#CD9B9B", "#FF8247", "#8B7355", "#8B3A62", "#68228B", "#CDB7B5", "#CD853F", "#6B8E23", "#E6E6FA", "#FFDAB9")
}

#### MODULE 1: 预处理原始seurat对象数据 ####

#' h5ad转为rds
to_rds <- function(h5ad_file = NULL, annot_key = NULL, batch = NULL, layer_name = '/layers/filtered') {
    layers_all <- h5read(file = h5ad_file, name = "/layers")
    if(! 'filtered' %in% names(layers_all) ){layer_name = '/layers/abnormal_filtered'}

    counts_layer <- h5read(file = h5ad_file, name = layer_name)
    var_genes <- h5read(file = h5ad_file, name = "/var")$'_index'
    obs_cells <- h5read(file = h5ad_file, name = "/obs")$'_index'

    counts_attr <- h5readAttributes(file = h5ad_file, name = layer_name)
    var_length <- counts_attr$shape[2]
    obs_length <- counts_attr$shape[1]

    mtx <- try(sparseMatrix(
        i = as.integer(counts_layer$indices),
        p = as.integer(counts_layer$indptr),
        x = as.numeric(counts_layer$data),
        dims = c(var_length, obs_length),
        index1 = F,
        repr = "C"
    ))
    if ("try-error" %in% class(mtx)){
        mtx <- t(sparseMatrix(i = as.integer(counts_layer$indices),
                      p = as.integer(counts_layer$indptr),
                      x = as.numeric(counts_layer$data),
                      dims = c(length(obs_cells), length(var_genes)),
                      index1 = F,
                      repr = "C"))
    }

    rownames(mtx) <- var_genes
    colnames(mtx) <- obs_cells

    PRO <- CreateSeuratObject(counts = mtx, project = "anno", min.cells = 0)
    PRO <- NormalizeData(PRO)

    sp_ind <- h5read(h5ad_file, name = "/obs/Sample ID")
    sp_str <- h5read(h5ad_file, name = "/obs/__categories/Sample ID")
    sp_ind_1 <- as.numeric(sp_ind) + 1
    sampleID <- sp_str[sp_ind_1]
    PRO$Sample.ID <- sampleID

    tryCatch({sp_ind_color <- h5read(h5ad_file, name = paste0("/obs/", "Sample ID", "_colors"))
              sp_str_color <- h5read(h5ad_file,name = paste0("/obs/__categories/", "Sample ID", "_colors"))
              sp_ind_color_1 <- as.numeric(sp_ind_color) + 1
              sample_color <- sp_str_color[sp_ind_color_1]
              PRO$sample_colors <- sample_color}
    ,error=function(e){print("NO SAMPLE COLOR")})

    sp_ind <- h5read(h5ad_file, name = paste0("/obs/", annot_key))
    sp_str <- h5read(h5ad_file, name = paste0("/obs/__categories/", annot_key))
    sp_ind_1 <- as.numeric(sp_ind) + 1
    annot_info <- sp_str[sp_ind_1]
    PRO$celltype <- annot_info

    tryCatch({sp_ind_color <- h5read(h5ad_file, name = paste0("/obs/", annot_key, "_colors"))
              sp_str_color <- h5read(h5ad_file,name = paste0("/obs/__categories/", annot_key, "_colors"))
              sp_ind_color_1 <- as.numeric(sp_ind_color) + 1
              cluster_color <- sp_str_color[sp_ind_color_1]
              PRO$cluster_colors <- cluster_color}
    ,error=function(e){print("NO CLUSTER COLOR")})

    if (batch != 'ignore') {
        g_ind  <- h5read(h5ad_file,name = paste0("/obs/", batch))
        g_str  <- h5read(h5ad_file,name = paste0("/obs/__categories/", batch))
        g_ind_1 <- as.numeric(g_ind) + 1
        group  <- g_str[g_ind_1]
        PRO$group <- group
        PRO$group <- factor(PRO$group, levels=g_str)
    }

    umapinfo <- as.data.frame(t(h5read(file = h5ad_file, name = '/obsm/X_umap')))
    rownames(umapinfo) <- colnames(PRO)
    colnames(umapinfo) <- c('UMAP_1', 'UMAP_2')

    PRO[['umap']] <- CreateDimReducObject(embeddings = as.matrix(umapinfo), key = 'UMAP_', assay = 'RNA', global = TRUE)

    return(PRO)
}

#' 获取subcelltype与subsample
GetTarget <- function(X, subcelltype, subsample) {
    X <- subset(X, Sample.ID %in% subsample)
    X <- subset(X, celltype %in% subcelltype)

    X
}

#### re-dimensionality reduction
re_reduction <- function(PRO = NULL, batch='ignore') {
    if (batch != 'ignore') {
        PRO <- NormalizeData(PRO)
        PRO <- FindVariableFeatures(PRO, selection.method = "vst", nfeatures = 2000)
        PRO <- ScaleData(PRO)
        PRO <- RunPCA(PRO, npcs = 20, verbose = FALSE)
        PRO <- RunHarmony(PRO, group.by.vars = 'group', plot_convergence = FALSE)
        PRO <- RunUMAP(PRO, reduction = "harmony", dims = 1:20)
    } else {
        PRO <- NormalizeData(PRO)
        PRO <- FindVariableFeatures(PRO, selection.method = "vst", nfeatures = 2000)
        PRO <- ScaleData(PRO)
        PRO <- RunPCA(PRO, npcs = 20, verbose = FALSE)
        PRO <- RunUMAP(PRO, reduction = "pca", dims = 1:20)
    }
    return(PRO)
}

#### MODULE 2: 建立monocle 3的CDS对象 ####

#' 构建monocle3的数据对象，根据提供的seurat数据对象进行构建
# PRO_seurat <- importCDS(PRO)
#' @param object, seurat数据对象
customed_createCDS <- function(PRO = NULL) {
    # pfile: 存储meta信息
    # @param PRO seurat data object
    PRO_meta_info <- PRO@meta.data[, setdiff(colnames(PRO@meta.data), 'seurat_clusters'), drop = F]
    # Cluster <- PRO$seurat_clusters
    barcode <- rownames(PRO_meta_info)
    # 固定最后两列列名
    pfile <- data.frame(barcode = barcode, PRO_meta_info, Size_Factor = NA, stringsAsFactors = FALSE)
    cat("meta info in original Seurat object\n", paste(colnames(pfile), collapse = ", "), "\n")

    # @param tmp_unique_cluster 用于热图函数 plot_pseudotime_heatmap，为参数num_clusters赋值
    # 取值角度：理想中，基因表达模式的聚类数目与数据中的细胞类型数目一致，即每种细胞类型都有特异的基因表达模式。
    # fdata：存储基因信息
    geneid <- rownames(PRO)
    fdata <- data.frame(id = geneid, gene_short_name = geneid, row.names = geneid)
    head(fdata)

    # matrix：使用原始表达的counts
    count <- GetAssayData(object = PRO[["RNA"]], slot = "counts")
    head(rownames(count))
    if (!identical(rownames(fdata), rownames(count))) stop(" 基因名称不一致 ")

    # 构建monocle3的CDS对象
    my_cds <- new_cell_data_set(
        expression_data = count,
        cell_metadata = pfile,
        gene_metadata = fdata
    )
    my_cds@int_colData$reducedDims$UMAP <- Embeddings(object = PRO, reduction = "umap")
    my_cds <- preprocess_cds(my_cds)
    cat("CDS数据对象构建完成", "\n")
    return(my_cds)
}

#### MODULE 3: 轨迹分析 ####

#' 利用umap维度进行图学习，建立轨迹
out_plot <- function(plot = NULL, filename = NULL, width = 7, height = 7) {
    CairoPDF(file = paste0(filename, ".pdf")); print(plot); dev.off()
    CairoPNG(file = paste0(filename, ".png")); print(plot); dev.off()
}

#' 取根节点：
#' 在官方参考示例中，作者的cds数据是胚胎的，且有每个细胞的胚胎发育时间，因此，以规定好的发育时间为指标，
#' 作者选取了 胚胎早期(130-170) 的所有细胞中，占比最大的细胞团作为起点
#' 在我们自己的代码中，则是从细胞类型出发，取每个细胞类型中细胞细胞数最多细胞团的编号，作为（该细胞类型的）轨迹起点
#' 这样做的好处在于，通过给不同细胞类型指定轨迹起点，可以一次分别得到多个细胞类型各自的伪时间分化轨迹
#' @param meta_col 选取那一列meta信息作为判断轨迹起始，默认是celltype，如果有其他信息也可（比如参考示例中的 实际胚胎发育时间）
#' @param supposed_start 在指定作为判断指标的meta_col中，哪些细胞可能是轨迹起点，这些细胞中占比最大的细胞团（以Y_xx命名），就是起点
get_earliest_principal_node <- function(cds, supposed_start = NULL, meta_col = "celltype") {
    cat("以", meta_col, "作为判断轨迹的指标\n")
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds),])
    root_pr_nodes <- c()
    for (tmp_start in supposed_start) {
        # 该细胞类型中占比最大的细胞团 Y_xx，设定为轨迹起点，并返回该起点细胞团名称，用于order_cells分析时的起点。
        cell_ids <- which(colData(cds)[, meta_col] == tmp_start)
        tmp_root_pr_node <-
            igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
        cat("以", tmp_start, "为假定起点，则细胞团", tmp_root_pr_node, "为轨迹起始节点细胞团\n")
        root_pr_nodes <- c(root_pr_nodes, tmp_root_pr_node)
    }
    return(root_pr_nodes)
}

#' 轨迹分析
trajectory_ana <- function(cds = NULL, outdir = NULL, prefix = NULL, set_root = FALSE, nodelist = NULL, knei = 20) {

    clustcol <- c("#0067AA", "#FF7F00", "#00A23F", "#FF1F1D", "#A763AC", "#B45B5D", "#FF8AB6", "#B6B800", "#01C1CC", "#85D5F8", "#FFC981", "#C8571B", "#C6CCC3", "#727272", "#EFC800", "#8A5626", "#502E91", "#59A4CE", "#344B2B", "#FBE29D", "#FDD6E6", "#849C8C", "#F07C6F", "#000101", "OrangeRed", "SlateBlue3", "DarkOrange", "GreenYellow", "Purple", "DarkSlateGray", "Gold", "DarkGreen", "DeepPink2", "Red4", "#4682B4", "#FFDAB9", "#708090", "#836FFF", "#CDC673", "#CD9B1D", "#FF6EB4", "#CDB5CD", "#008B8B", "#43CD80", "#483D8B", "#66CD00", "#CDC673", "#CDAD00", "#CD9B9B", "#FF8247", "#8B7355", "#8B3A62", "#68228B", "#CDB7B5", "#CD853F", "#6B8E23", "#E6E6FA", "#FFDAB9")

    if (!is.null(cds$cluster_colors)) {
        celltype_colors <- data.frame(
             celltype = cds$celltype,
             cluster_colors = cds$cluster_colors
        )
        unique_celltype_colors <- celltype_colors %>% distinct(celltype, cluster_colors)
        unique_celltype_colors$celltype <- factor(unique_celltype_colors$celltype,
                                          levels = levels(factor(cds$celltype)))
        sorted_colors <- unique_celltype_colors[order(unique_celltype_colors$celltype), ]
        celltypecol <- sorted_colors$cluster_colors
    }else{
        celltypecol <- clustcol
    }

    if (!set_root) {
        cat("未指定轨迹根节点，程序讲自动为您选择根节点进行分析")
        p <- plot_cells(cds, label_groups_by_cluster = TRUE, color_cells_by = "celltype") + scale_color_manual(values = clustcol)

        out_plot(plot = p, filename = paste0(outdir, "/", prefix, "_umap_celltype"))
        # 运行cluster_cells()，每个cell不仅被分配给一个cluster，还被分配给一个partition。
        # 当学习轨迹时，每一个partition最终会成为一个独立的trajactory。partition有点儿类似于monocle2中的state
        # 计算返回 clusters，partitions， 储存在cds@@clusters@listData$UMAP$clusters, cds@@clusters@listData$UMAP$partitions 中
        cds <- cluster_cells(cds, k=knei)
        cds@clusters@listData$UMAP$partitions <- as.factor(cds@clusters@listData$UMAP$partitions)
        cds <- learn_graph(cds)
        write.table(sort(unique(cds@clusters@listData$UMAP$partitions)), file = paste(outdir, "/", prefix, ".partitions.xls", sep = ""), quote = F, sep = "\t", row.names = F, col.names = 'partitions')
        for (color_by_meta_col in c("partition", "celltype")) {
            if (color_by_meta_col == "partition") {
                p <- plot_cells(cds, color_cells_by = color_by_meta_col, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE)
                p <- p + scale_color_manual(values = clustcol)
                out_plot(plot = p, filename = paste0(outdir, "/", prefix, "_trajectory_", color_by_meta_col))
                p <- p +
                    theme(
                        panel.background = element_rect(fill = "transparent"), # bg of the panel
                        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                        panel.grid.major = element_blank(), # get rid of major grid
                        panel.grid.minor = element_blank(), # get rid of minor grid
                        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                        legend.box.background = element_rect(color = NA) # get rid of legend panel bg
                    )

                svglite(paste0(outdir, "/partition.svg"), width = 10.5, height = 7)
                print(p)
                dev.off()
            }
            if (color_by_meta_col == "celltype") {
                p <- plot_cells(cds, color_cells_by = color_by_meta_col, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE)
                p <- p + scale_color_manual(values = clustcol)
                out_plot(plot = p, filename = paste0(outdir, "/", prefix, "_trajectory_", color_by_meta_col))

                p <- plot_cells(cds, color_cells_by = color_by_meta_col, label_groups_by_cluster = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 2.6, label_cell_groups = FALSE)
                p <- p + scale_color_manual(values = celltypecol)
                out_plot(plot = p, filename = paste0(outdir, "/", prefix, "_trajectory_node_", color_by_meta_col))
                p <- p +
                    theme(
                        panel.background = element_rect(fill = "transparent"), # bg of the panel
                        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                        panel.grid.major = element_blank(), # get rid of major grid
                        panel.grid.minor = element_blank(), # get rid of minor grid
                        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                        legend.box.background = element_rect(color = NA) # get rid of legend panel bg
                    )

                svglite(paste0(outdir, "/celltype.svg"), width = 10.5, height = 7)
                print(p)
                dev.off()
            }
        }
        # get all nodes
        lline <- t(cds@principal_graph_aux$UMAP$dp_mst)
        lline_m <- data.frame(ID = rownames(lline), UMAP_1 = lline[, "UMAP_1"], UMAP_2 = lline[, "UMAP_2"], stringsAsFactors = F)
        write.table(lline_m, file = paste(outdir, "/", prefix, ".plotlines.xls", sep = ""), quote = F, sep = "\t", row.names = F)

        library(ggrepel)
        # 计算每个细胞类型/聚类的其实节点
        node <- get_earliest_principal_node(cds = cds, supposed_start = unique(cds$celltype), meta_col = "celltype")
        cat("所有起始节点为：", node, "\n")
        pseudotime_file_name = "_pseudotime_auto_from_root_cluster"

        cds <- order_cells(cds, root_pr_nodes = unique(node))

        p <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 1.5)
        out_plot(plot = p, filename = paste0(outdir, "/", prefix, pseudotime_file_name))

        saveRDS(object = cds, file = paste0(outdir, "/", prefix, ".monocle3.rds"))
        # 输出轨迹分析完后的所有细胞的meta信息。
        Pseud <- data.frame("pseudotime" = cds@principal_graph_aux$UMAP$pseudotime, stringsAsFactors = F)
        Pseud$barcode <- rownames(Pseud)
        write.table(Pseud, file = paste(outdir, "/", prefix, "_metadata.txt", sep = ""), quote = F, sep = "\t", row.names = T)

        p <- p +
            theme(
                panel.background = element_rect(fill = "transparent"), # bg of the panel
                plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                panel.grid.major = element_blank(), # get rid of major grid
                panel.grid.minor = element_blank(), # get rid of minor grid
                legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                legend.box.background = element_rect(color = NA) # get rid of legend panel bg
            )

        svglite(paste0(outdir, "/pseudotime.svg"), width = 10.5, height = 7)
        print(p)
        dev.off()
        return(list(nodelist_umap = lline_m))
    } else {
        if (is.na(nodelist)) {
            stop("please input root nodes!")
        } else {
            node <- unlist(strsplit(nodelist, split = ","))
            cat("指定的轨迹根节点为：", node, "\n")
            print("use set roots to process!!!")
            pseudotime_file_name <- "_pseudotime_from_set_root_cluster"
            cds <- order_cells(cds, root_pr_nodes = node)

            Pseud <- data.frame("pseudotime" = cds@principal_graph_aux$UMAP$pseudotime, stringsAsFactors = F)
            write.table(Pseud, file = paste(outdir, "/", prefix, "_metadata.txt", sep = ""), quote = F, sep = "\t", row.names = T)

            p <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 1.5)
            out_plot(plot = p, filename = paste0(outdir, "/", prefix, pseudotime_file_name))
            p <- p +
                theme(
                    panel.background = element_rect(fill = "transparent"), # bg of the panel
                    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                    panel.grid.major = element_blank(), # get rid of major grid
                    panel.grid.minor = element_blank(), # get rid of minor grid
                    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                    legend.box.background = element_rect(color = NA) # get rid of legend panel bg
                )

            svglite(paste0(outdir, "/pseudotime_recal.svg"), width = 10.5, height = 7)
            print(p)
            dev.off()
        }
    }
}

#### MODULE 4: 差异基因分析 ####
select_partitions <- function(cds = NULL, partitions = NULL) {
    if (partitions != "all") {
        tmp_partitions <- unlist(strsplit(partitions, split = ","))
        cds <- cds[, is.element(cds@clusters$UMAP$partitions, tmp_partitions)]
    }
    return(cds)
}

replace_rBind_with_rbind <- function(func_name="calculateLW", package_name="monocle3") {
  recursive_replace <- function(expr) {
    if (is.call(expr)) {
      if (identical(expr[[1]], quote(Matrix::rBind))) {
        expr[[1]] <- as.name("rbind")
      }
      expr <- as.call(lapply(expr, recursive_replace))
    }
    return(expr)
  }

  f <- get(func_name, asNamespace(package_name))
  new_body <- recursive_replace(body(f))
  body(f) <- new_body
  assignInNamespace(func_name, f, ns = package_name)
}

#' graph_test()，可以找到不同trajectory 或不同 cluster之间的 "差异基因"。
#' @param partitions 是否选择指定的partitions进行差异分析
DE_ana <- function(cds = NULL, outdir = NULL, prefix = NULL, partitions = "all", resolution = 0.001, q_value = 0.001, nodelist = nodelist) {
    clustcol <- c("#0067AA", "#FF7F00", "#00A23F", "#FF1F1D", "#A763AC", "#B45B5D", "#FF8AB6", "#B6B800", "#01C1CC", "#85D5F8", "#FFC981", "#C8571B", "#C6CCC3", "#727272", "#EFC800", "#8A5626", "#502E91", "#59A4CE", "#344B2B", "#FBE29D", "#FDD6E6", "#849C8C", "#F07C6F", "#000101", "OrangeRed", "SlateBlue3", "DarkOrange", "GreenYellow", "Purple", "DarkSlateGray", "Gold", "DarkGreen", "DeepPink2", "Red4", "#4682B4", "#FFDAB9", "#708090", "#836FFF", "#CDC673", "#CD9B1D", "#FF6EB4", "#CDB5CD", "#008B8B", "#43CD80", "#483D8B", "#66CD00", "#CDC673", "#CDAD00", "#CD9B9B", "#FF8247", "#8B7355", "#8B3A62", "#68228B", "#CDB7B5", "#CD853F", "#6B8E23", "#E6E6FA", "#FFDAB9")

    if (nodelist != 'auto') {
        cat("Not use auto rootlist\n")
        node <- unlist(strsplit(nodelist, split = ","))
        time_file <- paste(outdir, "/", prefix, "_metadata.txt", sep = "")
        # 存在文件的不重新计算时间，理论上文件应该都存在
        if (file.exists(time_file)) {
            cat("Use pseudotimefile to DE_ana\n")
            timef <- read.table(time_file, sep = "\t", stringsAsFactors = F)
            cds@principal_graph_aux[['UMAP']]$root_pr_nodes <- node
            cds@principal_graph_aux[['UMAP']]$pseudotime <- timef$pseudotime
            names(cds@principal_graph_aux[['UMAP']]$pseudotime) <- row.names(colData(cds))
        }else {
            # 重新计算pseudotime,很费时
            cat("warning: pseudotimefile not find, to recalculate pseudotime\n")
            node <- unlist(strsplit(nodelist, split = ","))
            cds <- order_cells(cds, root_pr_nodes = node)
        }
    }
    cds <- select_partitions(cds = cds, partitions = partitions)
    save_file2 <- paste0(outdir, "/", prefix, ".partition_", gsub(',', '_', partitions), ".rds")
    saveRDS(cds, file = save_file2)
    # Track_genes colnames in sequence(7 in total): status, p_value, morans_test_statistic, morans_I, id, gene_short_name, q_value
    # @return morans_I 如果想按 effect大小 对基因进行排序，则由Morans_i column对此表进行排序。
    # 该表从-1到+1。值0表示无效果，而+1表示完美的阳性自相关，并表明附近的细胞具有基因表达的非常相似的值。大于零的重要值通常是罕见的。
    # 正值表示基因在UMAP空间的焦点区域中表达（例如，特定于一个或多个簇）。
    replace_rBind_with_rbind()
    Track_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 2)
    # 对Track_genes的列进行重排，并依据q_value值进行过滤，Track_gene自带行名，为基因名称，id列
    Track_genes <- Track_genes[, c("id", "p_value", "morans_test_statistic", "morans_I", "status", "gene_short_name")] %>%
        dplyr::filter(Track_genes$q_value <= q_value)

    genelist <- pull(Track_genes, gene_short_name) %>% as.character()
    gene_module <- find_gene_modules(cds[genelist,], resolution = resolution, cores = 8)

    cell_group <- tibble::tibble(cell = row.names(colData(cds)), cell_group = colData(cds)$celltype)
    # aggregate_gene_expression, 将基因表达矩阵换算成 基因模块表达矩阵，以便进行模块热图
    agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
    row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
    p <- pheatmap::pheatmap(agg_mat, scale = "column", clustering_method = "ward.D2", color = colorRampPalette(c("#DFDED9", "#FEFB8F", "#E58158", "#B20C1E"))(50))
    out_plot(plot = p, filename = paste0(outdir, "/", prefix, "_Genes_Module_heatmap"))

    svglite(paste0(outdir, "/heatmap_module.svg"), width = 7, height = 7)
    print(p)
    dev.off()

    rownames(gene_module) <- gene_module$id
    gene_module <- gene_module[Track_genes$id,]
    gene_module <- gene_module[, setdiff(colnames(gene_module), 'id')]

    gene_result <- cbind(Track_genes, gene_module)
    write.table(gene_result, paste(outdir, "/", prefix, "_Trajectory_genes.xls", sep = ""), sep = "\t", quote = F, row.names = F)

    # Choose top10 genes
    Track_genes_sig <- Track_genes %>%
        top_n(n = 10, morans_I) %>%
        pull(gene_short_name) %>%
        as.character()
    # 画图top 10 DE genes的伪时间表达变化抖动图
    p <- plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by = "celltype", min_expr = 0.5, ncol = 2) +
        scale_color_manual(values = clustcol)
    out_plot(plot = p, filename = paste0(outdir, "/", prefix, "_Genes_Jitterplot"))
    # 画图top 10 DE genes的在UAMP图中表达分布图，注意表达量取对数展示
    p <- plot_cells(cds, genes = Track_genes_sig, show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE) +
        scale_color_viridis_c(name = "log10(Expression)")
    p$facet$params$ncol <- round(sqrt(10))
    out_plot(plot = p, filename = paste0(outdir, "/", prefix, "_Genes_Featureplot"))
    return(list(gene_result = gene_result))
}


GeneModule_exp_UMAP <- function(cds = NULL, gene_module = NULL, outdir = NULL, prefix = NULL) {
    if (class(gene_module) == "character") gene_module <- read.csv(paste0(outdir, "/", prefix, "_Genes_Module.csv"))
    # 画出gene module在UMAP坐标中的表达情况分布图
    f1 <- function(cds = NULL, gene_module = NULL) {
        p <- plot_cells(cds = cds, genes = gene_module, label_cell_groups = FALSE, show_trajectory_graph = FALSE) +
            scale_color_viridis_c(name = "Expression Score")
        return(p)
    }

    # Draw all modules on one plot
    suboutdir <- paste(outdir, "/", "module", "/", sep = "")
    if (!dir.exists(suboutdir)) dir.create(suboutdir, recursive = TRUE)
    p <- f1(cds = cds, gene_module = gene_module)
    p$facet$params$ncol <- round(sqrt(length(levels(gene_module$module))))
    out_plot(plot = p, filename = paste0(suboutdir, "/", prefix, "_all_module_Featureplot"))
    # 画出每一个module在UAMP图上的表达分布图
    for (i in levels(gene_module$module)) {
        tmp_gene_module <- gene_module[gene_module$module == i,]
        p <- f1(cds = cds, gene_module = tmp_gene_module)
        out_plot(plot = p, filename = paste0(suboutdir, "/module_", i, "_Featureplot"))
    }
}


plotheatmap_genes_in_pseudotime <- function(cds_subset, min_expr = NULL, panel_order = NULL) {
    f_id <- NA
    Cell <- NA
    colData(cds_subset)$pseudotime <- pseudotime(cds_subset)
    cds_subset <- cds_subset[, is.finite(colData(cds_subset)$pseudotime)]
    cds_exprs <- SingleCellExperiment::counts(cds_subset)
    cds_exprs <- Matrix::t(Matrix::t(cds_exprs) / size_factors(cds_subset))
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    if (is.null(min_expr)) {
        min_expr <- 0
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_colData <- colData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_colData, by.x = "Cell", by.y = "row.names")
    cds_exprs$feature_label <- cds_exprs$f_id
    cds_exprs$f_id <- as.character(cds_exprs$f_id)
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    if (!is.null(panel_order)) {
        cds_exprs$feature_label <- factor(
            cds_exprs$feature_label,
            levels = panel_order
        )
    }

    data <- as.data.table(cds_exprs[, c('feature_label', 'pseudotime', 'expression')])
    newdata <- data.frame(Pseudotime = seq(min(data$pseudotime), max(data$pseudotime), length.out = 100)) ##100 bin
    data$pseudotime <- findInterval(data$pseudotime, newdata$Pseudotime)
    keys <- c('feature_label', 'pseudotime')
    data <- data[, list(expression = mean(expression)), keys]

    p_cor <- ggplot(as.data.frame(data), aes(pseudotime, feature_label)) +
        geom_tile(aes(fill = expression)) +
        theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        scale_fill_gradientn(colors = colorRampPalette(c("#DFDED9", "#FEFB8F", "#E58158", "#B20C1E"))(50)) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_discrete(expand = c(0, 0))
    # scale_fill_gradient2(low='blue',high='red',mid='yellow',na.value = "grey50",midpoint=midv)

    p_cor <- p_cor +
        theme(
            panel.background = element_rect(fill = "transparent"), # bg of the panel
            plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
            panel.grid.major = element_blank(), # get rid of major grid
            panel.grid.minor = element_blank(), # get rid of minor grid
            legend.background = element_rect(fill = "transparent"), # get rid of legend bg
            legend.box.background = element_rect(color = NA) # get rid of legend panel bg
        )

    svglite("pseudotime_change.svg", width = 10.5, height = 7)
    print(p_cor)
    dev.off()
}

###############################main function################################################################

#function step1, to get basic result,在outdir输出一系列结果，并返回网页展示所需的绘图dataframe
RunMonocle3 <- function(h5ad_file, annot_key, subcelltype, subsample, reduction, batch, knei, prefix = 'out', outdir) {
    load_env()
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    subcelltype <- unlist(strsplit(subcelltype, split = ','))
    subsample <- unlist(strsplit(subsample, split = ','))

    PRO <- to_rds(h5ad_file = h5ad_file, annot_key = annot_key, batch = batch)
    PRO <- GetTarget(X = PRO, subcelltype = subcelltype, subsample = subsample)

    PRO <- SetIdent(PRO, value = PRO@meta.data[, 'celltype'])
    DefaultAssay(PRO) <- "RNA"

    if (reduction)
       PRO <- re_reduction(PRO, batch=batch)

    PRO <- subset(x = PRO, downsample = 1000)

    cat("选择 ", subsample, " 样本进行后续分析\n")
    cat("选择 ", subcelltype, " 细胞类型进行后续分析\n")
    cat("选择数据子集后的细胞数为：", dim(PRO)[2], "cells\n")
    cat("具体样本信息\n")
    print(table(PRO$Sample.ID))
    cat("具体细胞类型信息\n")
    print(table(PRO$celltype))

    my_cds <- customed_createCDS(PRO = PRO)
    # return三张轨迹图及nodelist坐标表格
    result_list <- trajectory_ana(cds = my_cds, set_root = FALSE, knei = knei, outdir = outdir, prefix = prefix)
    return(result_list)
}

# 重新计算假时间
Recal <- function(monocle3file, nodelist, knei = 20, prefix = 'out') {
    load_env()
    my_cds <- readRDS(monocle3file)
    # outdir <- paste0(dirname(monocle3file),'/')
    outdir <- './'
    prefix <- paste0(prefix, '_', gsub(',', '_', nodelist))

    # return重新计算得到的轨迹图
    trajectory_ana(cds = my_cds, set_root = TRUE, nodelist = nodelist, outdir = outdir, prefix = prefix)
}

# 差异基因分析
diffana <- function(monocle3file, nodelist, partitions, resolution = 0.001, q_value = 0.05, prefix = 'out') {
    load_env()
    my_cds <- readRDS(monocle3file)
    # outdir <- paste0(dirname(monocle3file),'/')
    outdir <- './'
    prefix <- paste0(prefix, '_', gsub(',', '_', nodelist))

    # return差异基因表格及genemodule热图
    gene_result <- DE_ana(cds = my_cds, outdir = outdir, prefix = prefix, partitions = partitions, resolution = resolution, q_value = q_value, nodelist = nodelist)
    return(gene_result)
}

# 基因随假时间表达变化
genes_in_pseudotime <- function(monocle3file, genelist) {
    load_env()
    clustcol <- c("#0067AA", "#FF7F00", "#00A23F", "#FF1F1D", "#A763AC", "#B45B5D", "#FF8AB6", "#B6B800", "#01C1CC", "#85D5F8", "#FFC981", "#C8571B", "#C6CCC3", "#727272", "#EFC800", "#8A5626", "#502E91", "#59A4CE", "#344B2B", "#FBE29D", "#FDD6E6", "#849C8C", "#F07C6F", "#000101", "OrangeRed", "SlateBlue3", "DarkOrange", "GreenYellow", "Purple", "DarkSlateGray", "Gold", "DarkGreen", "DeepPink2", "Red4", "#4682B4", "#FFDAB9", "#708090", "#836FFF", "#CDC673", "#CD9B1D", "#FF6EB4", "#CDB5CD", "#008B8B", "#43CD80", "#483D8B", "#66CD00", "#CDC673", "#CDAD00", "#CD9B9B", "#FF8247", "#8B7355", "#8B3A62", "#68228B", "#CDB7B5", "#CD853F", "#6B8E23", "#E6E6FA", "#FFDAB9")
    genes <- unlist(strsplit(genelist, split = ','))
    my_cds <- readRDS(monocle3file)
    if (length(genes) == 1) {
        p <- plot_genes_in_pseudotime(my_cds[genes,], color_cells_by = "celltype", min_expr = 0.5, ncol = 2) + scale_color_manual(values = clustcol)
        p <- p +
            theme(
                panel.background = element_rect(fill = "transparent"), # bg of the panel
                plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                panel.grid.major = element_blank(), # get rid of major grid
                panel.grid.minor = element_blank(), # get rid of minor grid
                legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                strip.background = element_rect(fill = "transparent", color = NA),
                legend.box.background = element_rect(color = NA) # get rid of legend panel bg
            )
        svglite("pseudotime_change.svg", width = 10.5, height = 7)
        print(p)
        dev.off()
    }else {
        plotheatmap_genes_in_pseudotime(cds_subset = my_cds[genes,])
    }
}

