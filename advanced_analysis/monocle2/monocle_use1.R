##3.6 绘制density-流图
###3.6.1 按样本分组统计
sampleRatio <- function(cs) {
  sample_r <- group_by(cs, Pseudotime, sample, cluster) %>%
    summarise(count = n())
}

###3.6.2 按样本分组分组统计
groupRatio <- function(cg) {
  group_r <- group_by(cg, Pseudotime, group, cluster) %>%
    summarise(count = n())
}

###3.6.3 绘图函数
ggStream <- function(pictureData) {
  ggplot(pictureData, aes(x = Pseudotime, y = count, fill = cluster)) +
    geom_stream(extra_span = 0.1, true_range = "both") +
    scale_x_continuous(limits = c(0, 10), expand = c(0, 0), breaks = seq(0, 10, 2.5)) + 
    scale_fill_manual(values = color_use, name = "Cluster") +
    coord_flip() +
    #facet_wrap(vars(as.name(dplyr::intersect(colnames(pictureData), c("sample", "group")))), nrow = 1) + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.background = element_rect(fill = NA),
          axis.title = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = NA),
          #text = element_text(size = 15),
          legend.spacing.x = unit(0.25, "cm"),
          plot.title = element_text(face = "bold", size = 20),
          plot.margin = unit(c(1,1,1,1),"cm"))
}

###3.6.4 main
violinDensity <- function(data) {
  a <- dplyr::select(data, sample, cluster, Pseudotime)
  s <- dplyr::mutate(a, sample = "all_sample")
  cs <- rbind(a, s)
  
  #设置样本顺序
  if (class(a$sample) == "factor") {
    cs$sample <- factor(cs$sample, levels = c(levels(a$sample), "all_sample"))
  } else {
    cs$sample <- factor(cs$sample, levels = c(sort(unique(a$sample)), "all_sample"))
  }
  
  #分组统计
  sample_r <- sampleRatio(cs)
  
  #绘图
  ps <- ggStream(sample_r) + facet_wrap(~sample, nrow = 1)
  
  if (any(str_detect("^group$", colnames(data)))) {
    #筛选数据
    a <- dplyr::select(data, group, cluster, Pseudotime)
    g <- dplyr::mutate(a, group = "all_group")
    cg <- rbind(a, g)
    
    #设置样本分组顺序
    if (class(cg$group) == "factor") {
      cg$group <- factor(cg$group, levels = c(levels(a$group), "all_group"))
    } else {
      cg$group <- factor(cg$group, levels = c(sort(unique(a$group)), "all_group"))
    }
    
    #分组统计
    group_r <- groupRatio(cg)
    
    #绘图
    pg <- ggStream(group_r) + facet_wrap(~group, nrow = 1)
    
  } 
  
  if (any(str_detect("^group$", colnames(data)))) {
    density <- vector("list", 2)
    density[[1]]<-ps
    density[[2]]<-pg
    return(density)
  } else {
    density <- vector("list", 1)
    density[[1]]<-ps
    return(density)
  }
}

##3.7 轨迹降维图

###3.7.1 非分面轨迹降维图
reductionTinotfacet <- function(my_cds) {
  
  #cluster
  p1 <- plot_cell_trajectory(my_cds, color_by = "cluster", cell_size = cell_size, show_tree = show_tree_argv, show_branch_points = show_branch_points_argv) + 
    scale_color_manual(values = color_use, name = "Cluster") + 
    theme(legend.position='right', aspect.ratio=1,
          axis.line=element_line(color="black",size=2),
          axis.text=element_text(size=12,face = "bold",color="black")) + 
    guides(color = guide_legend(override.aes = list(size = 4), ncol = ncol))
  
  #Pseudotime
  p2 <- plot_cell_trajectory(my_cds, color_by = "Pseudotime", cell_size = cell_size, show_tree = show_tree_argv, show_branch_points = show_branch_points_argv) +
    scale_color_viridis_c() + 
    theme(legend.position='right', aspect.ratio=1,
          axis.line=element_line(color="black",size=2),
          axis.text=element_text(size=12,face = "bold",color="black"))
  
  #state
  p3 <- plot_cell_trajectory(my_cds, color_by = "State",cell_size = cell_size, show_tree = show_tree_argv, show_branch_points = show_branch_points_argv)  + 
    theme(legend.position='right', aspect.ratio=1,
          axis.line=element_line(color="black",size=2),
          axis.text=element_text(size=12,face = "bold",color="black")) + 
    guides(color = guide_legend(override.aes = list(size = 4), ncol=ncol))
  #group
  if (any(str_detect("^group$", colnames(pData(my_cds))))) {
  p5 <- plot_cell_trajectory(my_cds, color_by = "group",cell_size = cell_size, show_tree = show_tree_argv, show_branch_points = show_branch_points_argv) +
    scale_color_manual(values = groupcol_use) +
    theme(legend.position='right', aspect.ratio=1,
          axis.line=element_line(color="black",size=2),
          axis.text=element_text(size=12,face = "bold",color="black")) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol=ncol))}
  #sample
  p4 <- plot_cell_trajectory(my_cds, color_by = "sample",cell_size = cell_size, show_tree = show_tree_argv, show_branch_points = show_branch_points_argv) +
    scale_color_manual(values = samplecol_use) +
    theme(legend.position='right', aspect.ratio=1,
          axis.line=element_line(color="black",size=2),
          axis.text=element_text(size=12,face = "bold",color="black")) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol=ncol))
  if (any(str_detect("^group$", colnames(pData(my_cds))))) {
  reductionTinotfacet <- vector("list", 5)
  reductionTinotfacet[[1]] <- p1
  reductionTinotfacet[[2]] <- p2
  reductionTinotfacet[[3]] <- p3
  reductionTinotfacet[[4]] <- p4
  reductionTinotfacet[[5]] <- p5
  return(reductionTinotfacet)
 }else{
  reductionTinotfacet <- vector("list", 4)
  reductionTinotfacet[[1]] <- p1
  reductionTinotfacet[[2]] <- p2
  reductionTinotfacet[[3]] <- p3
  reductionTinotfacet[[4]] <- p4
  return(reductionTinotfacet)
 }
}
###3.7.2 分面轨迹降维图
reductionTifacet <- function(my_cds) {
  #facet_cluster
  p1 <- plot_cell_trajectory(my_cds, color_by = "cluster",cell_size = cell_size, show_tree = show_tree_argv, show_branch_points = show_branch_points_argv) + 
    scale_color_manual(values = color_use, name = "Cluster") +
    facet_wrap(~cluster) + 
    theme(legend.position='none', aspect.ratio=1,
          axis.line=element_line(color="black",size=2),
          axis.text=element_text(size=12,face = "bold",color="black"))
  
  #facet_sample
  p2 <- plot_cell_trajectory(my_cds, color_by = "cluster", cell_size = cell_size, show_tree = show_tree_argv, show_branch_points = show_branch_points_argv) + 
    scale_color_manual(values = color_use, name = "Cluster") +
    facet_wrap(~sample) + 
    theme(legend.position='right', aspect.ratio=1,
          axis.line=element_line(color="black",size=2),
          axis.text=element_text(size=12,face = "bold",color="black")) + 
    guides(color = guide_legend(override.aes = list(size = 4), ncol=ncol))
  
  #facet_State
  p3 <- plot_cell_trajectory(my_cds, color_by = "State") +
    facet_wrap(~State, nrow = 2) + 
    theme(axis.line=element_line(color="black",size=2),
          axis.text=element_text(size=12,face = "bold",color="black"))
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  #facet_group
  if (any(str_detect("^group$", colnames(pData(my_cds))))) {
    p4 <- plot_cell_trajectory(my_cds, color_by = "cluster", cell_size = cell_size, show_tree = show_tree_argv, show_branch_points = show_branch_points_argv) + 
      scale_color_manual(values = color_use, name = "Cluster") +
      facet_wrap(~group) + 
      theme(legend.position='right', aspect.ratio=1,
            axis.line=element_line(color="black",size=2),
            axis.text=element_text(size=12,face = "bold",color="black")) + 
      guides(color = guide_legend(override.aes = list(size = 4), ncol = ncol))
  }
  
  if (any(str_detect("^group$", colnames(pData(my_cds))))) {
    reductionTifacet <- vector("list", 4)
    reductionTifacet[[1]] <- p1
    reductionTifacet[[2]] <- p2
    reductionTifacet[[3]] <- p3
    reductionTifacet[[4]] <- p4
    return(reductionTifacet)
  } else {
    reductionTifacet <- vector("list", 3)
    reductionTifacet[[1]] <- p1
    reductionTifacet[[2]] <- p2
    reductionTifacet[[3]] <- p3
    return(reductionTifacet)
  }
}

###3.7.2 complex树状图
complexfacet <- function(my_cds) {
   p1 <- plot_complex_cell_trajectory(my_cds,color="cluster",cell_size = cell_size, show_tree = show_tree_argv, show_branch_points = show_branch_points_argv)+
          scale_color_manual(values = color_use, name = "Cluster") +         
          theme(legend.position='right', aspect.ratio=1,
          axis.line=element_line(color="black",size=1.5),
          axis.text=element_text(size=12,face = "bold",color="black")) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol=ncol))

   p2 <- plot_complex_cell_trajectory(my_cds,color="cluster",cell_size = cell_size, show_tree = show_tree_argv, show_branch_points = show_branch_points_argv)+
          scale_color_manual(values = color_use, name = "Cluster") +
          facet_wrap(~sample) +
          theme(legend.position='right', aspect.ratio=1,axis.line=element_line(color="black",size=1.5),
          axis.text=element_text(size=12,face = "bold",color="black")) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol=ncol))
   
   p3 <- plot_complex_cell_trajectory(my_cds,color="cluster",cell_size = cell_size, show_tree = show_tree_argv, show_branch_points = show_branch_points_argv)+
          scale_color_manual(values = color_use, name = "Cluster") +
          facet_wrap(~cluster) +
          theme(legend.position='right', aspect.ratio=1,axis.line=element_line(color="black",size=1.5),
          axis.text=element_text(size=12,face = "bold",color="black"))
   if (any(str_detect("^group$", colnames(pData(my_cds))))){ 
   p4 <- plot_complex_cell_trajectory(my_cds,color="cluster",cell_size = cell_size, show_tree = show_tree_argv, show_branch_points = show_branch_points_argv)+
          scale_color_manual(values = color_use, name = "Cluster") +
          facet_wrap(~group) +
          theme(legend.position='right', aspect.ratio=1,axis.line=element_line(color="black",size=1.5),
          axis.text=element_text(size=12,face = "bold",color="black")) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol=ncol))}
    if (any(str_detect("^group$", colnames(pData(my_cds))))) {
    complexfacet <- vector("list", 4)
    complexfacet[[1]] <- p1
    complexfacet[[2]] <- p2
    complexfacet[[3]] <- p3
    complexfacet[[4]] <- p4
    return(complexfacet)
  } else {
    complexfacet <- vector("list", 3)
    complexfacet[[1]] <- p1
    complexfacet[[2]] <- p2
    complexfacet[[3]] <- p3
    return(complexfacet)
  }
}

##3.8 差异分析top基因获取——gene_to_pseudotime（top8）：用于基因随细胞状态等的表达变化/monocle降维体系；gene_to_cluster（top30）：用于热图
topGene <- function(my_pseudotime_de, toPseudotime, toHeatmap) {
  if (toPseudotime == "T") {
    gene_to_pseudotime <- my_pseudotime_de %>%
      arrange(qval) %>%
      head(8) %>%
      pull(gene_short_name) %>%
      as.character() %>%
      sort()
    return(gene_to_pseudotime)
  }
  
  if (toHeatmap == "T") {
    gene_to_cluster <- my_pseudotime_de %>%
      arrange(qval) %>%
      head(30) %>%
      pull(gene_short_name) %>%
      as.character()
    
    return(gene_to_cluster)
  }
}

##3.9 top8 gene图片展示
pictureTop8.0 <- function(data) {
  
  #top8 pseudotime
  p1 <- plot_genes_in_pseudotime(data, ncol = 2, cell_size = cell_size, color_by = "cluster") +
    theme(legend.position = "top") +
    scale_color_manual(values = color_use, name = "Cluster") +
    theme(legend.position = 'right',
          axis.line=element_line(color="black",size=2),
          axis.text=element_text(size=12,face = "bold",color="black")) + 
    guides(color = guide_legend(override.aes = list(size = 4), ncol = ncol))
  
  #top8 trajectory
  ##去掉基因名称中的中文字符
  gene_to_pseudotime_use <- c()
  gene_to_pseudotime_tmp <- gene_to_pseudotime
  for (i in 1:length(gene_to_pseudotime)) {
    if (str_detect(gene_to_pseudotime[i], "-")) {
      print(str_c("gene of ", gene_to_pseudotime[i], " will be change as ", str_replace_all(gene_to_pseudotime[i], "-", "_")))
      gene_to_pseudotime[i] <- str_replace_all(gene_to_pseudotime[i], "-", "_")
    }
    gene_to_pseudotime_use <- c(gene_to_pseudotime_use, gene_to_pseudotime[i])
  }
  gene_to_pseudotime <- gene_to_pseudotime_tmp
  ##将top基因表达信息取log后加入pData(my_cds)——添加到pData(my_cds)中的基因为，去掉中文字符的基因名称，后面绘图时，使用没修改的基因名称对legend重新命名
  for (i in 1:length(gene_to_pseudotime)) {
    pData(my_cds)[[gene_to_pseudotime_use[i]]] = log10(exprs(my_cds)[gene_to_pseudotime[i],] + 0.1)
  }
  
  ##for循环绘图，并将图片保存如list
  tmp <- vector("list", length(gene_to_pseudotime_use))
  for (i in 1:length(gene_to_pseudotime_use)) {
    assign(paste("p", gene_to_pseudotime_use[i], sep = "_"), (plot_cell_trajectory(my_cds, color_by = gene_to_pseudotime_use[i]) + scale_color_gradient(low ="lightgrey", high = "red") + theme(legend.position='right', aspect.ratio=1) + labs(color = gene_to_pseudotime[i], "_", "-")))
    tmp[[i]] <- get(paste("p", gene_to_pseudotime_use[i], sep = "_"))
  }
  
  ##拼图
  p2 <- (tmp[[1]] / tmp[[3]] / tmp[[5]] / tmp[[7]]) | (tmp[[2]] / tmp[[4]] / tmp[[6]] / tmp[[8]])
  
  pictureTop8 <- vector("list", 3)
  pictureTop8[[1]] <- p1
  pictureTop8[[2]] <- p2
  pictureTop8[[3]] <- tmp
  
  return(pictureTop8)
}

pictureTop8 <- function(data) {
  
  #top8 pseudotime
  p1 <- plot_genes_in_pseudotime(data, ncol = 2, cell_size = cell_size, color_by = "cluster") +
    theme(legend.position = "top") +
    theme( strip.text.x = element_text( size = 12,face = "bold")) +
    scale_color_manual(values = color_use, name = "Cluster") +
    theme(legend.position = 'right',
          axis.line=element_line(color="black",size=2),
          axis.text=element_text(size=12,face = "bold",color="black")) + 
    guides(color = guide_legend(override.aes = list(size = 4), ncol = ncol))
  
  ##for循环绘图，并将图片保存如list
  tmp <- vector("list", length(gene_to_pseudotime))
  for (i in 1:length(gene_to_pseudotime)) {
    assign(paste("p", gene_to_pseudotime[i], sep = "_"), 
           (plot_cell_trajectory(my_cds, markers = gene_to_pseudotime[i], use_color_gradient = T, show_tree = T, show_branch_points = F) + 
              scale_color_gradient(low ="lightgrey", high = "red") + 
              theme(legend.position='right', aspect.ratio=1)))
    tmp[[i]] <- get(paste("p", gene_to_pseudotime[i], sep = "_"))
  }
  
  ##拼图
  p2 <- (tmp[[1]] / tmp[[3]] / tmp[[5]] / tmp[[7]]) | (tmp[[2]] / tmp[[4]] / tmp[[6]] / tmp[[8]])
  
  pictureTop8 <- vector("list", 3)
  pictureTop8[[1]] <- p1
  pictureTop8[[2]] <- p2
  pictureTop8[[3]] <- tmp
  
  return(pictureTop8)
}

##3.10 top 30热图
pictureTop30 <- function(data) {
  p <- plot_pseudotime_heatmap_sgr(data, num_clusters = length(unique(pData(my_cds)$State)), cores = 1, show_rownames = TRUE, return_heatmap = TRUE)
  return(p)
}

##热图修改字体颜色函数

table.ramp <- function (n, mid = 0.5, sill = 0.5, base = 1, height = 1)
{
    x <- seq(0, 1, length.out = n)
    y <- rep(0, length(x))
    sill.min <- max(c(1, round((n - 1) * (mid - sill/2)) + 1))
    sill.max <- min(c(n, round((n - 1) * (mid + sill/2)) + 1))
    y[sill.min:sill.max] <- 1
    base.min <- round((n - 1) * (mid - base/2)) + 1
    base.max <- round((n - 1) * (mid + base/2)) + 1
    xi <- base.min:sill.min
    yi <- seq(0, 1, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    xi <- sill.max:base.max
    yi <- seq(1, 0, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    height * y
}
rgb.tables <- function (n, red = c(0.75, 0.25, 1), green = c(0.5, 0.25, 1),
    blue = c(0.25, 0.25, 1))
{
    rr <- do.call("table.ramp", as.list(c(n, red)))
    gr <- do.call("table.ramp", as.list(c(n, green)))
    br <- do.call("table.ramp", as.list(c(n, blue)))
    rgb(rr, gr, br)
}
blue2green2red <- function (n){
    rgb.tables(n, red = c(0.8, 0.2, 1), green = c(0.5, 0.4, 0.8),
    blue = c(0.2, 0.2, 1))
}
plot_pseudotime_heatmap_sgr <- function (cds_subset, cluster_rows = TRUE, hclust_method = "ward.D2",
    num_clusters = 6, hmcols = NULL, add_annotation_row = NULL,
    add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE,
    norm_method = c("log", "vstExprs"), scale_max = 3, scale_min = -3,
    trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE,
    cores = 1)
{
    num_clusters <- min(num_clusters, nrow(cds_subset))
    pseudocount <- 1
    newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime),
        max(pData(cds_subset)$Pseudotime), length.out = 100))
    m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula,
        relative_expr = T, new_data = newdata)
    m = m[!apply(m, 1, sum) == 0, ]
    norm_method <- match.arg(norm_method)
    if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) ==
        FALSE) {
        m = vstExprs(cds_subset, expr_matrix = m)
    }
    else if (norm_method == "log") {
        m = log10(m + pseudocount)
    }
    m = m[!apply(m, 1, sd) == 0, ]
    m = Matrix::t(scale(Matrix::t(m), center = TRUE))
    m = m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] = 0
    m[m > scale_max] = scale_max
    m[m < scale_min] = scale_min
    heatmap_matrix <- m
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1
    if (is.null(hmcols)) {
        bks <- seq(-3.1, 3.1, by = 0.1)
        hmcols <- blue2green2red(length(bks) - 1)
    }
    else {
        bks <- seq(-3.1, 3.1, length.out = length(hmcols))
    }
    ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE,
        cluster_rows = cluster_rows, show_rownames = F, show_colnames = F,
        clustering_distance_rows = row_dist, clustering_method = hclust_method,
        cutree_rows = num_clusters, silent = TRUE, filename = NA,
        breaks = bks, border_color = NA, color = hmcols)
    if (cluster_rows) {
        annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row,
            num_clusters)))
    }
    else {
        annotation_row <- NULL
    }
    if (!is.null(add_annotation_row)) {
        old_colnames_length <- ncol(annotation_row)
        annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row),
            ])
        colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
    }
    if (!is.null(add_annotation_col)) {
        if (nrow(add_annotation_col) != 100) {
            stop("add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!")
        }
        annotation_col <- add_annotation_col
    }
    else {
        annotation_col <- NA
    }
    if (use_gene_short_name == TRUE) {
        if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
            feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix),
                "gene_short_name"])
            feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
            row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row),
                "gene_short_name"])
            row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
        }
        else {
            feature_label <- row.names(heatmap_matrix)
            row_ann_labels <- row.names(annotation_row)
        }
    }
    else {
        feature_label <- row.names(heatmap_matrix)
        if (!is.null(annotation_row))
            row_ann_labels <- row.names(annotation_row)
    }
    row.names(heatmap_matrix) <- feature_label
    if (!is.null(annotation_row))
        row.names(annotation_row) <- row_ann_labels
    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
    ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE,
        cluster_rows = cluster_rows, show_rownames = show_rownames,
        show_colnames = F, clustering_distance_rows = row_dist,
        clustering_method = hclust_method, cutree_rows = num_clusters,
        annotation_row = annotation_row, annotation_col = annotation_col,
        treeheight_row = 20, breaks = bks, fontsize = 10, color = hmcols,
        border_color = NA, silent = TRUE, filename = NA)
    grid::grid.rect(gp = grid::gpar("fill", col = NA))
    grid::grid.draw(ph_res$gtable)
    if (return_heatmap) {
        return(ph_res)
    }
}
##3.11 绘制state热图
substateheatmap <- function(cds,state,outdir,prefix,kgroups=3){
cds_sub <- cds[,pData(cds)$State %in% state]
my_pseudotime_de <- differentialGeneTest(cds_sub,fullModelFormulaStr = "~sm.ns(Pseudotime)")
#my_pseudotime_de <- my_pseudotime_de[my_pseudotime_de$pval<0.05,]
my_pseudotime_de <- subset(my_pseudotime_de, qval < 1e-4)
my_pseudotime_de %>% arrange(qval) %>% select(gene_short_name) -> gene_to_cluster
gene_to_cluster <- as.character(gene_to_cluster$gene_short_name)
print(length(gene_to_cluster))
my_pseudotime_cluster <- plot_pseudotime_heatmap_sgr(my_cds[gene_to_cluster,], num_clusters = kgroups,cores = 3,show_rownames = F,return_heatmap = TRUE)
clusterGene <- cutree(my_pseudotime_cluster$tree_row, k=kgroups)
clusterGene <- data.frame(clusterGene)
clusterGene$gene <- rownames(clusterGene)

cluster <-c()
gene <- c()
for(c in unique(clusterGene[,1])){
    clusterGene_sub <- clusterGene[clusterGene[,1]==c, 2]
    cluster <- c(cluster, paste0('cluster', c))
    gene <- c(gene, paste(clusterGene_sub, collapse=','))
}
df <- data.frame(cluster, gene)
write.table(df, paste(outdir,'/',prefix,'_state_heatmap.xls',sep=''), row.names=F,col.names=F,quote=F, sep='\t')
}
 
