#0 umap点的大小——暂时不使用
dotSize <- function() {
  if (dim(t@meta.data)[1] > 6500) {
    pt_use <- 0.1
  } else if (dim(t@meta.data)[1] > 5500) {
    pt_use <- 0.15
  } else if (dim(t@meta.data)[1] > 4000) {
    pt_use <- 0.2
  } else if (dim(t@meta.data)[1] > 2500) {
    pt_use <- 0.3
  } else if (dim(t@meta.data)[1] > 1000) {
    pt_use <- 0.4
  }
  
  return(pt_use)
}

#pt_use <- dotSize()

#1 umap上展示top克隆型（所有样本克隆型放一起筛选top）
topCloneUmap <- function() {
  
  #1.1 定义cloneDefine，用于后面筛选数据
  if (cloneDefine == "VDJ") {
    define = "CTgene"
  } else if(cloneDefine == "CDR3"){
    define = "CTaa"
  } else if(cloneDefine == "VDJ_CDR3"){
	define = "CTstrict"	
  } else if(cloneDefine == "VDJ_CDR3aa"){
	define = "CTstrictaa"
 }
  
  #1.2 对meta.data进行筛选
  meta.data <- dplyr::select(t@meta.data, get("define")) %>%
    rownames_to_column(var = "cellid") %>%
    dplyr::rename(cellid = cellid, clonetype = get("define")) %>%
    mutate(clonetype = if_else(is.na(clonetype), "No_clonetype", clonetype))
  
  #1.3 去除get("define")中含NA的行
  data <- dplyr::filter(meta.data, clonetype != "No_clonetype")
  
  #1.4 根据get("define")统计频率
  data_fre <- as.data.frame(table(data$clonetype))
  data_fre$Var1 <- as.character(data_fre$Var1)
  data_fre$Freq <- as.numeric(data_fre$Freq)
  
  #1.5 获取top克隆型
  topClone <- arrange(data_fre, desc(Freq)) %>%
    pull(Var1) %>%
    head(top)
    
  #1.6 top克隆命名
  num <- seq(top)
  clone <- rep("Topclone", top)
  cloneName <- str_c(clone, num)
  
  #1.7 克隆差集
  intersectClone <- dplyr::setdiff(data_fre$Var1, topClone)
  
  #1.8 将差集克隆命名为others
  for (i in intersectClone) {
    meta.data$clonetype[meta.data$clonetype == i] <- "Others_clonetype"
  }
  
  #1.9 重命名top克隆
  col_clone <- meta.data$clonetype
  for (i in 1:length(topClone)) {
    col_clone <- str_replace_all(col_clone, topClone[i], cloneName[i])
  }
  meta.data$clonetype <- col_clone
  
  #1.10 获取umap坐标
  umap_data <- as.data.frame(Embeddings(object = t, reduction = "umap")[,1:2])
  
  #1.11 将meta.data3与umap_data合并
  data <- full_join(rownames_to_column(umap_data, var = "cellid"), meta.data, by = "cellid")
  
  #1.12 颜色
  color1 <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B",
              "#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101",
              "#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B",
              "#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B",
              "#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060",
              "#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC",
              "#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8")
  color2 <- c("#CCFFFF", "#CCCCCC")
  
  color_ce <- c(color1[1:top], color2)
  names(color_ce) <- c(cloneName, "Others_clonetype", "No_clonetype")
  
  #1.13 画图
	########################3
  p <- ggplot(data) +
    geom_point(aes(UMAP_1, UMAP_2, color = clonetype)) +
    geom_point(data = dplyr::filter(data, clonetype == "Others_clonetype"), aes(UMAP_1, UMAP_2, color = clonetype)) +
    geom_point(data = dplyr::filter(data, clonetype %in% cloneName), aes(UMAP_1, UMAP_2, color = clonetype)) +
    scale_color_discrete(limits = c(cloneName, "Others_clonetype", "No_clonetype")) +
    scale_color_manual(values = color_ce) +
    labs(color = "") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme_half_open() +
    theme(legend.position = "right", aspect.ratio = 1)
  #输出表格
  write.table(data.frame(ClonoType = topClone,topClone = cloneName),paste0(outdir,"/07.sc.UmapCloneType/",prefix,".top_clonotype.xls"),row.names = F,quote = F,sep = "\t")
  pdf(str_c(outdir, "/07.sc.UmapCloneType/", prefix, ".top_umap.pdf"), width = 8,height = 9)
  print(p)
  dev.off()
  
  png(str_c(outdir, "/07.sc.UmapCloneType/", prefix,".top_umap.png"), width = 800,height = 900)
  print(p)
  dev.off()
}

#2 umap上展示含克隆型的细胞
CloneUmap <- function() {
  
  #1.1 定义cloneDefine，用于后面筛选数据
  if (cloneDefine == "VDJ") {
    define = "CTgene"
  } else {
    define = "CTaa"
  }
  
  #1.2 对meta.data进行筛选
  meta.data <- dplyr::select(t@meta.data, get("define")) %>%
    rownames_to_column(var = "cellid") %>%
    dplyr::rename(cellid = cellid, clonetype = get("define")) %>%
    mutate(clonetype = if_else(is.na(clonetype), str_c("No", get("vdjType"), sep = "_"), get("vdjType")))
  
  #1.3 获取umap坐标
  umap_data <- as.data.frame(Embeddings(object = t, reduction = "umap")[,1:2])
  
  #1.11 将meta.data3与umap_data合并
  data <- full_join(rownames_to_column(umap_data, var = "cellid"), meta.data, by = "cellid")
  
  #1.12 颜色
#  color_ce <- c("#CCFFFF", "#CCCCCC")
  color_ce <- c("#063364", "#CCCCCC")
  names(color_ce) <- c(get("vdjType"), str_c("No", get("vdjType"), sep = "_"))
  
  #1.13 画图
  p <- ggplot(data) +
#    geom_point(aes(UMAP_1, UMAP_2, color = clonetype),size = 2) +
#	scale_color_manual(values = color_re[2])+
	geom_point(aes(UMAP_1, UMAP_2, color = clonetype),size = 0.75)+
#    geom_point(data = dplyr::filter(data, clonetype == get("vdjType")), aes(UMAP_1, UMAP_2, fill = clonetype),size = 2,shape = 21,stroke = 0.5) +
	geom_point(data = dplyr::filter(data, clonetype == get("vdjType")), aes(UMAP_1, UMAP_2, color = clonetype),size = 0.75) +
#    scale_fill_discrete(limits = c(get("vdjType"), str_c("No", get("vdjType"), sep = "_"))) +
 #   scale_fill_manual(values = color_ce[1]) +
	scale_color_discrete(limits = c(get("vdjType"), str_c("No", get("vdjType"), sep = "_"))) +
	scale_color_manual(values = color_ce) +
    labs(color = "") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme_half_open() +
    theme(legend.position = "right", aspect.ratio = 1)
  
  pdf(str_c(outdir, "/07.sc.UmapCloneType/", prefix, ".umap.pdf"), width = 8,height = 9)
  print(p)
  dev.off()
  
  png(str_c(outdir, "/07.sc.UmapCloneType/", prefix,".umap.png"), width = 800,height = 900)
  print(p)
  dev.off()
}










topCloneUmap_Sample <- function() {
  
  #1.1 定义cloneDefine，用于后面筛选数据
  if (cloneDefine == "VDJ") {
    define = "CTgene"
  } else if(cloneDefine == "CDR3"){
    define = "CTaa"
  } else if(cloneDefine == "VDJ_CDR3"){
	define = "CTstrict"	
  } else if(cloneDefine == "VDJ_CDR3aa"){
	define = "CTstrictaa"
 }
  #1.2 对meta.data进行筛选
  t@meta.data[!is.na(t@meta.data$CTgene),define] <- paste0(t@meta.data[!is.na(t@meta.data$CTgene),define],":",t@meta.data[!is.na(t@meta.data$CTgene),"sample"])
  meta.data <- dplyr::select(t@meta.data, get("define")) %>%
    rownames_to_column(var = "cellid") %>%
    dplyr::rename(cellid = cellid, clonetype = get("define")) %>%
    mutate(clonetype = if_else(is.na(clonetype), "No_clonetype", clonetype))
  
  #1.3 去除get("define")中含NA的行
  data <- dplyr::filter(meta.data, clonetype != "No_clonetype")
  #1.4 根据get("define")统计频率
  data_fre <- as.data.frame(table(data$clonetype))
  data_fre$Var1 <- as.character(data_fre$Var1)
  data_fre$Freq <- as.numeric(data_fre$Freq)
  
  #1.5 获取top克隆型
  topClone <- arrange(data_fre, desc(Freq)) %>%
    pull(Var1) %>%
    head(top)
  
  #1.6 top克隆命名
  num <- seq(top)
  clone <- rep("Topclone", top)
  cloneName <- str_c(clone, num)
  
  #1.7 克隆差集
  intersectClone <- dplyr::setdiff(data_fre$Var1, topClone)
  #1.8 将差集克隆命名为others
  for (i in intersectClone) {
    meta.data$clonetype[meta.data$clonetype == i] <- "Others_clonetype"
  }
  
  #1.9 重命名top克隆
  col_clone <- meta.data$clonetype
  for (i in 1:length(topClone)) {
    col_clone <- str_replace_all(col_clone, topClone[i], cloneName[i])
  }
  meta.data$clonetype <- col_clone
  #1.10 获取umap坐标
  umap_data <- as.data.frame(Embeddings(object = t, reduction = "umap")[,1:2])
  
  #1.11 将meta.data3与umap_data合并
  data <- full_join(rownames_to_column(umap_data, var = "cellid"), meta.data, by = "cellid")
  
  #1.12 颜色
  color1 <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B",
              "#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101",
              "#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B",
              "#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B",
              "#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060",
              "#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC",
              "#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8")
  color2 <- c("#CCFFFF", "#CCCCCC")
  
  color_ce <- c(color1[1:top], color2)
  names(color_ce) <- c(cloneName, "Others_clonetype", "No_clonetype")
  
  #1.13 画图
  p <- ggplot(data) +
    geom_point(aes(UMAP_1, UMAP_2, color = clonetype)) +
    geom_point(data = dplyr::filter(data, clonetype == "Others_clonetype"), aes(UMAP_1, UMAP_2, color = clonetype)) +
    geom_point(data = dplyr::filter(data, clonetype %in% cloneName), aes(UMAP_1, UMAP_2, color = clonetype)) +
    scale_color_discrete(limits = c(cloneName, "Others_clonetype", "No_clonetype")) +
    scale_color_manual(values = color_ce) +
    labs(color = "") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme_half_open() +
    theme(legend.position = "right", aspect.ratio = 1)
  #输出表格
  write.table(data.frame(ClonoType = gsub("[:].*$","",topClone),topClone = cloneName,Sample = gsub("^.*[:]","",topClone)),paste0(outdir,"/07.sc.UmapCloneType/",prefix,".top_clonotype_sample.xls"),row.names = F,quote = F,sep = "\t")
  pdf(str_c(outdir, "/07.sc.UmapCloneType/", prefix, ".top_umap_sample.pdf"), width = 8,height = 9)
  print(p)
  dev.off()
  
  png(str_c(outdir, "/07.sc.UmapCloneType/", prefix,".top_umap_sample.png"), width = 800,height = 900)
  print(p)
  dev.off()
}

