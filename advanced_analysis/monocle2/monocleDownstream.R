#!/usr/bin/env Rscript
suppressMessages({
  library(argparser)
  library(Seurat)
  library(tidyverse)
  library(pheatmap)
  library(grid)
  library(cowplot)
  library(monocle)
  library(patchwork)
  library(ggsci)
  library(RhpcBLASctl)
  library(ggstream)
  library(reshape)
})

blas_set_num_threads(2)
source("./seurat_use.R")
source("./seurat_function.R")
source("./monocle_use1.R")

#1 创建参数解释器
argv <- arg_parser('')
argv <- add_argument(argv,"--inputdir", help = "the outdir of result of monocle2")
argv <- add_argument(argv,"--clusterlevels", help = "the cluster levels", default = 'F')
argv <- add_argument(argv,"--newcol", help="useing new color protocol", default = 'T')
argv <- add_argument(argv,"--colorset", help="useing  color protocol of metadata", default = 'T')
argv <- add_argument(argv,"--do", help="Three parameters for selection: setroot/ reverse/ beam/ ucell/state'setroot': set root, 'reverse' reverse Pseudotime: , 'beam': heatmap of brach,'ucell':ucell mapping Pseudotime,'state':heatmap of state",default = 'F')
argv <- add_argument(argv,"--rootcluster", help="the root cluster", default="F")
argv <- add_argument(argv,"--rootstate", help="the root state", default="F")
argv <- add_argument(argv,"--branch", help="heatmap of branch", default="F")
argv <- add_argument(argv,"--branch_cluster", help="heatmap of branch", default="default")
argv <- add_argument(argv,"--state", help="heatmap of state,split by ','", default="F")
argv <- add_argument(argv,"--genesetscore", help="the different geneset score result")
argv <- add_argument(argv,"--geneset", help="the different geneset,split by ','")
argv <- add_argument(argv,"--show_tree", help="Trajectories show tree", default = TRUE)
argv <- add_argument(argv,"--show_branch_points", help="Trajectories show tree", default = FALSE)
argv <- add_argument(argv,"--prefix", help="output prefix")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- add_argument(argv,"--core", help="number of threads core",default = 2)
argv <- parse_args(argv)

#2解析参数
inputdir <- argv$inputdir

if (str_detect(argv$clusterlevels, ",")) {
  clusterlevels <- unlist(strsplit(argv$clusterlevels, split = ","))
} else {
  clusterlevels <- argv$clusterlevels
}

newcol <- argv$newcol
colorset <- argv$colorset
do <- argv$do
rootcluster <- argv$rootcluster
rootstate <- argv$rootstate
branch <- as.numeric(argv$branch)
core <- as.numeric(argv$core)
#branch <- argv$branch
branch_cluster <- argv$branch_cluster
state <- argv$state
show_tree_argv <- as.logical(argv$show_tree)
show_branch_points_argv <- as.logical(argv$show_branch_points)
prefix <- argv$prefix
outdir <- argv$outdir
dir.create(outdir)

feedback <- str_c(outdir, "/", prefix)   #反馈目录
dir.create(feedback)

#3 加载数据
file_name <- list.files(inputdir)

if (str_detect(inputdir, "/$")) {
  inputdir <- str_sub(inputdir, 1, str_length(inputdir)-1)
} else {
  inputdir <- inputdir
}

pro <- readRDS(str_c(inputdir, "/", file_name[str_detect(file_name, ".cluster_as_state.rds$")]))
my_cds <- readRDS(str_c(inputdir, "/", file_name[str_detect(file_name, ".pseudotime.cds$")]))

pro@meta.data$cluster <- pro@meta.data$cluster_old
  Idents(pro) <- "cluster"

#3 分析

##3.1 根据monocle输出结果设置root或reverse Pseudotime——因为在monocle2中输入及参数完全一致，输出结果可能不一样，所以增加此脚本
if (do %in% c("setroot", "reverse")) {
  
 ###3.1.1修改cluster因子水平
  if (!(identical(clusterlevels, "F"))) {
    pro <- clusterLevels(pro, clusterlevels)
  }
  
  ###3.1.2 设置色板
 # clustcol <- colorSet(newcol)
 # color_use <- clustcol[1:length(levels(pro))]
 # names(color_use) <- levels(pro)
pro$cluster <- as.character(pro$cluster)
pro$sample <- as.character(pro$sample)
colorlist <- getColors(pro,meta=c("sample","group","cluster","raw_cluster"))
if (colorset != "F"){
   if (length(colorlist)!=0){
       clustcol <- colorlist$cluster_colors
       #clustcol <- as.data.frame(clustcol)
       #clustcol <- rownames_to_column(clustcol,'cluster')
       #clustcol <- clustcol %>%
       #mutate(cluster=fct_relevel(cluster,levels(pro))) %>%
       #arrange(cluster)
       #clustcol <- column_to_rownames(clustcol,'cluster')
       #a <- as.character(clustcol$clustcol)
       #names(a) <- rownames(clustcol)
       #clustcol <- a
  }else{
       clustcol <- colorSet(newcol)
  }
}else{
clustcol <- colorSet(newcol)
}
color_use <- clustcol[1:length(unique(pro@meta.data$cluster))]
names(color_use) <- levels(pro)

# my_cds <- orderCells(my_cds)
  
  ###3.1.3 要干啥
  #### 基于之前的cds设置起点
  if (do == "setroot") {
    if (rootcluster != "F") {
      print("set point cluster as root")
      my_cds <- orderCells(my_cds, root_state = pointRoot(my_cds, rootcluster))
    }
    
    if (rootstate != "F") {
      print("set point state as root")
      my_cds <- orderCells(my_cds, root_state = as.numeric(rootstate))
    }
  }
  
  #### 基于之前的cds翻转轨迹
  if (do == "reverse") {
    print("reverse Pseudotim")
    my_cds <- orderCells (my_cds, reverse = TRUE)
  }
  
  ##5.8 对Pseudotime列进行标准化（离差标准化）
  print("  5.8 normalize Pseudotime (min/max)#############################")
  pData(my_cds)$Pseudotime <- (pData(my_cds)$Pseudotime - min(pData(my_cds)$Pseudotime))*10 / (max(pData(my_cds)$Pseudotime) - min(pData(my_cds)$Pseudotime))
  pData(my_cds)$Pseudotime <- round(pData(my_cds)$Pseudotime, 1)
  
  ##6 将Seurat对象中cluster列替换成state
  print("6 reshape seurat rds #############################")
  state <- rownames_to_column(pData(my_cds), var = "cell_id") %>%
    dplyr::select(cell_id, State)
  
  pro@meta.data <- rownames_to_column(pro@meta.data, var = "cell_id") %>%
    dplyr::select(-State) %>%
    left_join(state, by = "cell_id") %>%
    column_to_rownames(var = "cell_id") %>%
    mutate(cluster_old = cluster, cluster = State)
  
  Idents(pro) <- "cluster"
  
  saveRDS(pro, str_c(outdir, '/', prefix, '.cluster_as_state.rds'))
  
  #7 可视化——轨迹图输出
  print("7 visualization of trajectory #############################")
  cell_number <- length(rownames(pro@meta.data))
  cell_size<-3

  if(cell_number>1000){
     cell_size<-2.5
     }
  if(cell_number>5000){
     cell_size<-1.25
     }
  if(cell_number>10000){
     cell_size<-0.6
     }
  ##7.1 根据图例设置图片尺寸
  print("  7.1 set picture size #############################")
  cluster_num <- length(unique(as.character(pData(my_cds)$cluster)))
  if(cluster_num > 10){
    ncol <- 2
    ratio <- 2
  }else{
    ncol <- 1
    ratio <- 16 / 9
  }
  ratio_1 <- ifelse(ceiling(sqrt(cluster_num)) > 2, ceiling(sqrt(cluster_num)) / 2, 1)
  
  ##7.2 降维轨迹图
  print("  7.2 visualization of trajectory of not facet #############################")
  ###7.2.1 cluster
  print("   7.2.1 cluster #############################")
    
  ###7.2.2 facet_cluster
  print("   7.2.2 visualization of trajectory of facet #############################")
    
  ##7.3 频率分布
  print("  7.3 frequency distribution #############################")
  
  ###7.3.1 密度图——折线
  print("   7.3.1 frequency line #############################")
  kk = pData(my_cds)
  k = kk[c('cluster','Pseudotime')]
   
  ###7.3.2 流图
  print("   7.3.2 frequency stream #############################")
   ####7.3.2.1 sample
  print("    7.3.2.1 frequency stream sample #############################")
   ####7.3.2.2 group 
  print("    7.3.2.2 frequency stream group #############################")
   
  #8 提取保存可用信息
  print("8 save some useful data #############################")
  
  ##8.1 pData信息
  pdata <- rownames_to_column(pData(my_cds), var = "barcode") %>%
    rename(c(cluster = "Cluster"))
  
  write.table(pdata, file = str_c(outdir, '/', prefix, '.monocle_Pseudotime.xls'), quote = F, sep = "\t", row.names = F)
  
  ##8.2 monocle降维信息
  ll <- t(my_cds@reducedDimS)
  cell_s <- data.frame(cellID=rownames(ll), ll)
  write.table(cell_s, file = str_c(outdir, '/', prefix, '.plotpseudotimedata.xls'), quote = F, sep = "\t", row.names=F)
  
  ##8.3 
  line <- t(my_cds@reducedDimK)
  lines <- data.frame(ID=rownames(line), line)
  write.table(lines, file=paste(outdir, '/', prefix, '.plotlines.xls', sep = ''), quote = F, sep = "\t", row.names = F)
  
  #9 拟时序差异分析
  print("9 diff analysis of Pseudotime #############################")
  
  ##9.1 差异分析
  print("  9.1 diff #############################")
  my_pseudotime_de <- differentialGeneTest(my_cds, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = core)
  
  ##9.2 获取top基因
  print("  9.2 get top8/30 gene #############################")
  
  ###9.2.1 用于基因随细胞状态等的表达变化/monocle降维体系
  print("   9.2.1 top8 gene #############################")
  gene_to_pseudotime <- topGene(my_pseudotime_de, toPseudotime = "T", toHeatmap = "F")
  
  ###9.2.2 用于热图展示
  print("   9.2.2 top30 gene #############################")
  gene_to_cluster <- topGene(my_pseudotime_de, toPseudotime = "F", toHeatmap = "T")
  
  ##9.3 数据保存
  print("  9.3 save cds data #############################")
  ###9.3.1 保存cds
  saveRDS(my_cds, file = str_c(outdir, '/', prefix, '.pseudotime.cds'))
  
  ###9.3.2 拟时序差异基因信息
  write.table(my_pseudotime_de, file = str_c(outdir, '/', prefix, '.pseudotime_gene.xls'), quote=F, sep="\t", row.names=F)
  
  ##9.4 可视化
  print("  9.4 visualization of diff #############################")
  
  ###9.4.1 基因随细胞状态等的表达变化——top8
  print("   9.4.1 visualization of diff top8 #############################")
    ##单张出图
  reductions_genes <- str_c(feedback, "/reductions_genes")
  dir.create(reductions_genes)
  
    
  ###9.4.2 热图——top30
  print("   9.4.2 visualization of diff top30 #######################")
    
  ####top30基因聚类情况
  gene_name_clust <- as.data.frame(sort(cutree(p$tree_row, k = length(unique(pData(my_cds)$cluster)))))
  colnames(gene_name_clust) <- "Cluster"
  write.table(gene_name_clust, file = str_c(outdir, "/", prefix, ".rowgene.xls"), sep='\t', quote=F, row.names=T, col.names = NA)
  
  #copy readme and finished
  #file.copy(from='/SGRNJ03/PiplineTest01/Pipline_test/xudaochao/13.monocle/monocle_README.pdf', to = feedback)
  file.copy(from='/Public/Script/shouhou/SCRIPT/monocle2.22.0/script/monocleUpstream.pdf', to = feedback)
} 
  if (do == "ucell") {
    #setscore <- read.table(argv$genesetscore,header=T,sep="\t",check.names=F,comment.char='')
    setscore <- readRDS(argv$genesetscore)
    setscore <- setscore@meta.data
    ll <- gsub(" ","_",colnames(setscore))
    ll <- gsub("-","_",ll)
    print(ll)
    colnames(setscore) <- ll
    pathway <- unlist(strsplit(argv$geneset,split = ","))
    #pathway <- paste0(pathway,"_UCell")
    pathway <- gsub(" ","_",pathway)
    pathway <- gsub("-","_",pathway)
    print(head(setscore)) 
    #pathway <- colnames(setscore)
    #pathway <- pathway[-(1:11)]
    pdata <- pData(my_cds)
    setscore <- setscore %>% select(pathway)
    pdata <- merge(pdata,setscore,by = 'row.names') %>% column_to_rownames('Row.names')
    pData(my_cds) <- pdata
    print(head(pData(my_cds)))
    if (do == "state"){
  state <- unlist(strsplit(state,split=","))
  substateheatmap(my_cds,state,feedback,prefix)}
  if (do == "beam"){
  #if (state!='F'){
  #state <- unlist(strsplit(state,split=","))
  #substateheatmap(my_cds,state,feedback,prefix)}
  if (branch!='F'){
  #ordergene <- load((str_c(inputdir, file_name[str_detect(file_name, "ordergene.Rdata$")])))
  BEAM_res <- BEAM(my_cds, branch_point = branch, cores = core)
  #BEAM_res <- BEAM(my_cds[ordergene,], branch_point = branch, cores = 2) #只使用高边基因，使用全部基因太耗时
  
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  #saveRDS(BEAM_res,str_c(outdir, '/', prefix,"_BEAM_res.rds")
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  
  if (branch_cluster == "default") {
    num_cluster = length(levels(pro))
  } else {
    num_cluster = as.numeric(branch_cluster)
  }
  
    
  write.table(branched_heatmap_gene, file = str_c(outdir, "/", prefix, "_branched_heatmap_gene.xls"), sep='\t', quote=F, row.names=T, col.names = NA)
 }
} 
