#!/usr/bin/env Rscript
suppressMessages({
  library(argparser)
  library(Seurat)
  library(tidyverse)
  library(grid)
  library(cowplot)
  library(monocle)
  library(patchwork)
  library(ggsci)
  #library(RhpcBLASctl)
  library(ggstream)
  library(reshape)
  library(pheatmap)
})


#1 创建参数解释器
argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help = "the cell rds")
argv <- add_argument(argv,"--subcluster", help="cluster list which use to analysis, split by ,", default = 'all')
argv <- add_argument(argv,"--subsample", help="sample list which use to analysis, split by ,", default = 'all')
argv <- add_argument(argv,"--subset", help="whether to build a subset of the rdsData, optional:'F','default:Default standard','number:How many cells are retained per cell type?'", default = "default")
argv <- add_argument(argv,"--clusterlevels", help = "the cluster levels", default = 'F')
argv <- add_argument(argv,"--sample", help="sample list which use to add group/platform", default = 'F')
argv <- add_argument(argv,"--group", help="group list which is match with sample", default = 'F')
argv <- add_argument(argv,"--orderinglist", help="'provide filelist' or 'seurat: use seurat to find variable features' or 'persion: use persion of monocle'", default = "seurat")
argv <- add_argument(argv,"--batch", help="remove batch effect, optional:'F', 'sample', 'platform'", default="F")
argv <- add_argument(argv,"--platform", help="platform which is match with sample", default = 'F')
argv <- add_argument(argv,"--newcol", help="useing new color protocol", default = 'T')
argv <- add_argument(argv,"--colorset", help="useing  color protocol of metadata", default = 'T')
argv <- add_argument(argv,"--hvgnum", help=" high variablegenes  number(default)\n#Ordering based on genes that differ between clusters\n#Selecting genes with high dispersion across cells\n#Ordering cells using known marker genes\n", default = 2000)
argv <- add_argument(argv,"--reversePseudotime", help=" reverse Pseudotime", default = "F")
argv <- add_argument(argv,"--reduction_method", help="DDRTree,ICA,tSNE,SimplePPT,L1-graph,SGL-tree", default = "DDRTree")
argv <- add_argument(argv,"--rootcluster", help="the root cluster", default="F")
argv <- add_argument(argv,"--show_tree", help="Trajectories show tree", default = TRUE)
argv <- add_argument(argv,"--show_branch_points", help="Trajectories show tree", default = FALSE)
argv <- add_argument(argv,"--prefix", help="output prefix")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- parse_args(argv)
#1.1source
dirnameScript <- function(){
    # get full directory path of current script located
    cmd = commandArgs(trailingOnly = FALSE)
    scriptName = sub('--file=', "", cmd[grep('^--file=', cmd)])
    if (length(scriptName) > 0) {
        path = normalizePath(scriptName)
        dirname = dirname(path)
    } else {
        print('Warning: not a runtime environment, using current directory instead.')
        dirname = getwd()
    }
    return(dirname)
}
sourceScript <- function(x, dir='') {
    # try source script from directory one by one: file itself, dir from argument, script directory, libR
    scriptdir = dirnameScript()
    libRdir = file.path(scriptdir,'..','libR')
    searchpath = c(x,
                    file.path(dir, x),
                    file.path(scriptdir, x),
                    file.path(libRdir, x))
    sourced = FALSE
    for (i in searchpath) {
        if (file.exists(i)) {source(i); print(i); sourced = TRUE; break}
    }
    if (!sourced) {stop(paste0('can not source: ', x))}
}
scriptdir = dirnameScript()
sourceScript("monocle_use1.R")
sourceScript("seurat_use.R")
sourceScript("seurat_function.R")
#2解析参数
rds <- argv$rds
subcluster <- unlist(strsplit(argv$subcluster, split = ","))
subsample <- unlist(strsplit(argv$subsample, split = ","))
subset <- argv$subset

if (str_detect(argv$clusterlevels, ",")) {
  clusterlevels <- unlist(strsplit(argv$clusterlevels, split = ","))
} else {
  clusterlevels <- argv$clusterlevels
}

sample <- unlist(strsplit(argv$sample, split = ","))
group <- unlist(strsplit(argv$group, split = ","))

if (argv$batch == "F") {
  batch <- argv$batch
} else {
  batch <- str_c("~", argv$batch)
}

platform <- unlist(strsplit(argv$platform, split = ","))
newcol <- argv$newcol
colorset <- argv$colorset
use_gnum <- as.numeric(argv$hvgnum)
reverse_pt <- argv$reversePseudotime
method_use <- as.character(argv$reduction_method)
rootcluster <- argv$rootcluster
show_tree_argv <- as.logical(argv$show_tree)
show_branch_points_argv <- as.logical(argv$show_branch_points)

prefix <- argv$prefix
outdir <- argv$outdir
dir.create(outdir)

feedback <- str_c(outdir, "/", prefix)   #反馈目录
dir.create(feedback)

#3 rds梳理
print("3 reshape rds #############################")

##3.1 读取数据
print("  3.1 read rds #############################")
pro <- readRDS(rds)
pro@meta.data$cluster <- pro@active.ident
#pro@meta.data[["cluster"]] = as.character(pro@meta.data[["cluster"]])
#pro@meta.data[["sample"]] = as.character(pro@meta.data[["sample"]])
print(str_c("the ident befor Process is ", str_c(sort(as.character(levels(pro))), collapse = "/ "), sep = ": "))
print(str_c("the sample befor Process is ", str_c(sort(as.character(unique(pro@meta.data$sample))), collapse = "/ "), sep = ": "))

##3.2 cluster/sample构建子集
print("  3.2 filter cluster/sample #############################")
#pro <- subsetSc(pro)
pro <- subsetSc(pro, subcluster = subcluster, subsample = subsample)

##3.3 减少数据量——对cluster保留指定数量细胞
print("  3.3 subset #############################")
pro <- subsetClustercount(pro)
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

##3.4 修改cluster因子水平
print("  3.4 set cluster levels #############################")
if (!(identical(clusterlevels, "F"))) {
  pro <- clusterLevels(pro, clusterlevels)
}

##3.5 添加样本分组信息
print("  3.5 add group #############################")
if (group != "F") {
  pro <- addColumn(pro, addgorup = "T", addplatform = "F")
}

##3.6 添加样本平台信息
print("  3.6 add platform #############################")
if (platform != "F") {
  pro <- addColumn(pro, addgorup = "F", addplatform = "T")
}

##3.7 默认数据
DefaultAssay(pro) <- "RNA"

#4 色板设置
print("4 color set #############################")
pro$cluster <- as.character(pro$cluster)
pro$sample <- as.character(pro$sample)
if (colorset != "F"){
    colorlist <- getColors(pro,meta=c("sample","group","cluster","raw_cluster"))
    if (length(colorlist)!=0){
        clustcol <- colorlist$cluster_colors
        clustcol <- as.data.frame(clustcol)
        clustcol <- rownames_to_column(clustcol,'cluster')
        clustcol <- clustcol %>%
        mutate(cluster=fct_relevel(cluster,levels(pro))) %>%
        arrange(cluster)
        clustcol <- column_to_rownames(clustcol,'cluster')
        a <- as.character(clustcol$clustcol)
        names(a) <- rownames(clustcol)
        clustcol <- a
        groupcol <- colorlist$group_colors
        groupcol <- groupcol[!is.na(names(groupcol))]
        samplecol <- colorlist$sample_colors
        samplecol <- samplecol[!is.na(names(samplecol))]
    }else{
       clustcol <- colorSet(newcol)
       groupcol <- clustcol
       samplecol <- clustcol
   }
}else{
clustcol <- colorSet(newcol)
groupcol <- clustcol
samplecol <- clustcol
}
color_use <- clustcol[1:length(unique(pro@meta.data$cluster))]
groupcol_use <- groupcol[1:length(unique(pro@meta.data$group))]
samplecol_use <- samplecol[1:length(unique(pro@meta.data$sample))]
names(color_use) <- levels(pro)

#5 monocle分析
print("5 monocle analysis start #############################")

##5.1 利用monocle2构建S4对象
print("  5.1 creat monoclo data set #############################")
###5.1.1 表达量矩阵exprs数值矩阵：行名是基因，列名是细胞编号
count <- GetAssayData(object = pro[["RNA"]], slot = "counts")
#使用counts矩阵效果最好，monocle内部会做所有标准化的分析，如果数据已经标准化可能会中断

###5.1.2 基因注释信息——featureData: 第一列是基因编号, 其他列是基因对应的信息，列名必须为“gene_short_name”
#fdata <- as.data.frame(rownames(pro)); colnames(fdata) <- "gene_short_name"; rownames(fdata) <- fdata$gene_short_name

fdata <- data.frame(gene_short_name = row.names(count), row.names = row.names(count))

###5.1.3 细胞的表型信息phenoData: 第一列是细胞编号，其他列是细胞的相关信息
pdata <- pro@meta.data

###5.1.4 构建S4对象
my_cds <- newCellDataSet(count,
                         featureData = new("AnnotatedDataFrame", data =fdata),
                         phenoData = new("AnnotatedDataFrame", data = pdata),
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

##5.2 标准化、scale及筛选
print("  5.2 normalize & scale & detectGenes #############################")
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)
my_cds <- detectGenes(my_cds, min_expr = 0.1)
pData(my_cds)$UMI <- Matrix::colSums(exprs(my_cds))

##5.3 高变基因
print("  5.3 find variable features #############################")
if(argv$orderinglist == "seurat"){
  pro <- FindVariableFeatures(pro, selection.method = "vst", nfeatures = 2000)
  ordergene <- head(VariableFeatures(pro), use_gnum)
} else if (argv$orderinglist == "persion") {
  disp_table <- dispersionTable(my_cds)
  ordergene <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
} else{
  my_pseudotime_GENE <- read.table(argv$orderinglist, header=F,sep="\t")
  ordergene <- as.character(my_pseudotime_GENE$V1)
}
print(length(ordergene))

##5.4 标记排序基因
print("  5.4 lable ordergene #############################")
my_cds <- setOrderingFilter(my_cds,ordergene)
#p <- plot_ordering_genes(my_cds)
#pdf(str_c(outdir, '/', prefix, '.ordering_genes.pdf'))
#print(p)
#dev.off()
#png(str_c(outdir, '/', prefix, '.ordering_genes.png'))
#print(p)
#dev.off()

##5.5 降维——是否去批次，及去批次类型（group/platform）
print("  5.5 rm batch or not  #############################")
if (batch == "F") {
  print("this analysis was not rm batch")
  my_cds <- reduceDimension(my_cds, reduction_method = method_use)
} else {
  print(str_c("rm batch of ", batch))
  my_cds <- reduceDimension(my_cds, reduction_method = method_use, residualModelFormulaStr = batch)
}

##5.6 拟时间轴轨迹构建和在拟时间内排列细胞
print("  5.6 order & reverse or not #############################")
if (reverse_pt == "F") {
  print("not reverse")
  my_cds <- orderCells(my_cds)
} else {
  print("reverse")
  my_cds <- orderCells (my_cds, reverse = TRUE)
}

##5.7 指定root
print("  5.7 set root or not #############################")
if (rootcluster != "F") {
  print("point root")
  my_cds <- orderCells(my_cds, root_state = pointRoot(my_cds, rootcluster))
} else {
  my_cds <- my_cds
}

##5.8 对Pseudotime列进行标准化（离差标准化）
print("  5.8 normalize Pseudotime (min/max)#############################")
pData(my_cds)$Pseudotime <- (pData(my_cds)$Pseudotime - min(pData(my_cds)$Pseudotime))*10 / (max(pData(my_cds)$Pseudotime) - min(pData(my_cds)$Pseudotime))
pData(my_cds)$Pseudotime <- round(pData(my_cds)$Pseudotime, 1)

##6 将Seurat对象中cluster列替换成state
print("6 reshape seurat rds #############################")
state <- rownames_to_column(pData(my_cds), var = "cell_id") %>%
  dplyr::select(cell_id, State, Pseudotime)

pro@meta.data <- rownames_to_column(pro@meta.data, var = "cell_id") %>%
  left_join(state, by = "cell_id") %>%
  column_to_rownames(var = "cell_id") %>%
  mutate(cluster_old = cluster, cluster = State)

Idents(pro) <- "cluster"

saveRDS(pro, str_c(outdir, '/', prefix, '.cluster_as_state.rds'))
###重新指定cds的clusterlevels
if (!(identical(clusterlevels, "F"))) {
pData(my_cds)$cluster <- factor(pData(my_cds)$cluster,levels=clusterlevels)
}
#7 可视化——轨迹图输出
print("7 visualization of trajectory #############################")

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
p <- reductionTinotfacet(my_cds)

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
my_pseudotime_de <- differentialGeneTest(my_cds, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 2)

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
gene_name_clust <- as.data.frame(sort(cutree(p$tree_row, k = length(unique(pData(my_cds)$State)))))
colnames(gene_name_clust) <- "Cluster"
write.table(gene_name_clust, file = str_c(outdir, "/", prefix, ".rowgene.xls"), sep='\t', quote=F, row.names=T, col.names = NA)
