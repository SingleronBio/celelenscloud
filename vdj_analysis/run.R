#1 加载软件包
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
scriptdir = dirnameScript()
sourceScript <- function(x, dir='/Public/Script/shouhou/Integration_Platform/libR') {
    # try source script from directory one by one: file itself, dir from argument, script directory, libR
    # scriptdir = dirnameScript()
    scriptdir = scriptdir # 多次source时 在函数内部执行dirnameScript()会报错，所以这里用全局变量scriptdir替换了>原来的dirnameScript()
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

#gwd <- getwd()
suppressMessages({
  for (i in dplyr::setdiff(list.files(scriptdir), c("run.R","CDR3_README.docx","CDR3_README.pdf","VDJ_README.docx","VDJ_README.pdf"))) {
	print(i)
    sourceScript(i)
  }
  library(Seurat)
  library(venn)
  library(tidyverse)
  library(ggpubr)
  library(reshape2)
  library(argparser)
  sourceScript("heatmap_function.R")
})
#setwd(gwd)

#2 创建参数解析器
argv <- arg_parser("")
argv <- add_argument(argv, "--rds", help = "the rds file")
argv <- add_argument(argv, "--infodata", help = "", default = "F")
argv <- add_argument(argv, "--fresplit", help = "", default = "0, 1, 10, inf")
argv <- add_argument(argv, "--fresplitname", help = "", default = "Single,Medium,Large")
argv <- add_argument(argv, "--spname", help = "split by ,")
argv <- add_argument(argv, "--cluster", help = "cluster list ,split by ,", default = "all")
argv <- add_argument(argv, "--sample", help="sample list ,split by ,")
argv <- add_argument(argv, "--group", help="group list ,split by ,")
argv <- add_argument(argv, "--vdjlist", help="vdj result file list ,split by ,")
argv <- add_argument(argv, "--prefix", help="group name prefix")
argv <- add_argument(argv, "--outdir", help="the output dir")
argv <- add_argument(argv, "--vdjType", help="the VDJ type,TCR/BCR")
argv <- add_argument(argv, "--cloneDefine", help="use VDJ gene or CDR3 define clone,chose VDJ, CDR3, ,VDJ_CDR3, or VDJ_CDR3aa", default = "VDJ")
argv <- add_argument(argv, "--topN", help="the num of display clonetype", default = 5)
argv <- add_argument(argv, "--chairPair", help="keep pair chain", default = TRUE)
argv <- add_argument(argv, "--deal_barcode",help = "if deal the -1 sufix in vdjlist",default = "no")
argv <- parse_args(argv)

#3 解析参数
rds <- argv$rds
cluster <-  unlist(strsplit(argv$cluster, split = ","))
fresplit <-  as.numeric(unlist(strsplit(argv$fresplit, split = ",")))
fresplitname <-  unlist(strsplit(argv$fresplitname, split = ","))
vdjType <- argv$vdjType
cloneDefine <-  argv$cloneDefine
chairPair <- argv$chairPair
pt_use <- 0.1
top <- argv$topN
prefix <- argv$prefix
outdir <- argv$outdir

if (str_detect(outdir, "/$")) {
  outdir <- str_sub(outdir, 1, str_length(outdir)-1)
} else {
  outdir <- outdir
}

if (identical(argv$infodata, "F")) {
  spname <- unlist(strsplit(argv$spname, split = ","))
  samples <-  unlist(strsplit(argv$sample, split = ","))
  groups <- unlist(strsplit(argv$group, split = ","))
  vdjlist <- unlist(strsplit(argv$vdjlist, split=","))
} else {
  infodata <- read.table(argv$infodata, sep = "\t", col.names = c("spname", "samples", "groups", "vdjlist"), stringsAsFactors = F)
  spname <- infodata$spname
  samples <- infodata$samples
  groups <- infodata$groups
  vdjlist <- infodata$vdjlist
}

#4 创建反馈目录
##4.1 QC目录
dir.create(paste0(outdir,"/00.QC/","01.QC.raw"),recursive =TRUE)
dir.create(paste0(outdir,"/00.QC/","02.QC.productive"),recursive =TRUE)
dir.create(paste0(outdir,"/00.QC/","03.QC.final"),recursive =TRUE)
dir.create(paste0(outdir,"/00.QC/","04.QC.chain"),recursive =TRUE)

##4.2 feedback
dir.create(paste0(outdir,"/00.Not_Feedback"),recursive =TRUE)

##4.3 分析结果目录
dir.create(paste0(outdir,"/01.sc.CloneBar"),recursive =TRUE)
dir.create(paste0(outdir,"/02.sc.Diversity"),recursive =TRUE)
dir.create(paste0(outdir,"/03.sc.Clustertree"),recursive =TRUE)
dir.create(paste0(outdir,"/04.sc.HaveOrNot"),recursive =TRUE)
dir.create(paste0(outdir,"/05.sc.Otherinfo"),recursive =TRUE)

if (cloneDefine == "VDJ"|cloneDefine == "VDJ_CDR3"|cloneDefine == "VDJ_CDR3aa") {
  dir.create(paste0(outdir,"/06.sc.CloneTypePie"), recursive =TRUE)
  dir.create(paste0(outdir,"/07.sc.UmapCloneType"),recursive =TRUE)
  dir.create(paste0(outdir,"/08.sc.UpsetR"),recursive =TRUE)
  dir.create(paste0(outdir,"/09.sc.VJproportionUsage"),recursive =TRUE)
  dir.create(paste0(outdir,"/10.sc.TopNgenexpTop"),recursive =TRUE)
  dir.create(paste0(outdir,"/11.sc.VJpairs"),recursive =TRUE)
}

#5 读取rds
print("step5 ######################################")
pro <- readRDS(rds)
pro@meta.data$cluster <- pro@active.ident
pro@meta.data$cluster <- as.character(pro@meta.data$cluster)
pro@meta.data$sample <- as.character(pro@meta.data$sample)
if(argv$cluster == "all"){
	cluster <- levels(pro)
}
#############从rds里读取颜色信息，如果没有就用默认值
color_list <- getColors(pro,meta=c("sample","group","cluster"))
if(!is.na(color_list[["cluster_colors"]][1])){
	cluster_colors <- color_list[["cluster_colors"]]
	print(cluster_colors)
	cluster_colors <- cluster_colors[cluster]
}else{
	print("No cluster colors information in rds.")
	cluster_colors <- NA
}
###
if(!is.na(color_list[["group_colors"]][1])){
    group_colors <- color_list[["group_colors"]]
}else{
	print("No group colors information in rds.")
    group_colors <- NA
}
###
if(!is.na(color_list[["sample_colors"]][1])){
    sample_colors <- color_list[["sample_colors"]]
}else{
	print("No sample colors information in rds.")
    sample_colors <- NA
}

############# 
#6 检查输入的样本、分组信息数量是否与免疫组数据数量一致，防止手误多写、少写(step1.judge.R-1)
print("step6 ######################################")
cheackSampleCount(pro@meta.data)

#7 rds构建子集及添加样本分组信息（step2.seuratProjectprepart.R）
print("step7 ######################################")

##7.1 筛选cluster/ sample
print(" step7.1 ######################################")
pro <- subsetSc()
cluster <- levels(pro)

##7.2 添加group列
print(" step7.2######################################")
pro <- addColumngroup()

#8 读取免疫组数据，并存入list（step3.vdjListprepart.R-1）
print("step8 ######################################")
contiglist <- mergeSample(vdjlist,deal_barcode = argv$deal_barcode)

#9 在免疫组标准分析文件baecode前加上文库前缀（step3.vdjListprepart.R-2）
print("step9 ######################################")
contiglist <- addBarcodeprefix(contiglist, spname)

#10 添加VDJgene列，以“-”连接v_gene、d_gene、j_gene（step3.vdjListprepart.R-3）
print("step10 ######################################")
newcontiglist <- addVDJgenecol(contiglist)

#11 免疫组数据筛选
print("step11 ######################################")

##11.1 原始数据克隆型频率展示_QA(step4.qc.R-1)
print(" step11.1 ######################################")
rawQC.picture <- dataQC(newcontiglist, vdjType)

for(i in 1:length(rawQC.picture)) {
  ggsave(rawQC.picture[[i]],file=paste0(outdir,"/","00.QC/01.QC.raw","/",samples[i],".raw_chain.freq.png"), height=8, width=14)
}

##11.2 按照我们的标准筛选，筛选标准见confluence(step5.fiterVDJ.R)
print(" step11.2 ######################################")

###11.2.1 数据筛选及Productive过滤（1、Productive过滤后绘图展示；2、按照我们的标准过滤）
print("  step11.2.1 ######################################")
newcontig.filter <- filterVDJ(newcontiglist, vdjType = vdjType, chairPair = chairPair)

###11.2.2 筛选后的数据-list
print("  step11.2.2 ######################################")
newcontiglist <- newcontig.filter[[1]]
save(newcontiglist,file = paste0(outdir,"/00.Not_Feedback/filter_10X_list.RData"))

###11.2.3 筛选Productive过滤后克隆型可视化结果
print("  step11.2.3 ######################################")
plist.productive <- newcontig.filter[[2]]
for(i in 1:length(plist.productive)) {
  ggsave(plot=plist.productive[[i]],file=paste0(outdir,"/","00.QC/02.QC.productive","/",samples[i],".productive_chain.freq.png"), height=8, width=14)
}

###11.2.4 筛选后QA(step4.qc.R-1)
print("  step11.2.4 ######################################")
plist.final <- dataQC(newcontiglist,vdjType)
for(i in 1:length(plist.final)) {
  ggsave(plist.final[[i]],file=paste0(outdir,"/","00.QC/03.QC.final","/",samples[i],".final_chain.freq.png"), height=8, width=14)
}

#12 数据整理——之前一个细胞的两条链信息分两行展示，此处将一个细胞两条链用一行展示，并添加一些行新列，后续调用(step6.combined.R-1)
print("step12 ######################################")
if(vdjType=="TCR"){
  combined <- combinedTCR(newcontiglist, samples, groups)
}else{
  combined <- combinedBCR(newcontiglist, samples, groups)
}

#13 保存经过塑型及过滤的免疫组数据(list)，方便后续报错排查原因，避免手动输入参数——脚本维护人员使用
print("step13 ######################################")
save(combined, file =paste0(outdir, "/00.Not_Feedback/combined.Rdata"))

#14 按样本绘制T/BCR单双链细胞数(step4.qc.R-2)
print("step14 ######################################")
data <- pairChainCells(combined, samples)
dataplot <- melt(data,id=c("sample"))
dataplot$value <- as.numeric(dataplot$value)

p <- ggplot(dataplot, aes( x = sample, y = value, fill = variable))+
  geom_col(position = 'stack', width = 0.6) + 
  theme_bw() +
  labs(x="Samples",y="Counts") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

write.table(data, file=paste0(outdir, "/00.QC/04.QC.chain/", prefix, ".pair.single.num.csv"), sep=",", row.names=F)
ggsave(p, file=paste0(outdir, "/00.QC/04.QC.chain/", prefix, ".pair.single.num.png"))
ggsave(p, file=paste0(outdir, "/00.QC/04.QC.chain/", prefix, ".pair.single.num.pdf"))

#15 检查输入样本名称是否与转录组barcode样本前缀一致(step1.judge.R-2)
print("step15 ######################################")
if (cheackSample(pro@meta.data)){
  print("samples which input was the same as prefix of barcode of scRNA")
} else {
  print("error: samples which input was not the same as prefix of barcode of scRNA")
  quit()
}

#16 将combined合并到meta.data(step6.combined.R-2)
print("step16 ######################################")
t <- combinedJoinTometa.data(combined, pro)

##步骤8处理时像数据框中加入了ID/ groups列，为了数据整洁度更高，此处删除
#t@meta.data <- dplyr::select(t@meta.data, -c("groups"))

#17 按样本和组统计克隆型频率(step7.CountClonalFrequency.R)
print("step17 ######################################")

##17.1 按样本统计克隆型频率——Frequency
print(" step17.1 ######################################")
sampleFrequency <- CountClonalFrequency(t, by = "sample", cloneDefine = cloneDefine)
sc.combinedSample <- sampleFrequency[[1]]
t <- sampleFrequency[[2]]

##17.2 按组统计克隆型频率——group_freq
print(" step17.2 ######################################")
if (identical(samples, groups)) {
  
} else {
  groupFrequency <- CountClonalFrequency(t, by = "group", cloneDefine = cloneDefine)
  sc.combinedGroup <- groupFrequency[[1]]
  t <- groupFrequency[[2]]
}

#t@meta.data <- mutate(t@meta.data, clontype = if_else(is.na(t@meta.data$CTgene), "na", "clontype"))

#18 保存数据(Rdata/rds)
print("step18 ######################################")
save(sc.combinedSample, file = str_c(outdir, "/00.Not_Feedback/sc.sampleList.Rdata"))

env <- ls()
if (!(identical(samples, groups))) {
  save(sc.combinedGroup, file = str_c(outdir, "/00.Not_Feedback/sc.groupList.Rdata"))
}

saveRDS(t, file = str_c(outdir, "/00.Not_Feedback/", prefix, ".rds"))

#19 可视化——cloneDefine为VDJ及CDR3都能输出的图
print("step19 ######################################")

##19.1 sc.ClonepBar_占比柱状图(picture1.clonBar.R)
print(" step19.1 ######################################")
###19.1.1 按样本绘图
cloneBarSample(sc.combinedSample)
statClonetype_percentage(sc.combinedSample)
###19.1.2 按样本分组画图
if (!(identical(samples, groups))) {
  cloneBarGroup(sc.combinedGroup)
}

##19.2 sc.Diversity_多态性柱状图(picture2_3.immunarch.R)
print(" step19.2 ######################################")
immunarch_data <- transform_data_for_immunarch(sc.combinedSample, vdjType, cloneDefine)

if(!all(as.character(immunarch_data$meta$Sample) == as.character(immunarch_data$meta$Group))){
  plot_diversity(immunarch_data, outdir = outdir, plot_by = "Group",color = group_colors)
}

plot_diversity(immunarch_data, outdir = outdir, plot_by = "Sample", test = FALSE,color = sample_colors)

##19.3 sc.Clustertree_样本聚类(picture2_3.immunarch.R)
print(" step19.3 ######################################")
if(length(immunarch_data$meta$Sample) > 2){
	plot_hclust(immunarch_data, prefix = prefix, outdir = outdir)
}

##19.4 sc.HaveOrNot_柱状图(picture4.HaveOrNot.R)

###19.4.1 按样本绘图
haveOrNotsample(sc.combinedSample)

###19.4.2 按样本分组画图
if (!(identical(samples, groups))) {
  haveOrNotGroup(sc.combinedGroup)
}

##19.5 sc.Otherinfo
print(" step19.5 ######################################")
getInfo(sc.combinedSample, split = "sample")
if (!(identical(samples, groups))) {
  getInfo(sc.combinedGroup, split = "group")
}
#stop("I'm testing script...")
#20 cloneDefine仅为VDJ时，输出的图
print("step20 ######################################")

if (cloneDefine == "VDJ"| cloneDefine == "VDJ_CDR3"| cloneDefine == "VDJ_CDR3aa") {
  
  ##20.1 sc.CloneType.Percent_点饼图(picture5.sc.clontype.multi.R)
  print(" step20.1 ######################################")
  
  ###20.1.1 样本出图
  print("  step20.1.1 ######################################")
  scClonetype(t, sample_p = "T", group_p = "F")
  
  ###20.1.2样本分组出图
  print("  step20.1.2 ######################################")
  if (identical(samples, groups)) {
    
  } else {
    t_dotpie <- t
    t_dotpie@meta.data <- rownames_to_column(t_dotpie@meta.data, var = "cellid") %>%
      mutate(sample_old = sample, sample_freq = Frequency,
             sample = group, Frequency = group_freq) %>%
      column_to_rownames(var = "cellid")
    
    scClonetype(t_dotpie, sample_p = "F", group_p = "T")
  }
  
  ##20.2 umap上展示top(默认展示top5)克隆型（所有样本克隆型放一起筛选top）(picture6.umapClonetype.R-1)
  print(" step20.2 ######################################")
  topCloneUmap()
  topCloneUmap_Sample()
  ##20.3 umap上展示含克隆型的细胞(picture6.umapClonetype.R-2)
  print(" step20.3 ######################################")
  CloneUmap()
  
  ##20.4 upset(picture7.upsetR.)
  print(" step20.4 ######################################")
  
  ###20.4.1 按样本出图
  print("  step20.4.1 ######################################")
  upsetR (t, outdir, num=16, width =14, height =7, cloneDefine = cloneDefine,color = cluster_colors,vdjtype = vdjType)
  
  ###20.4.2 按样本分组出图
  print("  step20.4.2 ######################################")
  if (identical(samples, groups)) {
    
  } else {
    pro <- t
    pro@meta.data <- dplyr::select(pro@meta.data, -sample, -Frequency) %>% dplyr::rename(sample = group, Frequency = group_freq)
    
    upsetR (pro, outdir, num = 16, width = 14, height = 7, cloneDefine = cloneDefine,color = cluster_colors, vdjtype = vdjType)
  }
  
  ##20.5 Usage热图(picture8.VJusage.R)
  print(" step20.5 ######################################")
#  size <- figureSize(t)
  #samples <- intersect(samples,unique(t@meta.data$sample))
  for(sample in samples){
    cell.use <- rownames(t@meta.data[t@meta.data$sample %in% sample,])
    plot.rds <- subset(t,cells=cell.use)
    if(sum(is.na(plot.rds@meta.data$CTgene))>0) {
      scVDJUsage(plot.rds,datatype=vdjType,outdir=outdir)}
  }
  
  if (identical(samples, groups)) {
    
  } else {
    tmp <- t@meta.data$sample
    for(i in 1:length(samples)){
      tmp <- gsub(paste0("^",samples[i],"$"),groups[i],tmp)	
    }
    t@meta.data$group <- tmp
    for(group in unique(t@meta.data$group)){
      cell.use <- rownames(t@meta.data[t@meta.data$group %in% group,])
      plot.rds <- subset(t,cells=cell.use)
      if(sum(is.na(plot.rds@meta.data$CTgene))>0) {
        scVDJUsageGroup(plot.rds,datatype=vdjType,outdir=outdir)}
    }
  }
  
  ##20.6 sc.TopNgenexpTop
  print(" step20.6 ######################################")
  
  ###20.6.1 按样本画图
  levels_sp <- samples
  subprefix <- "sample"
  plotFreqPerSample(sc.combinedSample,color = sample_colors)
  
  ###20.6.2 按组画图
  if (identical(samples, groups)) {
    
  } else {
    levels_sp <- unique(groups)
    subprefix <- "group"
    for (i in 1:length(sc.combinedGroup)) {
      sc.combinedGroup[[i]]$sample <- sc.combinedGroup[[i]]$group
    }
    
    plotFreqPerSample(sc.combinedGroup,color = group_colors)
  }
  
  
  ##20.7 sc.VJpairs
  print(" step20.7 ######################################")
  
  ###20.7.1 按样本画图
  print("  step20.7.1 ######################################")
  plotVJHeatmap(sc.combinedSample)
  
  ###20.7.2 按组画图
  print("  step20.7.2 ######################################")
  if (identical(samples, groups)) {
    
  } else {
    plotVJHeatmap(sc.combinedGroup)
  }
  
  ##VDJ_README
  scriptdir = dirnameScript()
  readmePath <- file.path(scriptdir, '/VDJ_README.pdf')
  file.copy(from=readmePath, to=outdir)
  #file.copy(from="/Public/Script/shouhou/SCRIPT/SigVDJ/SigVDJ/SigVDJ_2.0_20220622/20220621/VDJ_README.pdf", to = outdir)
} else {
  scriptdir = dirnameScript()
  readmePath <- file.path(scriptdir, '/CDR3_README.pdf')
  file.copy(from=readmePath, to=outdir)
  ##CDR3_README
  #file.copy(from="/Public/Script/shouhou/SCRIPT/SigVDJ/SigVDJ/SigVDJ_2.0_20220622/20220621/CDR3_README.pdf", to = outdir)
}
print("All job finished!!!!")
