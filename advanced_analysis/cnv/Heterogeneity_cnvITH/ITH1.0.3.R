#!/usr/bin/env Rscript
library(argparser)
library(tibble)
library(dplyr)
library(readr)
library(Seurat)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(svglite)
# source("/SGRNJ03/pipline_test/fangsheng/seurat/R/seurat_function.R")
# source("/SGRNJ01/Public/Script/shouhou/SCRIPT/scTools/boxplot_function/boxplot_function.r")
# source("/Public/Script/shouhou/SCRIPT/seuratSingleron/seuratSingleron/R/seuratSingleron.R")

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="Seurat rds")
argv <- add_argument(argv,"--ncell_min", help="minimum number of cell per identity. If certain identity has cell number lower than ncell_min, the identity will be discarded.",default=20)
argv <- add_argument(argv,"--ncell_max", help="maximum number of cell per identity.",default=500)
argv <- add_argument(argv,"--outdir", help="output filepath",default=".")
argv <- add_argument(argv,"--cluster", help="subcluster list ,split by ,", default='all')
argv <- add_argument(argv,"--sample", help="subsample list ,split by ,", default='all')
argv <- add_argument(argv,"--type", help="figure x-axis for type,choose celltype,sample,group or HCL",default="celltype")
argv <- add_argument(argv,"--color", help="set color for figure",default="default")
argv <- add_argument(argv,"--order",help="x axis order",default="F")
argv <- add_argument(argv,"--method", help="select the method of calculating ITH:'dist or iqr'",default="dist")
argv <- add_argument(argv,"--compare", help="AvsB,AvsC, required when has multi groups, split by ,",default="F")
argv <- add_argument(argv,"--cmethod", help="t.test, wilcox.test", default='t.test')
argv <- add_argument(argv,"--prefix", help="prefix",default="")
argv <- add_argument(argv,"--HCL", help="HCL.barcode/HCL.withanno.barcode file path, only need for type is HCL",default="F")
argv <- parse_args(argv)


# 测试参数
########################################
# argv$rds <- "/SGRNJ06/NJ06Aftersales/2022/Cancer/P22091302/B1/module_analysis/result/Modules/1.ClustAnno/CellLabel/Main/v3/dalei_v4_to_v3.rds"
# argv$outdir <- "/SGRNJ03/pipline_test/chenkai/ITH/HCL_result"
# argv$type <- "HCL"
# argv$prefix <- "P22050607" 
# argv$HCL <- "/SGRNJ06/NJ06Aftersales/2022/Cancer/P22091302/B1/module_analysis/result/Modules/5.Heterogeneity/inferCNV/02.P22091302.infercnv/infercnv.HCL.withanno.barcode.txt"
# argv$order <- "StromalCells,TCells,HCL_1,HCL_3,HCL_2,HCL_4"
# argv$method <- "iqr"
# argv$cluster <- "StromalCells,TCells,HCL_1,HCL_2,HCL_3,HCL_4"
########################################

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

sourceScript <- function(x, dir='/Public/Script/shouhou/Integration_Platform/libR') {
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


scriptDirname <- dirnameScript()
# readmePath <- file.path(scriptDirname, "readme", "ITH-score_README.pdf")
sourceScript("seurat_function.R")
sourceScript("boxplot_function.r")
sourceScript("seuratSingleron.R")

prefix<- argv$prefix
method <- argv$method
print(paste0("method:",method))
type <- argv$type
print(paste0("type:",type))
rds = argv$rds
ncell.min = as.numeric(argv$ncell_min)
ncell.max = as.numeric(argv$ncell_max)
print(paste0("ncell.min:",ncell.min))
print(paste0("ncell.max:",ncell.max))
outdir <- argv$outdir
dir.create(outdir)
set.seed(0)
print ("loading rds...")
rds = readRDS(rds)
#subset sample or cluster
sample <- unlist(strsplit(argv$sample,split=","))
cluster <- unlist(strsplit(argv$cluster,split=","))

if (type == "HCL") {
  HCL <- read.table(argv$HCL, col.names=c("cell_id","HCL"))
  rds <- subset(rds, cell=HCL$cell_id)
  rds <- AddMetaData(rds, metadata=HCL$HCL, col.name="HCL")
  Idents(rds) <- rds$HCL
}

rds <- subsetRDS(rds, sample=sample, cluster=cluster)

if (type == "group"){
  idents = unique(rds@meta.data$group)
} else if (type == "sample"){
  idents = unique(rds@meta.data$sample)
} else if (type == "celltype"){
  type <- "cluster"
  idents = unique(Idents(rds))
} else if (type == "HCL") {
  idents = unique(rds@meta.data$HCL)
}


print(table(rds$sample))
print(table(rds@active.ident))
print(table(rds$group))
print(idents)

if (argv$color=="default"){
  if (paste0(type, "_colors") %in% colnames(rds@meta.data)) {
    color1 <- getColors(rds, type)
	color1[[1]] <- color1[[1]][!is.na(color1[[1]])]
  } else {
    color_old <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold", "DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4",  "#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B",  "#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080",  "#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49",  "#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD",  "#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B",  "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9",  "#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000",  "#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
    color_new <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
    color1 <- c(color_new,color_old)
    print(paste(paste0(type, "_colors") , "is not in metadata. we'll use new color"))
  }
} else{
  color1 <- unlist(strsplit(argv$color, split = ","))
}

  if(length(idents)<=7){
  fig_width <- 6 
} else {
  fig_width <- 6 + length(idents)*0.3
}

#---------------------------从这里开始-----------------------------------------


#ITH_iqr function
ITH_iqr <- function(rds,idents,type,ncell.min,ncell.max,color,outdir,prefix,fig_width){
df = tibble(ITH_score=numeric(),type = character())
dat = tibble(type=character(),iqr=numeric())
for (ident in idents){
   print(ident)
  if (type == "group"){
    Idents(rds) <- rds$group
    rds.ident = subset(rds,ident= ident,downsample = ncell.max)
    ncell.ident = length(colnames(rds.ident))
    print(table(rds.ident$group))
  } else if (type == "sample"){
    Idents(rds) <- rds$sample
    rds.ident = subset(rds,ident= ident,downsample = ncell.max)
    ncell.ident = length(colnames(rds.ident))
    print(table(rds.ident$sample))
  } else if (type == "cluster"){
    rds.ident = subset(rds,ident = ident,downsample = ncell.max)
    ncell.ident = length(colnames(rds.ident))
    print(table(rds.ident@active.ident))
  } else if (type == "HCL"){
    
    Idents(rds) <- rds$HCL
    rds.ident = subset(rds,ident= ident,downsample = ncell.max)
    ncell.ident = length(colnames(rds.ident))
    print(table(rds.ident$HCL))
    
  }
  if ( ncell.ident < ncell.min ){
      log = paste(ident, "has",ncell.ident,"cells. Discarded.")
      print (log)
      next
    }
    tryCatch({
      data = GetAssayData(object = rds.ident)
      corr = cor(as.matrix(data))
      melt_cor <- melt(corr)
      IQR = IQR(melt_cor$value)
      dat = add_row(dat,type=ident,iqr=IQR)
      value = as.vector(melt_cor$value)
      df.ident = tibble(ITH_score=value,type=ident)
      df = rbind(df,df.ident)
      print(paste(ident,"done"))},
    error = function(e){print(e)}
    )    
}

#df.med = dplyr::group_by(df,type) %>% dplyr::summarize(ITH_median=median(ITH_score))
#write.csv(dat,paste0(outdir,"/",prefix,"_ITH_score_iqr.csv"))
write.csv(data.frame(df),paste0(outdir,"/",prefix,"_ITH_score_all.csv"), row.names = F, quote=F)

df1 <- as.data.frame(df)
df.mean=group_by(df1, type) %>% summarize_each(funs(mean))
df.med=group_by(df1, type) %>% summarize_each(funs(median))
df.mean <- as.data.frame(df.mean)
colnames(df.mean) <- c("type","mean")
df.med <- as.data.frame(df.med)
colnames(df.med) <- c("type","median")
df.stat <- merge(df.mean,df.med,by="type")
write.csv(df.stat, paste0(outdir,"/",prefix,"_ITH_score_stat.csv"), row.names = TRUE)
gcompare <- list()
if(argv$compare!="F"){                                           # iqr组间比较
  compare <- unlist(strsplit(argv$compare,split=","))
    for(i in 1:length(compare)){
  gcompare[[i]] <-unlist(strsplit(compare[i], split = "vs"))
  }
  print(gcompare)
  for(i in 1:length(compare)){
    com <- compare[i]
    com <- unlist(strsplit(com,split="vs"))
    sub_data <- df[df$type %in% com,]
    sub_test <- compare_means(ITH_score~type,data = sub_data,method = argv$cmethod)
    if(i==1){
      stat_test <- sub_test
    }else{
      stat_test <- rbind(stat_test, sub_test)
    }
  }
  stat_test <- data.frame(stat_test)
  write.csv(stat_test,paste0(outdir,"/",prefix,"_ITH_score_sig.csv"))
  
  
  type1<- type
  if (argv$order!="F"){
    df$type <- factor(df$type,levels=unlist(strsplit(argv$order,",")))
	dat$type <- factor(dat$type,levels=unlist(strsplit(argv$order,",")))
  } else {
    df$type <- factor(df$type, levels = levels(rds@meta.data[,type]))
  }
  
	
  if (is.list(color)) {
    df <- as.data.frame(df)
    df$color <- color[[paste0(type1,"_colors")]][match(df$type,names(color[[paste0(type1,"_colors")]]))]  # 传颜色参数至df中
    
	box_color <- data.frame(color)
	box_color <- rownames_to_column(box_color, var="ident")
  box_color$ident <- factor(box_color$ident,levels=levels(df$type))
	box_color <- box_color[box_color$ident %in% unique(df$type),]

                  
  p1 <- boxplot_function(data=df, x="type", y="ITH_score", group="type", filltype="fill", usecolor=box_color[order(box_color$ident),2], savefig=FALSE)
  p1 <- p1 + stat_compare_means(comparisons = gcompare, label = "p.signif", method = argv$cmethod)
  } else {
    p1 <- boxplot_function(data=df, x="type", y="ITH_score", group="type", filltype="fill", savefig=FALSE)
    p1 <- p1 + stat_compare_means(comparisons = gcompare, label = "p.signif", method = argv$cmethod)
  }
  pdf(paste0(outdir,"/",prefix,"_ITH_score_boxplot.pdf"), width=fig_width, height = 8)
  print(p1)
  dev.off()
  svglite(paste0(outdir,"/",prefix,"_ITH_score_boxplot.svg"), width=fig_width, height = 8)
  print(p1)
  dev.off()
  png(paste0(outdir,"/",prefix,"_ITH_score_boxplot.png"), width=fig_width, height = 8 , units='in',res=400)
  print(p1)
  dev.off()
}else{                                                         #### iqr无组间比较
type1<- type
if (argv$order!="F"){
  df$type <- factor(df$type,levels=unlist(strsplit(argv$order,",")))
} else {
  df$type <- factor(df$type, levels = levels(rds@meta.data[,type]))
}

if (is.list(color)) {
  df <- as.data.frame(df)
  df$color <- color[[paste0(type1,"_colors")]][match(df$type,names(color[[paste0(type1,"_colors")]]))]  # 传颜色参数至df中
  
	box_color <- data.frame(color)
	box_color <- rownames_to_column(box_color, var="ident")
  box_color$ident <- factor(box_color$ident,levels=levels(df$type))
	box_color <- box_color[box_color$ident %in% unique(df$type),]
  
  p1 <- boxplot_function(data=df, x="type", y="ITH_score", group="type", filltype="fill", usecolor=box_color[order(box_color$ident),2], rmoutlier=TRUE, fig_save_name=paste0(outdir,"/",prefix,"_ITH_score_boxplot"))
} else {
  df <- as.data.frame(df)
  p1 <- boxplot_function(data=df, x="type", y="ITH_score", group="type", filltype="fill", rmoutlier=TRUE, fig_save_name=paste0(outdir,"/",prefix,"_ITH_score_boxplot"))
}

}


#--------------------------------------------------条形图-----------------------------------------------
# if (is.list(color)) {
  # dat$color <- color[[paste0(type1,"_colors")]][match(dat$type,names(color[[paste0(type1,"_colors")]]))]
  # p2 <- ggplot(dat,aes(x=type,y=iqr,fill=type)) + 
   # geom_bar(stat="identity") +                    
   # theme(panel.grid = element_blank(),
         # panel.background = element_rect(color = 'black', fill = 'transparent'),
         # legend.key = element_rect(fill = 'transparent')
   # ) + 
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
   # scale_fill_manual(values=dat$color)  + labs(x=type1) 
# } else {
  # p2 <- ggplot(dat,aes(x=type,y=iqr,fill=type)) + 
   # geom_bar(stat="identity") +                   
   # theme(panel.grid = element_blank(),
         # panel.background = element_rect(color = 'black', fill = 'transparent'),
         # legend.key = element_rect(fill = 'transparent')
   # ) + 
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
   # scale_fill_manual(values=color_new)  + labs(x=type1)
# }
# pdf(paste0(outdir,"/",prefix,"_ITH_score_iqr.pdf"), width=fig_width)
# print(p2)
# dev.off()
# png(paste0(outdir,"/",prefix,"_ITH_score_iqr.png"), width=fig_width, height=8, units='in',res=400)
# print(p2)
# dev.off()
#--------------------------------------------------条形图------------------------------------------------

}

#ITH_dist function
ITH_dist <- function(rds,idents,type,ncell.min,ncell.max,color,outdir,prefix,fig_width){
df = tibble(ITH_score=numeric(),type = character())
dat = tibble(type=character(),iqr=numeric())
for (ident in idents){
   print(ident)
  if (type == "group"){
    Idents(rds) <- rds$group
    rds.ident = subset(rds,ident= ident,downsample = ncell.max)
    ncell.ident = length(colnames(rds.ident))
    print(table(rds.ident$group))
  } else if (type == "sample"){
    Idents(rds) <- rds$sample
    rds.ident = subset(rds,ident= ident,downsample = ncell.max)
    ncell.ident = length(colnames(rds.ident))
    print(table(rds.ident$sample))
  } else if (type == "cluster"){
    rds.ident = subset(rds,ident = ident,downsample = ncell.max)
    ncell.ident = length(colnames(rds.ident))
    print(table(rds.ident@active.ident))
  } else if (type == "HCL"){
    
    Idents(rds) <- rds$HCL
    rds.ident = subset(rds,ident= ident,downsample = ncell.max)
    ncell.ident = length(colnames(rds.ident))
    print(table(rds.ident$HCL))
    
  }
  if ( ncell.ident < ncell.min ){
    log = paste(ident, "has",ncell.ident,"cells. Discarded.")
    print (log)
    next
  }
  tryCatch({
    rds.ident = FindVariableFeatures(rds.ident)
    rds.ident = ScaleData(rds.ident)
    rds.ident = RunPCA(object = rds.ident, npcs = 20,pc.genes = rds.ident@var.genes, do.print = FALSE)
    data = (rds.ident[["pca"]]@cell.embeddings)
    DIST = dist(data)
    dat = as.vector(DIST)
    df.ident = tibble(ITH_score=dat,type=ident)
    df = rbind(df,df.ident)
    print(paste(ident,"done"))},
    error = function(e){print(e)}
  )    
}
write.csv(data.frame(df),paste0(outdir,"/",prefix,"_ITH_score_all.csv"), row.names=F, quote=F)

df1 <- as.data.frame(df)
df.mean=group_by(df1, type) %>% summarize_each(funs(mean))
df.med=group_by(df1, type) %>% summarize_each(funs(median))
df.mean <- as.data.frame(df.mean)
colnames(df.mean) <- c("type","mean")
df.med <- as.data.frame(df.med)
colnames(df.med) <- c("type","median")
df.stat <- merge(df.mean,df.med,by="type")
write.csv(df.stat ,paste0(outdir,"/",prefix,"_ITH_score_stat.csv"))
gcompare <- list()
if(argv$compare!="F"){                                          # 绘制有组间比较的小提琴图
  compare <- unlist(strsplit(argv$compare,split=","))
  for(i in 1:length(compare)){
  gcompare[[i]] <-unlist(strsplit(compare[i], split = "vs"))
  }
  print(gcompare)
  for(i in 1:length(compare)){
    com <- compare[i]
    com <- unlist(strsplit(com,split="vs"))                     # com存放需要进行比较的两组的组名
    sub_data <- df[df$type %in% com,]
    sub_test <- compare_means(ITH_score~type,data = sub_data,method = argv$cmethod)
    if(i==1){
      stat_test <- sub_test
    }else{
      stat_test <- rbind(stat_test, sub_test)
    }
  }
  stat_test <- data.frame(stat_test)                            # stat_test存放组间比较的p值
  print(head(stat_test))
  write.csv(stat_test,paste0(outdir,"/",prefix,"_ITH_score_sig.csv"))   # 输出p值统计文件
  type1<- type
  if (argv$order!="F"){
    df$type <- factor(df$type,levels=unlist(strsplit(argv$order,",")))
  } else {
    df$type <- factor(df$type, levels = levels(rds@meta.data[,type]))
  }
  
  
  if (is.list(color)) {
    df$color <- color[[paste0(type1,"_colors")]][match(df$type,names(color[[paste0(type1,"_colors")]]))]  # 传颜色参数至df中
    df <- as.data.frame(df)
	
    box_color <- data.frame(color)
    box_color <- rownames_to_column(box_color, var="ident")
    box_color$ident <- factor(box_color$ident,levels=levels(df$type))
    box_color <- box_color[box_color$ident %in% unique(df$type),]

    p1 <- boxplot_function(data=df, x="type", y="ITH_score", group="type", filltype="fill", usecolor=box_color[order(box_color$ident),2], savefig=FALSE)
    p1 <- p1 + stat_compare_means(comparisons = gcompare, label = "p.signif", method = argv$cmethod)

  } else {
    df <- as.data.frame(df)
    p1 <- boxplot_function(data=df, x="type", y="ITH_score", group="type", filltype="fill", savefig=FALSE)
    p1 <- p1 + stat_compare_means(comparisons = gcompare, label = "p.signif", method = argv$cmethod)               

  }
  pdf(paste0(outdir,"/",prefix,"_ITH_score_boxplot.pdf"),width=fig_width)
  print(p1)
  dev.off()
  svg(paste0(outdir,"/",prefix,"_ITH_score_boxplot.svg"),width=fig_width)
  print(p1)
  dev.off()
  png(paste0(outdir,"/",prefix,"_ITH_score_boxplot.png"),width=fig_width, height=8, units='in',res=400)
  print(p1)
  dev.off() 
} else {                                  # 绘制无组间比较的小提琴图
type1 <- type
if (argv$order!="F"){
  df$type <- factor(df$type,levels=unlist(strsplit(argv$order,",")))
} else {
  df$type <- factor(df$type, levels = levels(rds@meta.data[,type]))
}

if (is.list(color)) {
  df <- as.data.frame(df)
  df$color <- color[[paste0(type1,"_colors")]][match(df$type,names(color[[paste0(type1,"_colors")]]))]  # 传颜色参数至df中
  
	box_color <- data.frame(color)
	box_color <- rownames_to_column(box_color, var="ident")
  box_color$ident <- factor(box_color$ident,levels=levels(df$type))
	box_color <- box_color[box_color$ident %in% unique(df$type),]

  p <- boxplot_function(data=df, x="type", y="ITH_score", group="type", filltype="fill", usecolor=box_color[order(box_color$ident),2], rmoutlier=TRUE, fig_save_name=paste0(outdir,"/",prefix,"_ITH_score_boxplot")) 
  
} else {
	df <- as.data.frame(df)
  p <- boxplot_function(data=df, x="type", y="ITH_score", group="type", filltype="fill", rmoutlier=TRUE, fig_save_name=paste0(outdir,"/",prefix,"_ITH_score_boxplot")) 
}

}
}

#run
if (method=="dist"){
  ITH_dist(rds=rds,
           idents=idents,
           type=type,
           ncell.min=ncell.min,
           ncell.max=ncell.max,
           color=color1,
           outdir=outdir,
           prefix=prefix,
		   fig_width=fig_width)
}else if(method=="iqr"){
  ITH_iqr(rds=rds,
           idents=idents,
           type=type,
           ncell.min=ncell.min,
           ncell.max=ncell.max,
           color=color1,
           outdir=outdir,
           prefix=prefix,
		   fig_width=fig_width)
}
#file.copy ("/Public/Script/shouhou/README/tumor-readme/ITH-score_README.docx",outdir)
# file.copy (from=readmePath, to=outdir)
