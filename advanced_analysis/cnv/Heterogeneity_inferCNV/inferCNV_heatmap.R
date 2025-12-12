#!/usr/bin Rscript
# conda activate fsEnv
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
library(cowplot)
library(argparser)
library(Seurat)

argv <- arg_parser('')
argv <- add_argument(argv,"--input", help="the infercnv result path")
argv <- add_argument(argv,"--cluster", help="subcluster list ,split by ,", default='all')
argv <- add_argument(argv,"--sample", help="sample list ,split by ,", default='all')
argv <- add_argument(argv,"--group", help="group list ,split by ,", default='all')
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- add_argument(argv,"--prefix",help="prefix")
argv <- add_argument(argv,"--chr",help="gene location file")
argv <- add_argument(argv,"--kgroups",help="Number of groups in which to break the observations",default = "F")
argv <- add_argument(argv,"--rds",help="rds file",default="F")
argv <- add_argument(argv,"--type1",help="sample,group,celltype",default="subclone")
argv <- add_argument(argv,"--sort1",help="whether sort for type1",default="F")
argv <- add_argument(argv,"--type2",help="sample,group,celltype",default="celltype")
argv <- add_argument(argv,"--sort2",help="whether sort for type2",default="F")
argv <- add_argument(argv,"--rowgap",help="the gap of row type1",default="auto")
argv <- add_argument(argv,"--pheight",help="pdf height",default=8)
argv <- add_argument(argv,"--pwidth",help="pdf width",default=15)
argv <- parse_args(argv)

# argv$input <- "C:/Users/Fangsheng/Desktop/za/02.a.infercnv"
# argv$outdir <- "D:/流程转化/inferCNv/heatmap/test12345"
# argv$prefix <- "test11132"
# argv$chr <- "D:/流程转化/scCancer恶性肿瘤评估RD21082401/test/my_gene_no_sex_chr.txt"
# argv$grouping <- "C:/Users/Fangsheng/Desktop/za/02.a.infercnv/infercnv.observation_groupings.txt"
# argv$rds <- "C:/Users/Fangsheng/Desktop/za/gynaecomastia.diff_PRO.rds"
# argv$sort2 <- "huaian-1,gynaecomastia"
# argv$type1 <- "subclone"
# argv$type2 <- "sample"

#param
input <- argv$input
# obser_input <- paste0(input,"/infercnv.observations.txt")
# refer_input <- paste0(input,"/infercnv.references.txt")
grouping1 <-  paste0(input,"/infercnv.observation_groupings.txt")
outdir <- argv$outdir
prefix <- argv$prefix
chr <- argv$chr
type1 <- argv$type1
type2 <- argv$type2
sort1 <- argv$sort1
sort2 <- argv$sort2
pheight <- as.numeric(argv$pheight)
pwidth <- as.numeric(argv$pwidth)
dir.create(outdir)
if (argv$kgroups != "F"){
  clusternum <- as.numeric(argv$kgroups) 
}

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

#input
obj_1 <- paste0(input,"/run.final.infercnv_obj")
infercnv_obj <- readRDS(obj_1)
expr <- infercnv_obj@expr.data
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc <- unlist(test_loc)
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <-  unlist(normal_loc)
barcode_t <-colnames(expr)[test_loc]
barcode_n <-colnames(expr)[normal_loc ]
print(length(barcode_t))



print("------compressed data------")
x.range="auto"
plot_data = as.matrix(expr)
x.center=mean(plot_data)
quantiles = quantile(plot_data[plot_data != x.center], c(0.01, 0.99))

if ( (length(x.range) == 1) & (x.range[1] == "auto") ) {
  
  # examine distribution of data that's off-center, since much of the center could
  # correspond to a mass of data that has been wiped out during noise reduction
  quantiles = quantile(plot_data[plot_data != x.center], c(0.01, 0.99))
  
  # determine max distance from the center.
  delta = max( abs( c(x.center - quantiles[1],  quantiles[2] - x.center) ) )
  low_threshold = x.center - delta
  high_threshold = x.center + delta
  x.range = c(low_threshold, high_threshold)

  
} else {
  
  # use defined values
  low_threshold = x.range[1]
  high_threshold = x.range[2]
  
  if (low_threshold > x.center | high_threshold < x.center | low_threshold >= high_threshold) {
    stop(paste("Error, problem with relative values of x.range: ", x.range, ", and x.center: ", x.center))
  }
}

plot_data[plot_data < low_threshold] <- low_threshold
plot_data[plot_data > high_threshold] <- high_threshold
expr <- as.data.frame(plot_data)

#sub
print("----sub cell----")
rds_raw <- readRDS(argv$rds)
rds_raw$cluster <- rds_raw@active.ident
rds_raw_ref <- rds_raw
sourceScript("seuratSingleron.v1.R")
sample <- unlist(strsplit(argv$sample,split=","))
cluster <- unlist(strsplit(argv$cluster,split=","))
group <- unlist(strsplit(argv$group,split=","))
rds_raw  <- subsetRDS(rds_raw , sample=sample, cluster=cluster, group=group)

barcode_t_sub <- intersect(colnames(rds_raw),barcode_t)
print("observation")
print(length(barcode_t_sub))
barcode_n_sub <- intersect(colnames(rds_raw),barcode_n)
print("reference")
print(length(barcode_n_sub))
obser <- expr[,barcode_t_sub]
refer <- expr[,barcode_n_sub]
print("----sub cell done!----")
geneFile <- read.table(chr,header = F,sep = "\t",stringsAsFactors = F)

grouping <- read.table(grouping1,header = T,sep="")

##infercnv raw heatmap position for order
meta_obser <- rds_raw@meta.data[colnames(obser),]
meta_obser$order <- c(1:nrow(meta_obser))

meta_refer <- rds_raw_ref@meta.data[colnames(refer),]
meta_refer$order <- c(1:nrow(meta_refer))

##gene check for 
print(setdiff(rownames(obser),rownames(refer)))
print(setdiff(rownames(refer),rownames(obser)))
gn <- rownames(obser)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
obser=obser[intersect(gn,geneFile$V1),]
refer=refer[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
print(obser[1:4,1:4])
print(refer[1:4,1:4])

#color for heatmap
#f2 = colorRamp2(c(min(obser), 1 ,max(obser)), c("#00008B", "white", "#8B0000"))
#f2 = colorRamp2(c(min(obser), 1 ,max(obser)), c("darkblue", "white", "darkred"))
f2 = colorRamp2(c(min(obser), 1 ,max(obser)), c("#0000CD", "white", "#800000"))

##check gene position(for test)
if (0){
  geneFile123 <-geneFile
  colnames(geneFile123) <- c("gene","chr","start","end")
  geneFile1 <- data.frame()
  for (i in unique(geneFile123$chr)){
    genefile <- geneFile123[which(geneFile123$chr==i),]
    genefile <- genefile[order(genefile$start,decreasing = F),]
    geneFile1 <- rbind(geneFile1,genefile)
  }
  geneFile1$order <- c(1:nrow(geneFile1))
  geneFile2 <- geneFile1[intersect(gn,geneFile123$gene),]
}

#add grouping cluster(when infercnv cluster is F,grouping cluster is subclone)
print(nrow(grouping)==nrow(meta_obser))
grouping <- grouping[rownames(meta_obser),]
meta_obser$grouping <- grouping$Dendrogram.Group


if (argv$kgroups != "F"){
  #k-means
  tumor_expr <- as.data.frame(t(obser))
  #k-means cluster
  print(nrow(tumor_expr))
  print(clusternum)
  kmeans.result <- kmeans(tumor_expr, clusternum)
  kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
  meta_obser <- cbind(meta_obser,kmeans_df)
  write.table(kmeans_df,paste0(outdir,"/",prefix,"_kmeans.xls"),sep="\t",row.names=T,col.names=F,quote=F)
}else{
  meta_obser$kmeans_class <- grouping$Dendrogram.Group
}


#type1
if (type1=="celltype"){
  anno_obser <- data.frame(CB=rownames(meta_obser),type=meta_obser$cluster,order=meta_obser$order)
}

if (type1=="sample"){
  anno_obser <- data.frame(CB=rownames(meta_obser),type=meta_obser$sample,order=meta_obser$order)   
}

if (type1=="group"){
  anno_obser <- data.frame(CB=rownames(meta_obser),type=meta_obser$group,order=meta_obser$order)   
}

if (type1=="subclone"){
  anno_obser <- data.frame(CB=rownames(meta_obser),type=paste0("clone",meta_obser$kmeans_class),order=meta_obser$order)  
  
}

if (sort1!="F"){
  print("----sort type1----")
  sort1 <- unlist(strsplit(sort1,split = ","))
  print(sort1)
  anno_obser$type   <- factor(anno_obser$type,levels = sort1)
}else{
  anno_obser$type   <- factor(anno_obser$type,levels = unique(anno_obser$type))
}

#type
if (type2=="celltype"){
  anno_obser$type1 <- meta_obser$cluster
}

if (type2=="sample"){
  anno_obser$type1 <- meta_obser$sample
}

if (type2=="group"){
  anno_obser$type1 <- meta_obser$group
}

if (type2=="subclone"){
  anno_obser$type1 <- paste0("clone",meta_obser$kmeans_class)
}

if(type2=="F"){
  anno_obser$type1 <- 1
}
anno_obser$type1 <- as.character(anno_obser$type1)



if (sort2!="F"){
  print("----sort type2----")
  sort2 <- unlist(strsplit(sort2,split = ","))
  print(sort2)
  anno_obser1 <- data.frame()
  for (i in unique(anno_obser$type)){
    anno_obser_tmp1 <-  anno_obser[which(anno_obser$type==i),]
    sort_i <- intersect(sort2,unique(anno_obser_tmp1$type1)) 
    print(i)
    print(sort_i)
    for (j in sort_i){
      print(j)
      anno_obser_tmp2 <- anno_obser_tmp1[which(anno_obser_tmp1$type1==j),]
      anno_obser_tmp2 <- anno_obser_tmp2[order(anno_obser_tmp2$order,decreasing = F),]
      anno_obser1 <- rbind(anno_obser1,anno_obser_tmp2)
    }
  }
  anno_obser1$type1 <- factor(anno_obser1$type1,levels=sort2)
}else{
  anno_obser1 <- anno_obser
}

#gap set
if (argv$rowgap == "auto"){
  if (pheight <= 3){
    rowgap <- 1
  } else if (pheight <= 8){
    rowgap <- 2
  }else{
    rowgap <- 3
  }
} else {
  rowgap <- as.numeric(argv$rowgap)
}

print("----rowgap----")
print(rowgap)


plot_type1 <- function(tumoranoo,
                        tumor_expr,
                        sub_geneFile,
                        prefix,
                        f2,
                        type1,
                        pheight,
                        pwidth,
                        rowgap,
                        outdir
){ 
tumor_expr$chr <- as.factor(sub_geneFile$V2)
print(table(tumor_expr$chr))
new_cluster <- unique(tumor_expr$chr)
get_group_color_palette <- function () {
  return(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3")))
}
color <- get_group_color_palette()(length(unique(tumor_expr$chr)))
color1 <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B",  "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9",  "#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000",  "#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
top_color <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = color),
                       labels = levels(new_cluster), 
                       labels_gp = gpar(cex = 0.9, col = "black"))) 

tumora <- unique(tumoranoo$type)
tumorb <- color1[1:length(tumora)]
print(unique(tumoranoo$type))
names(tumorb) <- unique(tumoranoo$type)
tumorc <- list(type=tumorb)

tumor_left_anno <- rowAnnotation(type = tumoranoo$type ,col = tumorc,show_annotation_name=F,annotation_name_gp = gpar(fontsize = 20),annotation_legend_param=list( type =list(simple_anno_size=unit(4,"mm"),legend_width=unit(pwidth*0.5, "cm"),legend_height = unit(pheight, "cm"))))
tumor_left_split <- data.frame(type=tumoranoo$type)
n <- t(tumor_expr[,-ncol(tumor_expr)])
n <- n[tumoranoo$CB,]

ht_opt$message = FALSE

print("----plotheatmap with obser by type1----")
ht_tumor = Heatmap(as.matrix(n),
                   col = f2, 
                   cluster_rows = F,
                   cluster_columns = F,
                   show_column_names = F,
                   show_row_names = F,
                   #column_split = new_cluster,
                   column_split = factor(sub_geneFile$V2),
                   row_split = tumor_left_split,
                   border = TRUE,
                   column_gap = unit(0, "mm"),
                   row_gap = unit(rowgap, "mm"),
                   heatmap_legend_param = list(
                     title = "Modified Expression",
                     title_position = "leftcenter-rot",
                     annotation_name_gp = gpar(fontsize = 20), 
                     # at=c(0.4,1.6), #图例范围
                     legend_height = unit(pheight*0.8, "cm") #图例长度
                   ),
                   left_annotation= tumor_left_anno,
                  
                   top_annotation = top_color,
                   row_title = "Observations (Cells)",
                   row_title_side = c("right"),
                   column_title = "Genomic Region",
                   column_title_side = c("bottom"))
# saveRDS(ht_tumor,paste0(outdir,"/",prefix,"heatmap.rds"))
p1 <- paste0(outdir,"/",prefix,"_observations_heatmap_by_",type1,".pdf")
pdf(p1,width =pwidth ,height =pheight,bg="white")
draw(ht_tumor, heatmap_legend_side = "right")
dev.off()

p1 <- paste0(outdir,"/",prefix,"_observations_heatmap_by_",type1,".png")
png(p1,width =pwidth ,height =pheight,bg="white",units='in',res=450)
draw(ht_tumor, heatmap_legend_side = "right")
dev.off()
# pngfile <- sub('pdf$', 'png', p1)
# system(paste0('convert ', p1, ' ', pngfile))
  
}

plot_type12 <- function(tumoranoo,
                        tumor_expr,
                        sub_geneFile,
                        prefix,
                        f2,
                        type1,
                        type2,
                        pheight,
                        pwidth ,
                        rowgap,
                        outdir
                        ){

#tumoranoo <- anno_obser1
#tumor_expr <- obser
tumor_expr$chr <- as.factor(sub_geneFile$V2)
table(tumor_expr$chr)
new_cluster <- unique(tumor_expr$chr)
get_group_color_palette <- function () {
  return(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3")))
}
color <- get_group_color_palette()(length(unique(tumor_expr$chr)))


color1  <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold", "DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4",  "#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B",  "#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080",  "#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49",  "#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD")
color2 <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B",  "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9",  "#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000",  "#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
top_color <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = color),
                       labels = levels(new_cluster), 
                       labels_gp = gpar(cex = 0.9, col = "black"))) 

tumora <- unique(tumoranoo$type)
tumorb <- color1[1:length(tumora)]
names(tumorb) <- unique(tumoranoo$type)
tumord <- color2[1:length(unique(tumoranoo$type1))]
names(tumord) <- unique(tumoranoo$type1)
tumorc <- list(type1=tumorb,type2=tumord)

tumor_left_anno <- rowAnnotation(type1= tumoranoo$type ,type2=tumoranoo$type1,col = tumorc,show_annotation_name=T,annotation_name_gp=gpar(fontsize = 15),annotation_legend_param=list( type1 =list(simple_anno_size=unit(4,"mm"),legend_width=unit(pwidth*0.5, "cm"),legend_height = unit(pheight, "cm")), type2 =list(simple_anno_size=unit(4,"mm"),legend_width=unit(pwidth*0.5, "cm"),legend_height = unit(pheight, "cm"))   ))
tumor_left_split <- data.frame(type1=tumoranoo$type)
n <- t(tumor_expr[,-ncol(tumor_expr)])
n <- n[tumoranoo$CB,]

print("----plotheatmap with obser by type1 and type2----")
ht_opt$message = FALSE

ht_tumor = Heatmap(as.matrix(n),
                   col = f2, 
                   cluster_rows = F,
                   cluster_columns = F,
                   show_column_names = F,
                   show_row_names = F,
                   #column_split = new_cluster,
                   row_split = tumor_left_split,
                   border = TRUE,
                   column_gap = unit(0, "mm"),
                   row_gap = unit(rowgap, "mm"),
                   column_split = factor(sub_geneFile$V2),
                   heatmap_legend_param = list(
                     title = "Modified Expression",
                     annotation_name_gp = gpar(fontsize = 20),
                     title_position = "leftcenter-rot", # ͼ������λ��
                     # at=c(0.4,1.6), #ͼ����Χ
                     legend_height = unit(pheight*0.8, "cm") #ͼ������
                   ),
                   left_annotation=tumor_left_anno ,
                   top_annotation = top_color,
                   row_title = "Observations (Cells)",
                   row_title_side = c("right"),
                   column_title = "Genomic Region",
                   column_title_side = c("bottom"))
p1 <- paste0(outdir,"/",prefix,"_observations_heatmap_by_",type1,"_",type2,".pdf")
pdf(p1,width =30 ,height =pheight,bg="white")
draw(ht_tumor, heatmap_legend_side = "right")
dev.off()

p1 <- paste0(outdir,"/",prefix,"_observations_heatmap_by_",type1,"_",type2,".png")
png(p1,width =30 ,height =pheight,bg="white",units = 'in',res=450)
draw(ht_tumor, heatmap_legend_side = "right")
dev.off()
# pngfile <- sub('pdf$', 'png', p1)
# system(paste0('convert ', p1, ' ', pngfile))

}

anno_refer <- data.frame(CB=rownames(meta_refer),type=meta_refer$cluster)


plot_ref <- function( normalanoo,
                     normal_expr,
                     sub_geneFile,
                     f2,
                     prefix,
                     outdir){
 # normal_expr <- refer
#  normalanoo <-anno_refer
  
  normal_expr$chr <- as.factor(sub_geneFile$V2)
  table(normal_expr$chr)
  new_cluster <- unique(normal_expr$chr)
  get_group_color_palette <- function () {
    return(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3")))
  }
  color <- get_group_color_palette()(length(unique(normal_expr$chr)))
  
  
   color1  <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold", "DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4",  "#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B",  "#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080",  "#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49",  "#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD")
color2 <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B",  "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9",  "#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000",  "#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
  top_color <- HeatmapAnnotation(
    cluster = anno_block(gp = gpar(fill = color),
                         labels = levels(new_cluster), 
                         labels_gp = gpar(cex = 0.9, col = "black"))) 
  
  normala <-unique(normalanoo$type)
  normalb <- color1[1:length(normala)]
  names(normalb) <- unique(normalanoo$type)
  normalc <- list(type=normalb)
  
  normal_left_anno <- rowAnnotation(type = normalanoo$type,col =normalc,show_annotation_name=F,annotation_name_gp=gpar(fontsize = 10))
  normal_left_split <- data.frame(type = normalanoo$type)
  m <- t(normal_expr[,-ncol(normal_expr)])
  m <- m[normalanoo$CB,]
  print("----plotheatmap with refer----")
  ht_opt$message = FALSE
  ht_normal = Heatmap(as.matrix(m),
                      col = f2, 
                      cluster_rows = F,
                      cluster_columns = F,
                      show_column_names = F,
                      show_row_names = F,
                      #column_split = new_cluster,
                      row_split = normal_left_split,
                      border = TRUE,
                      column_gap = unit(0, "mm"),
                      row_gap = unit(rowgap, "mm"),
                      column_split = factor(sub_geneFile$V2),
                      heatmap_legend_param = list(
                        title = "Modified Expression",
                        title_position = "leftcenter-rot", # ͼ������λ��
                        # at=c(0.4,1.6), #ͼ����Χ
                        legend_height = unit(pheight*0.8, "cm") #ͼ������
                      ),
                      left_annotation=normal_left_anno ,
                      top_annotation = top_color,
                      row_title = "References (Cells)",
                      row_title_side = c("right"),
                      column_title = "Genomic Region",
                      column_title_side = c("bottom"))
  p1 <- paste0(outdir,"/",prefix,"_references_heatmap.pdf")
  pdf(p1,width =15 ,height = 8,bg="white")
  draw(ht_normal, heatmap_legend_side = "right")
  dev.off()

  p1 <- paste0(outdir,"/",prefix,"_references_heatmap.png")
  png(p1,width =15 ,height = 8,bg="white",units = 'in',res = 450)
  draw(ht_normal, heatmap_legend_side = "right")
  dev.off()
  # pngfile <- sub('pdf$', 'png', p1)
  # system(paste0('convert ', p1, ' ', pngfile))
  
}

if (type2!="F"){
  plot_type12(tumoranoo=anno_obser1,tumor_expr=obser,sub_geneFile=sub_geneFile,prefix=prefix,f2=f2,type1=type1,type2=type2,pheight=pheight,pwidth =pwidth ,rowgap=rowgap,outdir=outdir)
}else{
  plot_type1(tumoranoo=anno_obser1,tumor_expr=obser,sub_geneFile=sub_geneFile,prefix=prefix,f2=f2,type1=type1,pheight=pheight,pwidth =pwidth ,rowgap=rowgap,outdir=outdir)
}

plot_ref( normalanoo=anno_refer,normal_expr=refer,sub_geneFile=sub_geneFile,prefix=prefix,f2=f2,outdir=outdir)
