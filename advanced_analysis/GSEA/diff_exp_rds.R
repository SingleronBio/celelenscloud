#!/usr/bin/env Rscript
suppressMessages({
library(ggplot2)
library(reshape2)
library(argparser)
library(Cairo)
library(Seurat)
library(corrplot)
library(dplyr)
library(grid)
library(cowplot)
library(dplyr)
library(viridis)
library(grid)
library (stringr)
})

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="the rds file")
argv <- add_argument(argv,"--subcluster", help="subcluster list ,split by ,", default="all")
argv <- add_argument(argv,"--sample", help="sample list ,split by ,", default="all")
argv <- add_argument(argv,"--gname", help="sample list ,split by ,", default="F")
argv <- add_argument(argv,"--DefaultAssay", help="sub cluster seurat[integrated|]",default="RNA")
argv <- add_argument(argv,"--prefix", help="group name prefix")
argv <- add_argument(argv,"--diffcluster", help="compare cluster  name:1:4vs2,3:4vs4,or T:NKvsB,NKvsB; split by ,",default="F")
argv <- add_argument(argv,"--average",help = "Average gene expression method, mean or Seurat", default = "Seurat")
argv <- add_argument(argv,"--group", help="compare group name:G1vsG2; split by ,",default="F")
argv <- add_argument(argv,"--logfcthreshold", help="Limit testing to genes which show,on average, at least X-fold difference (log-scale) between the two groups of cells.",default=0.25)
argv <- add_argument(argv,"--minpct", help="only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations.",default=0.1)
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- parse_args(argv)

print(argv$rds)
print(argv$subcluster)
do_Assay<-as.character(argv$DefaultAssay)
rdsfile<-argv$rds
group <- unlist(strsplit(argv$group,split=","))
diffcluster <- unlist(strsplit(argv$diffcluster,split=","))
samplelist <- unlist(strsplit(argv$sample,split=","))
gnamelist <- unlist(strsplit(argv$gname,split=","))
outdir <- argv$outdir
dir.create(outdir)
compare <- argv$prefix
logfc <- as.numeric(argv$logfcthreshold)
minpct_u <- as.numeric(argv$minpct)
#stat
col1 <- colorRampPalette(c("#7F0000","red","red","#FF7F00","#FF7F00","yellow","yellow","cyan", "#007FFF", "blue", "#00007F"))
corrcol <- colorRampPalette(c("red","orange","blue","white","white"))
clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
col2<-colorRampPalette(clustcol)
###Mean Expression Function###
MeanExp <- function(object, genes){
	exp <- object@assays$RNA@data
	exp<- data.frame(exp[genes,])
  exp <- as.data.frame(t(exp))
	exp$Term <- object@active.ident
	
	table <- aggregate(exp[,genes],by = list(exp$Term), FUN = mean)
  rn <- table[,1]
  table <- as.data.frame(table[,-1])
  colnames(table) <- genes
  rownames(table) <- rn
  table <- t(table)
  return(table)
}


######
print("loading RDS!!!")
PRO<-readRDS(rdsfile)
subc<-levels(PRO)

print("cluster ID now is:")
print(subc)
substr<-c()
if (grepl(',',argv$subcluster)){
        substr<-unlist(strsplit(argv$subcluster,split=","))
}else{
        substr<-c(argv$subcluster)
}
print(substr)
if (subc[1] == 0){
        print(subc[1])
        substr<-as.character(as.numeric(substr)-1)
}
print(substr)
DefaultAssay(PRO) <- "RNA"
LL<-cbind(data.frame(names(PRO@active.ident)),data.frame(PRO@meta.data$sample))
print(head(LL))
cell_use<-c()
print(samplelist)
print("sub cluster:")
if ( substr=='all' ){
  print('all idents used') 
}else{
  PRO  <- subset(PRO, idents = substr)
}

if ( samplelist=='all' ){
  print('all idents used') 
}else{
  PRO  <- subset(PRO, subset = sample %in% samplelist)
}
PRO$sample <- as.character(PRO$sample)


all.genes <- rownames(PRO)
PRO <- ScaleData(PRO, features = all.genes)
subc<-levels(PRO)
print(subc)
if (grepl('vs',argv$diffcluster)){
	for (cm in diffcluster){
		if (grepl (':',cm)) {
			id <- str_replace_all (cm,":",".")
		}else {
			id <- cm
		}
		compareC <- unlist(strsplit(cm,split="vs"))
		if(grepl(':',compareC[1])){
			c1<-unlist(strsplit(compareC[1],split=":"))
		}else{
			c1<-compareC[1]
		}
		if(grepl(':',compareC[2])){
			c2<-unlist(strsplit(compareC[2],split=":"))
		}else{
			c2<-compareC[2]
		}
    
		diff <- FindMarkers(PRO, ident.1 = c1, ident.2 = c2,min.pct=minpct_u,logfc.threshold = logfc)
		write.table(data.frame(gene_id=rownames(diff),diff),file=paste(outdir,'/',compare,'_',id,'.diffexpressed.xls',sep=''),sep='\t',quote=F,row.names=F)
		LL<-cbind(data.frame(rownames(diff)),data.frame(diff))
		markergenetop1<-LL %>% arrange(avg_log2FC) %>%  head(20)
		markergenetop2<-LL %>% arrange(desc(avg_log2FC)) %>%  head(20)
		markergene<-c(as.character(markergenetop2$rownames.diff.),as.character(markergenetop1$rownames.diff.))
    markergene <- unique(markergene)
		print(markergene)
		#markergene<-grep(pattern = "ENS", x = markergene,invert=TRUE,value = TRUE)
    PRO_sub <- subset(PRO, idents = c(c1,c2))
    levels(PRO_sub) <- c(c1,c2)
		PRO_sub<-subset(PRO,downsample=300)
    
		all.genes <- rownames(PRO_sub)
		PRO_sub <- ScaleData(PRO_sub, features = all.genes)
		PRO_sub$celltype <- Idents(PRO_sub)
		Idents(PRO_sub) <- "celltype"
		heatname=paste(cm,' genes heatmap',sep='')
				#markergene<-c(as.character(markergenetop2$rownames.diff.)[1:20],as.character(markergenetop1$rownames.diff.)[1:20])
		#markergene<-grep(pattern = "ENS", x = markergene,invert=TRUE,value = TRUE)
			}
 
  save_file=paste(outdir,'/',compare,'.diff_PRO.rds',sep='')
  PRO <- DietSeurat(PRO, counts = TRUE, data = TRUE, scale.data = FALSE,dimreducs = c('pca','tsne','umap'))
  #saveRDS(PRO,file=save_file)
  if (argv$average == "Seurat"){
  cluster.averages1 <- AverageExpression(object = PRO)
  write.table(cluster.averages1[["RNA"]],file=paste(outdir,'/',compare,'.cluster_averages.xls',sep=''),sep='\t',quote=F,row.names=T,col.names = NA)
  }else if (argv$average == "mean"){
  table <- MeanExp(PRO, genes = rownames(PRO))
  write.table(table,file=paste(outdir,'/',compare,'.cluster_mean.xls',sep=''),sep='\t',quote=F,row.names=T,col.names = NA)
  }
  #dir=paste0(outdir,'/',compare,'/')
  #dir.create(dir)
  #cmd1 <- paste0('mv ',outdir,'/','*DOHeatmapplot* ', dir)
  #system(cmd1)
  #cmd3 <- paste0('mv ',outdir,'/','*diffexpressed.xls ', dir)
  #system(cmd3)
}

if (gnamelist != "F"){
  if ( !('group' %in% colnames(PRO@meta.data))){
  PRO$group = ''
  }else{PRO$group <- as.character(PRO$group)}
  for (i in 1:length(samplelist)){
  if (samplelist[i] %in% unique(PRO$sample)){PRO$group[PRO$sample == samplelist[i]] <- gnamelist[i]}
  
  }

}



PRO$celltype.sample <- PRO$group
PRO$celltype <- Idents(PRO)
print("change ID")
Idents(PRO) <- "celltype.sample"
print("change ID done")
subc<-levels(PRO)
print(subc)
if (grepl('vs',argv$group)){
for (cm in group){
	print(cm)
	if (grepl (':',cm)) {
		id <- str_replace_all (cm,":",".")
	}else {
		id <- cm
	}
	comparegroup <- unlist(strsplit(cm,split="vs"))
	if(grepl(':',comparegroup[1])){
		c1<-unlist(strsplit(comparegroup[1],split=":"))
	}else{
		c1<-comparegroup[1]
	}
	if(grepl(':',comparegroup[2])){
		c2<-unlist(strsplit(comparegroup[2],split=":"))
	}else{
		c2<-comparegroup[2]
	}
	cel.1 = rownames(PRO@meta.data)[ PRO$celltype.sample %in% c1 ]
	cel.2 = rownames(PRO@meta.data)[ PRO$celltype.sample %in% c2 ]
	if ( length(cel.1)<3 | length(cel.2)<3 ){ print( paste0(cm,' not enough cells') ) ; next }
	c1 = c1[c1 %in% Idents(PRO)] ; print(paste0(c1,':',length(cel.1)))
	c2 = c2[c2 %in% Idents(PRO)] ; print(paste0(c2,':',length(cel.2)))
	diff <- FindMarkers(PRO,  ident.1= c1, ident.2 = c2,min.pct=minpct_u,logfc.threshold = logfc)
#	write.table(diff,file=paste(outdir,'/',compare,'_',cm,'.diffexpressed.xls',sep=''),sep='\t',quote=F,row.names=T,col.names = NA)
	write.table(data.frame(gene_id=rownames(diff),diff),file=paste(outdir,'/',compare,'_',id,'.diffexpressed.xls',sep=''),sep='\t',quote=F,row.names=F)
	LL<-cbind(data.frame(rownames(diff)),data.frame(diff))
	markergenetop1<-LL %>% arrange(avg_log2FC) %>%  head(20)
	markergenetop2<-LL %>% arrange(desc(avg_log2FC)) %>%  head(20)
	markergene<-c(as.character(markergenetop2$rownames.diff.),as.character(markergenetop1$rownames.diff.))
  markergene <- unique(markergene)
	print(markergene)
	#markergene<-grep(pattern = "ENS", x = markergene,invert=TRUE,value = TRUE)
	#PRO_sub<-SubsetData(PRO,max.cells.per.ident=300)
  Idents(PRO) <- 'sample'
  PRO_sub <- subset(PRO, downsample = 300)
  Idents(PRO) <- 'celltype.sample'
  Idents(PRO_sub) <- 'celltype.sample'
  PRO_sub <- subset(PRO_sub,subset = group %in% c(c1,c2))
	all.genes <- rownames(PRO_sub)
	PRO_sub <- ScaleData(PRO_sub, features = all.genes)
	heatname=paste(cm,' genes heatmap',sep='')
  levels(PRO_sub) <- c(c1,c2)
  splevel <- c(unique(PRO_sub$sample[PRO_sub$group == c1]),unique(PRO_sub$sample[PRO_sub$group == c2]))
  PRO_sub$sample <- factor(PRO_sub$sample,levels = splevel)
 
	#markergene<-c(as.character(markergenetop2$rownames.diff.)[1:20],as.character(markergenetop1$rownames.diff.)[1:20])
  #markergene <- unique(markergene)
	#markergene<-grep(pattern = "ENS", x = markergene,invert=TRUE,value = TRUE)
}

save_file=paste(outdir,'/',compare,'.diff_PRO.rds',sep='')
PRO <- DietSeurat(PRO, counts = TRUE, data = TRUE, scale.data = FALSE,dimreducs = c('pca','tsne','umap'))
saveRDS(PRO,file=save_file)
if (argv$average == "Seurat"){
  cluster.averages1 <- AverageExpression(object = PRO)
  write.table(cluster.averages1[["RNA"]],file=paste(outdir,'/',compare,'.cluster_averages.xls',sep=''),sep='\t',quote=F,row.names=T,col.names = NA)
  }else if (argv$average == "mean"){
  table <- MeanExp(PRO, genes = rownames(PRO))
  write.table(table,file=paste(outdir,'/',compare,'.cluster_mean.xls',sep=''),sep='\t',quote=F,row.names=T,col.names = NA)
  }

#dir=paste0(outdir,'/',compare,'/')
#dir.create(dir)
#cmd1 <- paste0('mv ',outdir,'/','*DOHeatmapplot_group* ', dir)
#system(cmd1)
#cmd2 <- paste0('mv ',outdir,'/','*DOHeatmapplot_sample* ', dir)
#system(cmd2)
#cmd3 <- paste0('mv ',outdir,'/','*diffexpressed.xls ', dir)
#system(cmd3)
}
