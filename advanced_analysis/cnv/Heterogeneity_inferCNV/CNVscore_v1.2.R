#! /usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(scales)
library(argparser)
library(ggpubr)
library (dplyr)
library(infercnv)
library(BiocParallel)
library(scales)

argv <- arg_parser('')
argv <- add_argument(argv,"--inputdir", help="the infercnv directory including run.final.infercnv_obj")
argv <- add_argument(argv,"--objlist", help="column1 sample,column2 infercnv result obj,table split",default="F")
argv <- add_argument(argv,"--split_ref", help = "whether to split sample for ref celltype ,",default = "T")
argv <- add_argument(argv,"--rds", help="the input rds")
argv <- add_argument(argv,"--sample", help="the sample name for running CNVscore",default="all")
argv <- add_argument(argv,"--including_clust", help="the celltype  for ploting CNVscore",default="all")
argv <- add_argument(argv,"--ref_ct", help="the reference name for running inferCNVscore",default="F")
argv <- add_argument(argv,"--group", help="the group name for running inferCNVscore",default="all")
argv <- add_argument(argv,"--ident",help="the ident for x-alix,choose from sample,group,cluster,HCL",default="HCL")
argv <- add_argument(argv,"--order",help="the x-ais order",default="F")
argv <- add_argument(argv,"--color",help="the color splitby , ",default="default")
argv <- add_argument(argv,"--vlnplot", help=" whether to plot volinplot", default= "T")
argv <- add_argument(argv,"--boxplot", help="whether to plot boxplot ", default= "T")
argv <- add_argument(argv,"--comparison",help="if add p-vaule in the plot, provide list for comparison",default="F")
argv <- add_argument(argv,"--splitby",help="use sample/group split ident to plot,only works for ident=cluster",default="F")
argv <- add_argument(argv,"--save", help="save rds with score",default=FALSE)
argv <- add_argument(argv,"--prefix", help="file name prefix")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- parse_args(argv)


#source
#-------------------------------------------------------------------------------------------------------------------------
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
scriptdir = dirnameScript()
sourceScript('seurat_function.R')
sourceScript('boxplot_function.r')
sourceScript('plot_Vln.R')
sourceScript("seuratSingleron.v1.R")

#--------------------------------------------------------------------------------------------------------------------------
CNVscore <- function(input) {
    data <- as.matrix(input)
    cnv.value <- data.frame(CNVscore= colMeans(abs(data-1)), CNVscore_SS = apply(abs(data-1), 2, sum), CNVscore_SD = apply(abs((data-1)), 2, sd)) # gene expression of cell was re-standardized and values were limited as -1 to 1
    return(cnv.value)
}
#--------------------------------------------------------------------------------------------------------------------------
split_data.matrix <- function(matrix, chunk.size=1000) {
    ncols <- dim(matrix)[2]
    nchunks <- (ncols-1) %/% chunk.size + 1
    
    split.data <- list()
    min <- 1
    for (i in seq_len(nchunks)) {
        if (i == nchunks-1) {  #make last two chunks of equal size
            left <- ncols-(i-1)*chunk.size
            max <- min+round(left/2)-1
        } else {
            max <- min(i*chunk.size, ncols)
        }
        split.data[[i]] <- matrix[,min:max]
        min <- max+1    #for next chunk
    }
    return(split.data)
}
#--------------------------------------------------------------------------------------------------------------------------
runcnvscore <- function(cnvmatrix,rds,chunk.size=1000,BPPARAM=NULL,ncores=1){
    cnvmatrixlist <- split_data.matrix(cnvmatrix,chunk.size=chunk.size)
    if (is.null(BPPARAM)) {
            BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
        }
    meta.list <- BiocParallel::bplapply(
        X = cnvmatrixlist, 
        BPPARAM =  BPPARAM,
        FUN = function(x) {
            cnvscore <- CNVscore(x)
            return(list(cnvscore=cnvscore))
        }
    )
    meta.merge <- lapply(meta.list,function(x) cbind(x[["cnvscore"]]))
    meta.merge <- Reduce(rbind, meta.merge)
    intercb <- intersect(rownames(meta.merge),colnames(rds))
    rds <- rds[,intercb]
    meta.merge <- meta.merge[intercb,]
    rds <- Seurat::AddMetaData(rds, as.data.frame(meta.merge))
    return(rds)
}
#merge
#-------------------------------------------------------------------------------------------------------------------------
mergeobj <- function(objdata,ref="same"){
    print("-----merge cnv matrix-----")
    i=1
    sample_1 <- objdata$sample[i]
    print(sample_1)
    obj_1 <- objdata$obj[i]
    infercnv_obj <- readRDS(obj_1)
    expr <- infercnv_obj@expr.data
    test_loc <- infercnv_obj@observation_grouped_cell_indices
    test_loc <- unlist(test_loc)
    normal_loc <- infercnv_obj@reference_grouped_cell_indices
    normal_loc <-  unlist(normal_loc)
    barcode_t <-colnames(expr)[test_loc]
    barcode_n <-colnames(expr)[normal_loc ]
    print(length(barcode_t))
    print("----------")

    if (nrow(objdata) >1){
        if (ref != "same"){
            for (i in 2:nrow(objdata)){
                sample_tmp <- objdata$sample[i]
                print(sample_tmp)
                obj_tmp <- objdata$obj[i]
                infercnv_obj <- readRDS(obj_tmp)
                expr_tmp <- infercnv_obj@expr.data
                test_loc <- infercnv_obj@observation_grouped_cell_indices
                test_loc <- unlist(test_loc)
                normal_loc <- infercnv_obj@reference_grouped_cell_indices
                normal_loc <-  unlist(normal_loc)
                barcode_t_tmp <-colnames(expr_tmp)[test_loc]
                barcode_n_tmp <-colnames(expr_tmp)[normal_loc ]
                print(length(barcode_t_tmp))
                barcode_t <- c(barcode_t,barcode_t_tmp)
                barcode_n <- c(barcode_n,barcode_n_tmp)
                #check gene and merge cnv matrix
                gene_tmp <- intersect(rownames(expr),rownames(expr_tmp))
                expr_tmp <- expr_tmp[gene_tmp,]
                expr <- expr[gene_tmp,]
                expr <- cbind(expr,expr_tmp)
                print("----------")
                }
        } else {
            for (i in 2:nrow(objdata)){
                sample_tmp <- objdata$sample[i]
                print(sample_tmp)
                obj_tmp <- objdata$obj[i]
                infercnv_obj <- readRDS(obj_tmp)
                expr_tmp <- infercnv_obj@expr.data
                test_loc <- infercnv_obj@observation_grouped_cell_indices
                test_loc <- unlist(test_loc)
                normal_loc <- infercnv_obj@reference_grouped_cell_indices
                normal_loc <-  unlist(normal_loc)
                barcode_t_tmp <-colnames(expr_tmp)[test_loc]
                print(length(barcode_t_tmp))
                barcode_t <- c(barcode_t,barcode_t_tmp)
                #check gene and merge cnv matrix
                gene_tmp <- intersect(rownames(expr),rownames(expr_tmp))
                expr_tmp <- expr_tmp[gene_tmp,barcode_t_tmp]
                expr <- expr[gene_tmp,]
                expr <- cbind(expr,expr_tmp)
                print("----------")
                }
        }
    }
    return(list(expr=expr,barcode_t=barcode_t,barcode_n=barcode_n))
}

cnvmatrix_scale <- function(expr,x.range="auto"){
    print("------compressed data------")
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
    print("------compressed data done!------")
    return(plot_data)
}

#vlnplot
#-------------------------------------------------------------------------------------------------------------------------
vlnplot_function <- function(data,x,y, group,filltype='all', 
                            title_plot=NULL,title_x=NULL,title_y=y,title_legend=NULL,x_angle=60,usecolor=NULL){
if (filltype == "fill"){
        fillgroup=data[,group]
        colorgroup=NULL
        #alpha=1
        alpha=0.7
    }else if (filltype == "color"){
        fillgroup=NULL
        colorgroup=data[,group]
        alpha=1
    }else if(filltype == "all"){
        fillgroup=data[,group]
        colorgroup=data[,group]
        alpha=0.7
    }

p1=ggplot(data,aes(x=data[,x],y=data[,y],fill=fillgroup,color=colorgroup))+
            geom_violin(alpha = alpha,position = position_dodge(0.85)) +
            geom_boxplot(lwd=1,width=0.1,outlier.colour = NA,outlier.shape = NA,position = position_dodge(0.85),alpha=alpha)+
            labs(x=title_x,y=title_y,title=title_plot)+
            theme_classic()+
            scale_fill_manual(values= usecolor)+
            #scale_color_manual(values= usecolor)+
            theme (plot.title = element_text(hjust = 0.5,face="bold",size=18),
                    axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
                    axis.line.x=element_line(color="black",size=0.5),axis.line.y=element_line(color="black",size=0.5),
                    axis.text.x = element_text(angle = x_angle, vjust = 1, hjust = 1,size=13,colour ="black"),
                    axis.text.y = element_text(size=13,colour="black"),
                    #panel.border=element_rect(fill=NA,color="black",size=0.6),
                    plot.margin = margin(1,1,1,1,'cm') )       
if (group!=x){
        p1<-p1 + guides(fill = guide_legend(title=title_legend),color = guide_legend(title=title_legend))+
                theme(legend.key.size=unit(2,"line"),legend.text=element_text(size=18))
    }else{ 
        p1<-p1+theme(legend.position='NA')
    }
print("*** violinplot done ***")
return(p1)
}

#analysis width and height
#-----------------------------------------------------------------------------------------------
get_params <- function(data,splitby,ident){
    width <- 8
    height <- 8
    if (splitby != "F"){
      if (length(unique(data[,"ident"]))<=4){
        if (length(unique(data[,splitby]))<=4){
         width=(0.35*(length(unique(data[,splitby]))))*(length(unique(data[,"ident"])))+2
        }else if (length(unique(data[,splitby]))<=15){
           width=(0.3*(length(unique(data[,splitby]))))*(length(unique(data[,"ident"])))
        }else {
             width=(0.25*(length(unique(data[,splitby]))))*(length(unique(data[,"ident"])))
        }        
        }else if(length(unique(data[,"ident"]))<=15){
            if (length(unique(data[,splitby]))<=3){
                width=(0.35*(length(unique(data[,splitby]))))*(length(unique(data[,"ident"])))+2
            }else if (length(unique(data[,splitby]))<=15){
                width=(0.3*(length(unique(data[,splitby]))))*(length(unique(data[,"ident"])))
            }else {
                width=(0.25*(length(unique(data[,splitby]))))*(length(unique(data[,"ident"])))
            }
        }else{
            width=(0.25*(length(unique(data[,splitby]))))*(length(unique(data[,"ident"]))) 
        }
        lengroup<-nchar(as.character(unique(data[,splitby])))
        if (max(lengroup)>7){
            width<-width+(max(lengroup)-7)*0.1 }
    }else{
        width=0.23*(length(unique(data[,"ident"])))+4
    #横坐标x名称长短与图片加宽情况：
        lenx<-nchar(as.character(unique(data[,"ident"])))
        if (max(lenx)>15){
        width=((max(lenx)-10)*0.05)+width
    }}
    #横坐标各x名称长短与图片加长情况：
    #原 +5.5
    lenx<-nchar(as.character(unique(data[,"ident"])))
    if (max(lenx)>10){
        height=(max(lenx)-10)*0.1+5.5 }
        else{
            height=5.5
        }
    return(list(width,height))
}
        
#plot
#--------------------------------------------------------------------------------------------------------------------------
pdf2png <- function(pdf){
    png <- gsub('pdf$', 'png', pdf)
    system(paste('convert -density 450 -background white', pdf, png, sep=' '))
}

#--------------------------------------------------------------------------------------------------------------------------
plotFunc <- function(pic, pdf, width=8, height=8, png=FALSE, convert=FALSE){
    pdf(pdf, width=width, height=height, bg='white')
    print(pic)
    dev.off()

    svg <- gsub('pdf$', 'svg', pdf)
    ggsave(file=svg, plot=pic, width=width, height=height, bg='white')

    if(png){
        png <- gsub('pdf$', 'png', pdf)
        png(png, width=width, height=height, units="in", res=400)
        print(pic)
        dev.off()
    }
    if(convert){
        pdf2png(pdf)
    }
}

#--------------------------------------------------------------------------------------------------------------------------
rds <- readRDS(argv$rds)
input <-argv$inputdir
sample <- unlist(strsplit(argv$sample,split=","))
cluster <- unlist(strsplit(argv$including_clust,split=","))
reference <- unlist(strsplit(argv$ref_ct,split=","))
group <- unlist(strsplit(argv$group,split=","))
boxplot <- argv$boxplot
vlnplot<- argv$vlnplot
prefix <- argv$prefix
ident <- argv$ident
comparison <- unlist(strsplit(argv$comparison,split=","))
splitby <- argv$splitby
saveRDS <- argv$save
outdir <- argv$outdir
dir.create(outdir)
outdir2 <- paste0(outdir,"/",prefix)
dir.create(outdir2)


#read obj
#--------------------------------------------------------------------------------------------------------------------------
if (argv$objlist !="F"){
    objdata <- read.table(argv$objlist,sep="\t",header=F) 
    colnames(objdata) <- c("sample","input")
    name  <- objdata$sample
    if (argv$split_ref != "T"){
        refsame <- "same"
    }else{
        refsame <- "diff"  
    }
    objdata$obj <- paste0(objdata$input,"/run.final.infercnv_obj")
    objdata$hcl <- paste0(objdata$input,"/infercnv.HCL.barcode.txt")
    objdata$grouping <- paste0(objdata$input,"/infercnv.observation_groupings.txt")
    print(name)
    print(objdata)
    expmlist <- mergeobj(objdata=objdata,ref=refsame)
    expr <- expmlist$expr
    barcode_t <- expmlist$barcode_t
    barcode_n <- expmlist$barcode_n
}else{
    infercnv_obj <- readRDS(paste0(input,"/run.final.infercnv_obj"))
    expr <- infercnv_obj@expr.data
    rds@meta.data$cluster<-rds@active.ident
    HCL <- read.table(paste0(input,"/","infercnv.HCL.withanno.barcode.txt"),sep="\t",row.names=1)
    colnames(HCL)<-"subclone"
    rds <- SetIdent(rds,cells=rownames(HCL),value=HCL$subclone)
    rds@meta.data$subclone<-rds@active.ident
    Idents(rds)<-rds@meta.data$cluster
}

#run CNV score
#------------------------------------------------------------------------------------------------------------------------
rds <- runcnvscore(cnvmatrix=expr,rds=rds,chunk.size=1000)
rds@meta.data$cluster<-rds@active.ident
write.table(rds@meta.data,paste0(outdir,"/",prefix,"_cnv_score.xls"),sep="\t",row.names=T,col.names=NA)
if (saveRDS){
    saveRDS(rds,paste0(outdir,"/",prefix,"_all.cnv_score.rds"))
}

#subset cnv score
#--------------------------------------------------------------------------------------------------------------------------
PRO<- subsetRDS(rds,sample=sample, cluster=cluster,group=group)
name1 <- c("sample","group","subclone","cluster","CNVscore","CNVscore_SS","CNVscore_SD")
name2 <- intersect(name1,colnames(PRO@meta.data))
cnv<-PRO@meta.data[,name2]
cnv[,]<-lapply(cnv[,],as.character)
cnv$CNVscore<-as.numeric(cnv$CNVscore)
write.table(cnv,paste(outdir2,"/",prefix,"_cnv_score.xls",sep=""),row.names=T,col.names=NA,sep="\t")


#如果某个分组里只有一个细胞，删掉，不然小提琴图会报错
#---------------------------------------
filter <- function(data,splitby,ident){
    tab <- table(data[,ident],data[,splitby])
    del <- which(tab==1,arr.ind=TRUE)
    if(nrow(del != 0)){
        for (i in 1:nrow(del)){
                spl <- colnames(tab)[del[i,2]]
                ide <- rownames(del)[i]
                row_del <- which(data[,ident]==ide&data[,splitby]==spl)
                if (length(row_del) > 0){
                        data <- data[-row_del,]
                }
        }
    print("splitby plot have filtered one cell for each group")
  }
  return(data)
}

#add p value
if (class(PRO$group) == "factor" & class(PRO$sample) == "factor"){
    cnv$group <- factor(cnv$group,levels=levels(PRO$group))
    cnv$sample <- factor(cnv$sample,levels=levels(PRO$sample))
}else{
    cnv$group <- factor(cnv$group,levels=unique(rds$group))
    cnv$sample <- factor(cnv$sample,levels=unique(rds$sample))
}
if (ident == "group"){
    if (reference == "F"){
        cnv$ident <- cnv$group
    }else{
        cnv$ident <- cnv$group
        cnv$ident <- as.character(cnv$ident)
        cnv[which(cnv$cluster %in% reference),]$ident <- cnv[which(cnv$cluster %in% reference),]$cluster
    }
}else if (ident =="sample"){
    if (reference == "F"){
        cnv$ident <- cnv$sample
    }else{
        cnv$ident <- cnv$sample
        cnv$ident <- as.character(cnv$ident)
        cnv[which(cnv$cluster %in% reference),]$ident <- cnv[which(cnv$cluster %in% reference),]$cluster
    }
}else if (ident =="cluster"){
    cnv$ident <- cnv$cluster
}else{
    cnv$ident <- cnv$subclone
}
if(argv$order != "F"){
    order <- as.character(unlist(strsplit(argv$order,split=",")))
    cnv$ident <- factor(cnv$ident,levels=order)
}else{
    if(cluster == "all"){
    cnv$ident <- factor(cnv$ident,levels=unique(cnv$ident))
    } else if (cluster != "all"& ident == "cluster"){
    cnv$ident <- factor(cnv$ident,levels=cluster)
    }
}
print("now data use for plot:")
print(head(cnv))

num <- length(unique(cnv$ident))
if (num > 1){
    p_value <- compare_means(CNVscore~ident,data=cnv,method="wilcox.test")
    write.table(p_value,paste(outdir2,"/",prefix,"_p_value.xls",sep=""),row.names=F,sep="\t",quote=F)
}
comparegroup <- list()
for(i in 1:length(comparison)){
comparegroup[i] <- strsplit(comparison[i],split="vs")
}

#color
#-------------------------------------------------------------------------------------------------------------------
if (argv$color=="default"){
  if (reference == "F"){
    if (paste0(ident, "_colors") %in% colnames(PRO@meta.data)) {
        color1 <- getColors(PRO, ident)
        color1[[1]] <- color1[[1]][!is.na(color1[[1]])]
        clustcol <- color1[[1]]
        num <- length(unique(cnv$ident))
        color_num <- length(unique(clustcol))
            if (color_num != num){
                color_old <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold", "DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4",  "#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B",  "#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080",  "#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49",  "#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD",  "#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B",  "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9",  "#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000",  "#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
                color_new <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
                clustcol <- c(color_new,color_old)
                print(paste(paste0(ident, "_colors") , "is not matched with cluster (number). we'll use new color"))
            }
    } else {
        color_old <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold", "DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4",  "#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B",  "#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080",  "#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49",  "#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD",  "#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B",  "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9",  "#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000",  "#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
        color_new <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
        clustcol <- c(color_new,color_old)
        print(paste(paste0(ident, "_colors") , "is not in metadata. we'll use new color"))
    }
  } else {
        if (paste0(ident, "_colors") %in% colnames(PRO@meta.data)) {
            color1 <- getColors(PRO, ident)
            color1[[1]] <- color1[[1]][!is.na(color1[[1]])]
            observe <- setdiff(cluster,reference)
            obs_color <- color1[[1]][observe]
            color2 <-  hue_pal()(15)
            ref_color <- color2[1:length(reference)]
            names(ref_color) <- reference
            obs_color <- color1[[1]][observe]
            clustcol <- c(ref_color,obs_color)
            num <- length(unique(cnv$ident))
            color_num <- length(clustcol)
              if (color_num != num){
                color_old <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold", "DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4",  "#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B",  "#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080",  "#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49",  "#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD",  "#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B",  "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9",  "#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000",  "#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
                color_new <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
                clustcol <- c(color_new,color_old)
                print(paste(paste0(ident, "_colors") , "is not matched with cluster (number). we'll use new color"))
            }
    } else {
        color_old <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold", "DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4",  "#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B",  "#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080",  "#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49",  "#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD",  "#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B",  "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9",  "#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000",  "#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
        color_new <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
        clustcol <- c(color_new,color_old)
        print(paste(paste0(ident, "_colors") , "is not in metadata. we'll use new color"))
        }
    }
} else{
  clustcol <- unlist(strsplit(argv$color, split = ","))
}

if (splitby != "F"){
    if (paste0(splitby, "_colors") %in% colnames(PRO@meta.data)) {
        color1 <- getColors(PRO, splitby)
        color1[[1]] <- color1[[1]][!is.na(color1[[1]])]
        clustcol_splitby <- color1[[1]]
        num <- length(unique(cnv[,splitby]))
        print(num)
        color_num <- length(clustcol_splitby)
        print(color_num)
            if (color_num != num){
                color_old <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold", "DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4",  "#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B",  "#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080",  "#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49",  "#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD",  "#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B",  "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9",  "#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000",  "#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
                color_new <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
                clustcol_splitby <- c(color_new,color_old)
                print(paste(paste0(splitby, "_colors") , "is not matched with splitby (number). we'll use new color"))
            }
    } else {
        color_old <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold", "DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4",  "#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B",  "#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080",  "#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49",  "#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD",  "#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B",  "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9",  "#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000",  "#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
        color_new <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
        clustcol_splitby <- c(color_new,color_old)
        print(paste(paste0(splitby, "_colors") , "is not in metadata. we'll use new color"))
    }
}


#caculate the outliners
#------------------------------------------------------------------------------------------------------------
#y <- c()
# for (i in unique(cnv$ident)){
#     data <-subset(cnv,ident=i)
#     y_tmp <- quantile(data$CNVscore)
#     iqr <- IQR(data$CNVscore)
#     ymax <- y_tmp[4]+1.5*iqr
#     y <- c(ymax,y)
# }
# ylim <- max(y)
# ymax <- ylim*1.2

#plot
#--------------------------------------------------------------------------------------------------------------
if (boxplot == "T"){
        p <-boxplot_function(data=cnv,x="ident",y="CNVscore",group="ident",filltype="fill",usecolor=clustcol,
                title_x =ident,title_y ="CNV score",rmoutlier =TRUE,savefig =FALSE)
        if (comparison != "F"){
                p <- boxplot_function(data=cnv,x="ident",y="CNVscore",group="ident",filltype="fill",usecolor=clustcol,
                title_x =ident,title_y ="CNV score",savefig =FALSE) +
                stat_compare_means(comparisons =   comparegroup , label = "p.signif") #add p-vaule
        }
        file <- paste(outdir2,"/",prefix,"_cnv_score_boxplot.pdf",sep="")
        params <- get_params(data=cnv,splitby="F",ident=ident)
        width <- params[[1]];height <- params[[2]]
        plotFunc(p,file,width=width,height=height,png=TRUE)
        if (splitby != "F" & ident =="cluster"){
                p <- boxplot_function(data=cnv,x="ident",y="CNVscore",group=splitby,filltype="fill",#usecolor=clustcol_splitby,
                    title_x ="cluster",title_y ="CNV score",rmoutlier =TRUE,savefig =FALSE,title_legend =splitby)
                p <- p + scale_fill_manual(values= clustcol_splitby)
                file <- paste(outdir2,"/",prefix,".","splitby",splitby,"_cnv_score_boxplot.pdf",sep="")
                params <- get_params(data=cnv,splitby=splitby,ident=ident)
                width <- params[[1]];height <- params[[2]]
                plotFunc(p,file,width=width,height=height,png=TRUE)
        }
}
if (vlnplot == "T"){
        p <- vlnplot_function(data=cnv,x="ident",y="CNVscore",group="ident",filltype="fill",usecolor=clustcol,
                                title_x =ident,title_y ="CNV score")
        if (comparison != "F"){
            p <- p + stat_compare_means(comparisons = comparegroup, label = "p.signif") #add p-vaule
        }  
        file <-paste(outdir2,"/",prefix,"_cnv_score_vlnplot.pdf",sep="")
        params <- get_params(data=cnv,splitby="F",ident=ident)
        width <- params[[1]];height <- params[[2]]
        plotFunc(p,file,width=width,height=height,png=TRUE)
        if (splitby != "F" & ident =="cluster"){
            cnv <- filter(cnv,splitby=splitby,ident="ident")
            p <- vlnplot_function(data=cnv,x="ident",y="CNVscore",group=splitby,filltype="fill",usecolor=clustcol_splitby,
                    title_x ="cluster",title_y ="CNV score",title_legend =splitby)
            file <- paste(outdir2,"/",prefix,".","splitby",splitby,"_cnv_score_vlnplot.pdf",sep="")
            params <- get_params(data=cnv,splitby=splitby,ident=ident)
            width <- params[[1]];height <- params[[2]]
            plotFunc(p,file,width=width,height=height,png=TRUE)
        }
}

# file.copy(from=paste0(scriptdir,"/inferCNV_CNVscore_READEME.pdf"), to = outdir2)


