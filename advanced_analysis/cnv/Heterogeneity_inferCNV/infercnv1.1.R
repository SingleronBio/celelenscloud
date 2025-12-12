library(Seurat)
library(Matrix)
library(ggplot2)
library(infercnv)
library(argparser)
library (dplyr)
library (stringr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)
library(parallelDist)
set.seed(123)


argv <- arg_parser('')
argv <- add_argument(argv,"--including_clust", help="cluster you want to included in your cnv. like Tcell, Bcell, seperate by ','",default = 'all')
argv <- add_argument(argv,"--rds", help="the absulote path of seurat rds file")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- add_argument(argv,"--version", help="seurat veresion",default = "V3")
argv <- add_argument(argv,"--subset",help="subset data, eg: extract 20% cell from total 10000 cell, print 0.2",default=1)
argv <- add_argument(argv,"--prefix",help="output prefix")
argv <- add_argument (argv,"--spname", help = "sample used to analysis,split by ,",default = "all")
argv <- add_argument (argv,"--split_ref", help = "whether to split sample for ref celltype ,",default = "T")
#step1:creat infercnv object
#argv <- add_argument(argv,"--mat",help = "the name of expression matrix")
#argv <- add_argument(argv,"--anno", help= "the annotation file")
argv <- add_argument(argv,"--ref_ct",help="the cell type as a reference")
argv <- add_argument(argv,"--genefile",help="the gene order file")
#step2: infercnv::run (perform infercnv operations to reveal cnv signal)
##cutoff, windows, min_cells_per_gene, outdir, cluster_by_groups
argv <- add_argument(argv,"--cutoff",help="1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics",default = 0.1)
argv <- add_argument(argv,"--window",help = "Length of the window for the moving average (smoothing)",default=101)
argv <- add_argument(argv,"--mingene",help="minimum number of reference cells requiring expression measurements to include the corresponding gene",default=3)
#argv <- add_argument(argv,"--outdir",help="the output dir")
argv <- add_argument(argv,"--cluster",help ="whether cluster each group, default=FALSE, instead will use k setting cluster number",default = FALSE)
argv <- add_argument(argv,"--kgroups",help="Number of groups in which to break the observations",default = 8)
##HMM
argv <- add_argument(argv,"--hmm",help="when set to True, runs HMM to predict CNV level",default = FALSE)
argv <- add_argument(argv,"--hmmtype",help="HMM model type.Options: (i6 or i3)",default = 'i6')
argv <- add_argument(argv,"--BayesMaxP",help="maximum P(Normal) allowed for a CNV prediction according to BayesNet. note zero turns it off",default = 0.5)
argv <- add_argument(argv,"--analysis_mode",help="Grouping level for image filtering or HMM predictions,options(samples|subclusters|cells)",default = "samples")
argv <- add_argument(argv,"--tumor_subcluster_partition_method",help="method for defining tumor subclusters. Options('random_trees','qnorm')",default = "qnorm")
##denoise
argv <- add_argument(argv,"--denoise",help = "de-noising filters",default = TRUE)
argv <- add_argument(argv,"--nfilter",help="de-nosing method",default = NA)
argv <- add_argument(argv, "--sd_amplifier",help = "Noise is defined as mean(reference_cells) +- sdev(reference_cells) * sd_amplifier",default=1.5)
argv <- add_argument(argv,"--logistic",help="de-noise method",default = FALSE)
argv <- parse_args(argv)

prefix <- argv$prefix
rdsfile <- argv$rds
#clusterin <- unlist(strsplit(argv$including_clust,split=","))
version <- argv$version
outdir <- argv$outdir;dir.create(outdir)
subcell <- as.numeric(argv$subset)
spname <- unlist (strsplit (argv$spname,split = ","))
ref_ct <- argv$ref_ct
split_ref <- unlist (strsplit(argv$split_ref,","))

##cutoff, windows, min_cells_per_gene, outdir, cluster_by_groups
cutoff <- as.numeric(argv$cutoff)
window <- as.numeric(argv$window)
mingene <- as.numeric(argv$mingene)
# outdir <- argv$outdir
cluster <- argv$cluster
k <- as.numeric(argv$kgroups)
##HMM
hmm <- argv$hmm
hmmtype <- argv$hmmtype
BayesMaxP <- as.numeric(argv$BayesMaxP)
analysis_mode <- argv$analysis_mode
tumor_subcluster_partition_method <- argv$tumor_subcluster_partition_method
##denoise
denoise <- argv$denoise
filter <- argv$nfilter
sd_amplifier <- as.numeric(argv$sd_amplifier)
logistic <- argv$logistic

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

######prepare the matrix and annotation file##########
if (version == 'V3') {
    rds <- readRDS (rdsfile)
    expm <- rds@assays$RNA@data
}

if(argv$including_clust == "all" ) {
    clusterin <- levels (rds)
}else{
    clusterin <- unlist(strsplit(argv$including_clust,split=","))
}
# if (spname == "all") {
#     spname <- unique (rds@meta.data$sample)
# }
if ("all" %in% spname) {
    spname <- unique (rds@meta.data$sample)
}

all <- unlist (strsplit (clusterin,","))
ref_ct <- unlist (strsplit (ref_ct,","))
obs <- setdiff (all,ref_ct)
sub <- subset (rds,idents = ref_ct)
## get all sample ref  #### 
if (split_ref == "F") {                                  ##不筛选样本，选取所有样本指定ref 细胞类型
    if (nrow (sub@meta.data) <5000) {
        ref.meta <- sub@meta.data
        ref.meta$celltype <- as.character (sub@active.ident)
        ref.anno <- data.frame (barcode = rownames (ref.meta),celltype = ref.meta$celltype)
        ref.anno <- ref.anno[order(ref.anno$celltype),]
        ref.barcode <- ref.anno$barcode
    }else {
        sub <- subset (sub,downsample = 5000/(length(levels (sub))))
        ref.meta <- sub@meta.data
        ref.meta$celltype <- as.character (sub@active.ident)
        ref.anno <- data.frame (barcode = rownames (ref.meta),celltype = ref.meta$celltype)
        ref.anno <- ref.anno[order(ref.anno$celltype),]
        ref.barcode <- ref.anno$barcode
    }
    print (unique (sub$sample))
}else if (split_ref == "T") {
    # if (spname != "all") {                               ## ref 选取与obs 细胞样本来源保持一致
        sub <- subset (sub,subset = sample %in% spname)
        if (nrow (sub@meta.data) <5000) {
            ref.meta <- sub@meta.data
            ref.meta$celltype <- as.character (sub@active.ident)
            ref.anno <- data.frame (barcode = rownames (ref.meta),celltype = ref.meta$celltype)
            ref.anno <- ref.anno[order(ref.anno$celltype),]
            ref.barcode <- ref.anno$barcode
        }else {
            sub <- subset (sub,downsample = 5000/(length(levels (sub))))
            ref.meta <- sub@meta.data
            ref.meta$celltype <- as.character (sub@active.ident)
            ref.anno <- data.frame (barcode = rownames (ref.meta),celltype = ref.meta$celltype)
            ref.anno <- ref.anno[order(ref.anno$celltype),]
            ref.barcode <- ref.anno$barcode
        }
        print (unique (sub$sample))        
    # }else {
    #     print ("--spname need to select sample")    
    # }
}else {
    sub <- subset (sub,subset = sample %in% split_ref)      #### 选取特定样本ref 细胞类型为ref
    if (nrow (sub@meta.data) <5000) {
        ref.meta <- sub@meta.data
        ref.meta$celltype <- as.character (sub@active.ident)
        ref.anno <- data.frame (barcode = rownames (ref.meta),celltype = ref.meta$celltype)
        ref.anno <- ref.anno[order(ref.anno$celltype),]
        ref.barcode <- ref.anno$barcode
    }else {
        sub <- subset (sub,downsample = 5000/(length(levels (sub))))
        ref.meta <- sub@meta.data
        ref.meta$celltype <- as.character (sub@active.ident)
        ref.anno <- data.frame (barcode = rownames (ref.meta),celltype = ref.meta$celltype)
        ref.anno <- ref.anno[order(ref.anno$celltype),]
        ref.barcode <- ref.anno$barcode
    }
    print (unique (sub$sample))       
}

#########################################
print ("############get obs information###############")
pro <- rds
## get obs information
# if (spname != "all") {
if (all(spname != "all")) {
    rds <- subset (rds,subset=sample %in% spname)
    print (spname)
    print (unique(as.character(rds@meta.data$sample)))
    celltypes <- levels (rds)

}else {
    celltypes <- levels (rds)
}
## 判断res  细胞类型是否存在
ref.use <- intersect (ref_ct,celltypes)
if ( ! any(ref.use  %in% levels (rds))) {
    print ("Error:No ref celltype for infercnv")
    quit()
}
## 判断obs 细胞类型是否存在
obs.use <- intersect (obs,celltypes)
if ( ! any(obs.use %in% levels (rds))) {
    print ("Error:No obs celltype for infercnv")
    quit()
}
sub.obs <- subset (rds,idents = obs.use)
if (subcell !=  1 ){
    sub.obs <- subset (sub.obs, downsample = max (table (sub.obs@active.ident)) * subcell)
}

##多个 obs  细胞类型重命名为一个
if(cluster == FALSE) {
    if(length(obs.use) > 1) {
        sub.obs <- SetIdent (sub.obs,cells = colnames (sub.obs),value = "observations")
    }
}


obs.meta <- sub.obs@meta.data
obs.meta$celltype <- as.character (sub.obs@active.ident)

obs.anno <- data.frame (barcode = rownames (obs.meta),celltype = obs.meta$celltype)
obs.anno <- obs.anno[order(obs.anno$celltype),]
obs.barcode <- obs.anno$barcode


df <- rbind (ref.anno,obs.anno)
barcode  <- c(ref.barcode,obs.barcode)
mat <- as.matrix(expm[,barcode])

### save predate
predio <- file.path(outdir,"01.predata")
dir.create (predio)

write.table(mat,file.path(predio,paste0(prefix,"_expm.txt")),quote=F, sep="\t")
colnames(df) <- NULL
write.table(df, file.path(predio,paste0(prefix,"_anno.txt")),quote=F, sep="\t",row.names=F)
## anno meta
anno.meta <- pro@meta.data[df[,1],] %>% rownames_to_column (var = "barcode")
write.table(anno.meta, file.path(predio,paste0(prefix,"_anno.metadata.xls")),quote=F, sep="\t",row.names=F,col.names = T)
##creat infercnv object
mat <- file.path(predio,paste0(prefix,"_expm.txt"))
anno <- file.path(predio,paste0(prefix,"_anno.txt"))
genefile <- argv$genefile
## 判断输入ref obs 细胞是否存在
# all <- unlist (strsplit (clusterin,","))
# ref <- unlist(strsplit(argv$ref_ct,split=","))
# obs <- setdiff(all,ref)
##判断写在了上面
use.ref <- intersect (ref_ct,unique (df[,2]))
# if (length(use.ref) < 1) {
#     print ("Error:No Ref celltype for infercnv")
#     quit()
# }
obs_ct <- intersect (obs,unique (df[,2]))
# if (length(obs_ct) < 1) {
#     print ("Error:No obs celltype for infercnv")
#     quit()
# }
## 判断obs 细胞数目是否>50,低于50 不进行infercnv 分析
# data.obs <- df[which(df[,2]%in% obs_ct),]
if (nrow(obs.anno) < 50) {
    print (obs.anno)
    print ("Error:obs cell number <50")
    quit()
}

out.infercnv <- file.path (outdir,paste0("02.",prefix,".infercnv"))
dir.create(out.infercnv,showWarnings = T)
#setp1:create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=mat,
                                    annotations_file=anno,
                                    gene_order_file=genefile,
                                    ref_group_names=c(use.ref),
                                    delim = "\t")
                                                                
# step2:perform infercnv operations to reveal cnv signal (首先选择聚类的方式，然后选择一种降噪方式进行降噪)
if (cluster != FALSE) {
    if (is.na(filter) != TRUE) {    
        infercnv_obj = infercnv::run(infercnv_obj,
                                    cutoff=cutoff,
                                    min_cells_per_gene = mingene,
                                    window_length = window,  
                                    out_dir=out.infercnv,
                                    cluster_by_groups = TRUE,
                                    HMM = hmm,
                                    HMM_type = hmmtype,
                                    BayesMaxPNormal = BayesMaxP,
                                    analysis_mode = analysis_mode,
                                    tumor_subcluster_partition_method = tumor_subcluster_partition_method,
                                    denoise=denoise,
                                    png_res = 500,
                                    noise_filter = filter)
    }else if(logistic != FALSE){
        infercnv_obj = infercnv::run(infercnv_obj,
                                    cutoff=cutoff,
                                    min_cells_per_gene = mingene,
                                    window_length = window,  
                                    out_dir=out.infercnv,
                                    cluster_by_groups = TRUE,
                                    HMM = hmm,
                                    HMM_type = hmmtype,
                                    BayesMaxPNormal = BayesMaxP,
                                    analysis_mode = analysis_mode,
                                    tumor_subcluster_partition_method = tumor_subcluster_partition_method,
                                    denoise=denoise,
                                    sd_amplifier = sd_amplifier,
                                    png_res = 500,
                                    noise_logistic = TRUE)
    }else {
        infercnv_obj = infercnv::run(infercnv_obj,
                                    cutoff=cutoff,
                                    min_cells_per_gene = mingene,
                                    window_length = window,  
                                    out_dir=out.infercnv,
                                    cluster_by_groups = TRUE,
                                    HMM = hmm,
                                    HMM_type = hmmtype,
                                    BayesMaxPNormal = BayesMaxP,
                                    analysis_mode = analysis_mode,
                                    tumor_subcluster_partition_method = tumor_subcluster_partition_method,
                                    denoise=denoise,
                                    png_res = 500,
                                    sd_amplifier = sd_amplifier)
    }                         
}
if (cluster == FALSE) {
    if (is.na(filter) != TRUE) {    
        infercnv_obj = infercnv::run(infercnv_obj,
                                    cutoff=cutoff,
                                    min_cells_per_gene = mingene,
                                    window_length = window,  
                                    out_dir=out.infercnv,
                                    k_obs_groups = k,
                                    HMM = hmm,
                                    HMM_type = hmmtype,
                                    BayesMaxPNormal = BayesMaxP,
                                    analysis_mode = analysis_mode,
                                    tumor_subcluster_partition_method = tumor_subcluster_partition_method,
                                    denoise=denoise,
                                    png_res = 500,
                                    noise_filter = filter)
    }else if(logistic != FALSE){
        infercnv_obj = infercnv::run(infercnv_obj,
                                    cutoff=cutoff,
                                    min_cells_per_gene = mingene,
                                    window_length = window,  
                                    out_dir=out.infercnv,
                                    k_obs_groups = k,
                                    HMM = hmm,
                                    HMM_type = hmmtype,
                                    BayesMaxPNormal = BayesMaxP,
                                    analysis_mode = analysis_mode,
                                    tumor_subcluster_partition_method = tumor_subcluster_partition_method,
                                    denoise=denoise,
                                    sd_amplifier = sd_amplifier,
                                    png_res = 500,
                                    noise_logistic = TRUE)
    }else {
        infercnv_obj = infercnv::run(infercnv_obj,
                                    cutoff=cutoff,
                                    min_cells_per_gene = mingene,
                                    window_length = window,  
                                    out_dir=out.infercnv,
                                    k_obs_groups = k,
                                    HMM = hmm,
                                    HMM_type = hmmtype,
                                    BayesMaxPNormal = BayesMaxP,
                                    analysis_mode = analysis_mode,
                                    tumor_subcluster_partition_method = tumor_subcluster_partition_method,
                                    denoise=denoise,
                                    png_res = 500,
                                    sd_amplifier = sd_amplifier)
    }                     
}
## plot pdf

final.obj <- readRDS (file.path(out.infercnv,"run.final.infercnv_obj"))
if (cluster == FALSE)  {
    # final.obj@tumor_subclusters$hc[["all_observations"]] =  final.obj@tumor_subclusters$hc[["CancerCells_Like"]]
    if (length(obs.use) <= 1) {
    final.obj@tumor_subclusters$hc[["all_observations"]] =  final.obj@tumor_subclusters$hc[[obs.use]]
    }else{
        final.obj@tumor_subclusters$hc[["all_observations"]] =  final.obj@tumor_subclusters$hc$observations
    }
}

if (cluster == FALSE) {
    infercnv::plot_cnv(final.obj,
         k_obs_groups=as.numeric(k),
         cluster_by_groups=FALSE,
         out_dir=out.infercnv,
         x.center=1,
         write_expr_matrix=TRUE,
         write_phylo=TRUE,
         x.range= "auto",
         title="inferCNV",
         output_filename="infercnv",
         output_format="pdf")

    infercnv::plot_cnv(final.obj,
         k_obs_groups=as.numeric(k),
         cluster_by_groups=FALSE,
         out_dir=out.infercnv,
         x.center=1,
         write_expr_matrix=FALSE,
         write_phylo=FALSE,
         x.range= "auto",
         title="inferCNV",
         output_filename="infercnv",
         output_format="png")
}
if (cluster == TRUE) {
    infercnv::plot_cnv(final.obj,
         k_obs_groups=8,
         cluster_by_groups=TRUE,
         out_dir=out.infercnv,
         x.center=1,
         x.range= "auto",
         write_expr_matrix=TRUE,
         write_phylo=TRUE,
         title="inferCNV",
         output_filename="infercnv",
         output_format="pdf")

    infercnv::plot_cnv(final.obj,
         k_obs_groups=8,
         cluster_by_groups=TRUE,
         out_dir=out.infercnv,
         x.center=1,
         x.range= "auto",
         write_expr_matrix=FALSE,
         write_phylo=FALSE,
         title="inferCNV",
         output_filename="infercnv",
         output_format="png")
}

#else {
#     infercnv::plot_cnv(final.obj,
#          k_obs_groups=as.numeric(k),
#          cluster_by_groups=FALSE,
#          out_dir=out.infercnv,
#          x.center=1,
#          write_expr_matrix=TRUE,
#          write_phylo=TRUE,
#          x.range= "auto",
#          title="inferCNV",
#          output_filename="infercnv",
#          output_format="pdf")

#     infercnv::plot_cnv(final.obj,
#          k_obs_groups=as.numeric(k),
#          cluster_by_groups=FALSE,
#          out_dir=out.infercnv,
#          x.center=1,
#          write_expr_matrix=FALSE,
#          write_phylo=FALSE,
#          x.range= "auto",
#          title="inferCNV",
#          output_filename="infercnv",
#          output_format="png")
# }
# #####
if (file.exists (file.path(out.infercnv,"General_HCL_1_members.txt"))) {
    fb <- data.frame ()
    files <- list.files (out.infercnv,"*members.txt")
    for (fi in files) {
        id <- str_replace (fi,"General_","") %>% str_replace ("_members.txt","")
        print (id)
        data <- read.table (file.path(out.infercnv,fi),header = T,sep = ' ')
        cellname <- data.frame (barcode = rownames (data),celltype=id)
        fb <- bind_rows (fb,cellname)
    }    
    write.table (fb,file.path(out.infercnv,paste0("infercnv.HCL.barcode.txt")),row.names = F,col.names = F,quote = F,sep='\t')
    fb.ref <- bind_rows (ref.anno,fb)
    write.table (fb.ref,file.path(out.infercnv,paste0("infercnv.HCL.withanno.barcode.txt")),row.names = F,col.names = F,quote = F,sep='\t')    
}else{
    fb <- read.table (file.path(predio,paste0(prefix,"_anno.txt")),header = F,sep ='\t')
    colnames (fb) <- c("barcode","celltype")
    fb1 <- fb[which(fb$celltype %in% obs_ct),]
    write.table (fb1,file.path(out.infercnv,paste0("infercnv.HCL.barcode.txt")),row.names = F,col.names = F,quote = F,sep='\t')
    file.copy (file.path(predio,paste0(prefix,"_anno.txt")),file.path(out.infercnv,paste0("infercnv.HCL.withanno.barcode.txt")),overwrite = TRUE)
}

out.rs <- file.path(out.infercnv,prefix)
dir.create(out.rs)

scriptDirname <- dirnameScript() # scriptDirname = '/SGRNJ03/PiplineTest01/Pipline_test/qianzhengzong/module/v2.0/Signature_GSVA'
print (scriptDirname)
# readmePath <- file.path(scriptDirname, 'inferCNV_heatmap_READEME.pdf')
# file.copy(readmePath, out.rs)


setwd (out.infercnv)
system("cut -d ' ' -f 2,3 infercnv.observation_groupings.txt | uniq |grep -v 'Group'| sed '1!G;h;$!d'  > infercnv.HCLorder.txt")

#file.copy ("/Public/Script/shouhou/README/tumor-readme/inferCNV_heatmap_READEME.docx",out.rs)
# file.copy ("inferCNV_heatmap_READEME.pdf",out.rs)
file.copy ("infercnv.png",out.rs,overwrite = TRUE)
file.copy ("infercnv.svg",out.rs,overwrite = TRUE)
file.copy ("infercnv.HCLorder.txt",out.rs,overwrite = TRUE)

## replot cnv heatmap
print ("replot CNV heatmap")
out.heatmap <- file.path(out.infercnv,'replot.heatmap')
dir.create (out.heatmap)
final.obj <- file.path(out.infercnv,"run.final.infercnv_obj")
fgroup <- file.path(out.infercnv,"infercnv.observation_groupings.txt")

# system (paste0("Rscript /SGRNJ03/pipline_test/fangsheng/inferCNV/heatmap/inferCNV_heatmap/inferCNV_heatmap.R --input ",out.infercnv,"  --outdir ",out.heatmap," --prefix ",prefix," --chr ",genefile," --rds ",rdsfile,""))
heatmapPath <- file.path(scriptDirname, 'inferCNV_heatmap.R')
system (paste0("Rscript ",heatmapPath, "  --input ",out.infercnv,"  --outdir ",out.heatmap," --prefix ",prefix," --chr ",genefile," --rds ",rdsfile,""))
