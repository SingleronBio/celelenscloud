library(gtools)
library(ggplot2)
library(argparser)
library(infercnv)

argv <- arg_parser('')
argv <- add_argument(argv, "--type", help="inferCNV result or inferCNVpy result", default="inferCNV")
argv <- add_argument(argv, "--cnvdir", help="infercnv result dir")
argv <- add_argument(argv, "--ref", help="the ref celltype")
argv <- add_argument(argv, "--cluster", help="the cluster need to analysis, split by ','", default='all')   # 如果是HCL的，用HCL_1的形式
argv <- add_argument(argv, "--outdir", help="out dir", default='./')
argv <- add_argument(argv, "--genesfile", help="gene positation information", default='/Public/Script/shouhou/CNV_gene_Enrich/my_gene_no_sex_chr.txt')
argv <- add_argument(argv, "--cutoff", help="the cutoff line show on the plot", default=0.02, type='float')
argv <- add_argument(argv, "--kgroups", help="k-means number", default= "FALSE")
argv <- add_argument(argv, "--splitby", help = "use sample/group split ident to plot", default = "FALSE")
argv <- add_argument(argv, "--groupfile", help = "column1: group, column2: samples -- split by ',', table split, ", default = "FALSE")
argv <- add_argument(argv, "--objlist", help = "column1: sample, column2: infercnv result obj, table split", default="FALSE")  # 目前只能用于inferCNV结果
argv <- add_argument(argv, "--metalist", help = "column1: sample, column2: meta.xls, table split", default="FALSE")

argv <- parse_args(argv)

#测试默认参数
######################################################
# argv$type <- "inferCNV"
# argv$cnvdir <- "/SGRNJ03/pipline_test/chenkai/CNV_region_py/inferCNV_R_1/02.test.infercnv/"
# argv$outdir <- "/SGRNJ03/pipline_test/chenkai/union_result/r_groupfile"
# argv$ref <- "KCs,TCells"
# argv$cluster <- "Hepatocytes1"
# argv$genesfile <- "/Public/Script/shouhou/CNV_gene_Enrich/my_gene_no_sex_chr.txt"
# argv$cutoff <- 0.05
# argv$groupfile <- "/SGRNJ03/pipline_test/chenkai/group"
#######################################################

type <- argv$type
cnvdir <- argv$cnvdir
ref <- unlist(strsplit(argv$ref, split = ','))
cluster <- unlist(strsplit(argv$cluster,split = ','))
outdir <- argv$outdir
genesfile <- argv$genesfile
cnv_cutoff <- as.numeric(argv$cutoff)
k_nums <- argv$kgroups
splitby <- argv$splitby
objlist <- argv$objlist
metalist <- argv$metalist

dir.create(outdir)

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

scriptDirname = dirnameScript()
readmePath <- file.path(scriptDirname, "readme", "CNV_region-README.pdf")


# function 1: cnv_region function
cnv_region <- function(expr, cnv_cutoff=0.02, genepos, outdir, type) {
    cnvs <- colMeans(expr)
    cnv_sd <- sd(cnvs)
    cnv_prop <- apply(expr, 2, function(x){sum(abs(x)>cnv_cutoff)/length(x)})
  
	
    if (type == "inferCNV") {
        cnv_mean = data.frame(index=c(1:ncol(expr)), cnv=cnvs, prop=cnv_prop, check.names=F)
        cnv_mean$chr = genepos[rownames(cnv_mean),1]
        write.table(cnv_mean,paste0(outdir,'/CNV_average_', i, '.xls'), row.names = T, sep='\t', col.names = NA)
    } else if (type == "inferCNVpy") {
        cnv_mean <- data.frame(index=c(0:(ncol(expr)-1)), cnv=cnvs, prop=cnv_prop, check.names = F)
        cnv_mean <- merge(cnv_mean, genepos, by.x = 'index', by.y = 'pos',)
        cnv_mean$chr <- as.numeric(cnv_mean$chr)			# merge之后chr列貌似变成了字符型,强制转换一下,避免在头尾增加0值时出现NA的情况
        write.table(cnv_mean,paste0(outdir,'/CNV_average_', i, '.xls'), row.names = F, sep='\t', col.names = NA)
	}
    
    
    cnv_mean$pos <- c(0, cnv_mean$cnv[1:(nrow(cnv_mean)-1)] * cnv_mean$cnv[2:nrow(cnv_mean)])      # 相邻基因一正一负会影响画图，通过相邻基因cnv相乘判断
    # 向cnv为一正一负的相邻两行基因中插入一行0值，-1.5是因为scanpy的index从0开始
    gap <- switch(type, inferCNV=0.5, inferCNVpy=1.5)
    if (any(cnv_mean$pos<0)) {
        tmp = data.frame(index=which(cnv_mean$pos<0)-gap,cnv=0,prop=NA,chr=NA, pos=0, check.names=F)   
        cnv_mean = rbind(cnv_mean,tmp); cnv_mean = cnv_mean[order(cnv_mean$index), ]
    }  
    cnv_mean = rbind( c(0,0,0,0,0), cnv_mean, c(max(cnv_mean$index)+1,0,0,0,0) )                   # 一头一尾添加0值,确保geom_polygon()在fill时会填充x轴以上的部分
    xind = cnv_mean$index[sapply(1:22, function(x){min(which(cnv_mean$chr==x))} )]
    cnv_max=max(abs(cnv_mean$cnv))
  
  
    pic_mean = ggplot(cnv_mean, aes(x=index,y=cnv)) + 
        geom_polygon(dat=cnv_mean[cnv_mean$cnv>=0,],aes(x=index,y=cnv),fill='#CC0033',alpha=0.9) +
        geom_polygon(dat=cnv_mean[cnv_mean$cnv<=0,],aes(x=index,y=cnv),fill='#0099CC',alpha=0.9)  +
        geom_hline(yintercept=0,linetype='dashed') + 
        geom_hline(yintercept=c(-cnv_cutoff,cnv_cutoff),linetype='dashed') +
        geom_vline(xintercept= xind,col='grey75') + 
        scale_x_continuous(expand=c(0,0),breaks=xind,labels=as.character(c(1:22))) + 
        theme_classic()  + labs(x='CNV',y='') +
        annotate('text',x=c(-120,-120),y=c(-cnv_max,cnv_max),label=c('Del','Amp') ,angle=90) +
        coord_cartesian(clip = "off",xlim=c(0,max(cnv_mean$index)), ylim=c(-cnv_max,cnv_max))
    
    pic_prop = ggplot(na.omit(cnv_mean), aes(x=index,y=prop)) + geom_line() + 
        geom_polygon(fill='#bfcde0',alpha=0.8) +
        scale_x_continuous(expand=c(0,0),breaks=xind,labels=as.character(c(1:22))) + scale_y_continuous(expand=c(0,0)) +
        geom_vline(xintercept= xind,col='grey75') + 
        theme_classic()  + labs(x='',y='Proportion')
    
    p <- cowplot::plot_grid(pic_mean, pic_prop, nrow=2)
    png(paste0(outdir,"/CNV_infor_", i, ".png"),width=16,height=8,units='in',res=400)
    print(p)
    dev.off()
    pdf(paste0(outdir,"/CNV_infor_", i, ".pdf"),width=16,height=8)
    print(p)
    dev.off()
    ggsave(file=paste0(outdir,"/CNV_infor_", i, ".svg"), plot=p, width=16, height=8)
}


# function 2: chr_pos function
chr_pos <- function(chr) {
    dat <- data.frame(chr=NA,pos=NA)
    dat <- dat[-1,]
    for (i in 1:(nrow(chr)-1)){
        chr_i <- chr[i,1]
        chr_i <- gsub("chr","",chr_i)
        b <- chr[i+1,2]- chr[i,2]
        dat_i <- data.frame(chr=rep(chr_i,b),pos=seq(chr[i,2],chr[i+1,2]-1))
        dat_i$pos <- dat_i$pos
        dat <- rbind(dat,dat_i)
    }
    return(dat)
}


# function 3: 合并expr
merge_expr <- function(objlist) {
    colnames(objlist) <- c("sample", "obj")
	sample <- objlist$sample
    expr_list <- list()
    for (i in 1:nrow(objlist)) {
        obj <- readRDS(objlist$obj[i])
        expr_list[[sample[i]]] <- t(obj@expr.data) - 1 # expr_list 行为细胞，列为基因
    }
    # 对基因取交集
    gene <- colnames(expr_list[[1]])
    for (i in expr_list) {
        gene <- intersect(gene, colnames(i))
    }
    # 筛选每个expr
    for (i in 1:length(expr_list)) {
        expr_list[[i]] <- expr_list[[i]][, gene]
    }
    # 合并expr，保存至expr_list中
    expr_list[["all"]] <- Reduce(rbind, expr_list)

	return(expr_list[["all"]])
}


# function 4: 合并meta表
merge_meta <- function(metalist) {
    colnames(metalist) <- c("sample", "meta")
    meta_list <- list()
    for (i in 1:nrow(metalist)) {
        meta_list[[i]] <- read.table(metalist$meta[i], sep = "\t", comment.char = "", header = T)
    }
    meta <- Reduce(rbind, meta_list)
    return(meta)
}



print("---------------------loading CNV matrix---------------------")
expr_list <- list()

if (objlist == "FALSE") {
	if (type == "inferCNV") {
        #expr_list[["all"]] <- t(read.table(paste0(cnvdir,'/infercnv.observations.txt'),header=T,row.names=1,sep=" ",check.names=F)) - 1
        obj <- readRDS(paste0(cnvdir, "/run.final.infercnv_obj"))
        exp <- t(obj@expr.data[,unlist(obj@observation_grouped_cell_indices)]) - 1   # 只保留anno/observation细胞
        meta <- read.table( Sys.glob( paste0(cnvdir, "/../01.predata/*anno.metadata.xls") ), sep = "\t", comment.char="", header = T) # , row.names = 1 
        rownames(meta) <- meta$barcode		
        expr_list[["all"]] <- exp
        # ------genepos
        genepos <- read.table(genesfile, header = F, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)
    } else {
        exp <- read.table(paste0(cnvdir,"/CNV_matrix.txt"), header=TRUE, row.names=1, check.names = F)
        meta <- read.table(paste0(cnvdir,"/metadata.txt"),header=TRUE, row.names=1, comment.char="",sep="\t")
        exp <- exp[rownames(meta)[!(meta$cluster %in% ref) ],]  
        print(meta[1:3,1:3])
        expr_list[["all"]] <- exp
        # ------chrfile
        chrFile <- read.table(paste0(cnvdir, "/chr_pos.txt"), header = T, col.names = c("chr", "pos"), stringsAsFactors = F)
        chrFile[nrow(chrFile)+1,] <- c("end", ncol(expr_list[['all']]))
        chrFile$pos <- as.numeric(chrFile$pos)
        genepos <- chr_pos(chr=chrFile)
        chrorder <- mixedsort(unique(genepos$chr))
        genepos$chr <- factor(genepos$chr,levels=chrorder)
    }
} else if (objlist != "FALSE") {
	objlist <- read.table(objlist, header = F, sep = "\t")
	metalist <- read.table(metalist, header = F, sep = "\t")
	expr_list[["all"]] <- merge_expr(objlist = objlist)
	meta <- merge_meta(metalist = metalist)
	meta <- meta[!meta$cluster %in% ref,]				# 删除合并meta表中的ref细胞，避免分样本且指定ref时，合并的meta表中有重复的ref细胞
    if("barcode" %in% colnames(meta)) {rownames(meta) <- meta$barcode}
    genepos <- read.table(genesfile, header = F, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)
}

cat("Cells number of anno: ", nrow(expr_list[["all"]]), "\n")


# 是否要kmeans聚类
if (k_nums != "FALSE") {
    set.seed(1)
    kmeans_result <- kmeans(expr_list[["all"]], as.numeric(k_nums))
    kmeans_df <- data.frame(kmeans_class = kmeans_result$cluster)
    kmeans_df$group <- meta[rownames(kmeans_df), gname]
    print(head(kmeans_df))
    table(kmeans_df$kmeans_class, kmeans_df$group)
}

# 是否按样本或组拆分cnv矩阵，分样本存放至expr_list中
if (argv$splitby != "FALSE") {
    splitby <- argv$splitby
    split <- unique(meta[, splitby])
    print(split)
    for (i in split) {
        cell <- rownames(meta)[meta[, splitby] == i]
        print(paste0("Cells number of ", i, ":", length(cell)))
        expr_list[[i]] <- expr_list$all[rownames(expr_list$all) %in% cell, ]
    }
}

# 是否重新组合样本
if (argv$groupfile != "FALSE") {
    groupfile <- read.table(argv$groupfile, header = F, sep = "\t")
	colnames(groupfile) <- c("group", "sample")
    split <- groupfile$group
	print(head(meta))
	print(groupfile)
    for(i in split) {
        print(i)
        cell <- rownames(meta)[meta$sample %in% unlist(strsplit(as.character(groupfile$sample[which(groupfile == i)]), split = ","))]
        print(head(cell))
        expr_list[[i]] <- expr_list$all[rownames(expr_list$all) %in% cell,]
    }
}


# 用于循环
ifelse (argv$splitby != "FALSE" | argv$groupfile != "FALSE", loop <- split, loop <- "all")
cat("loop: ", loop, "\n")

for (i in loop) {
    print(i)
    expr_list[[i]][1:3,1:3]
    expr <- expr_list[[i]]
    print(expr[1:3,1:3])

  # 筛选细胞
    if (all(cluster != "all")) {
    if (k_nums != "FALSE") {													                      # 按kmeans筛选
        cell_id <- rownames(kmeans_df)[kmeans_df$kmeans_class %in% cluster]
        expr <- expr[rownames(expr) %in% cell_id, ] 
    } else if (all(cluster %in% unique(meta$cluster))) {						        # 按celltype筛选
        cell_id <- rownames(meta)[meta$cluster %in% cluster]
        expr <- expr[rownames(expr) %in% cell_id, ]
    } else if (argv$type == "inferCNV" & all(grepl("^HCL_", cluster))) {		# 按HCL筛选
        HCL <- read.table(paste0(cnvdir, "/infercnv.HCL.barcode.txt"), header = F, row.names = 1, sep = "\t", check.names = F, col.names = c("barcode","HCL"))
        cell_id <- rownames(HCL)[HCL$HCL %in% cluster]
        expr <- expr[rownames(expr) %in% cell_id, ]
    }
    }

    print(paste0("Used cells number: ", nrow(expr)))
    if (nrow(expr) != 0) {
        cnv_region(expr=expr, cnv_cutoff=cnv_cutoff, genepos=genepos, outdir=outdir, type=type)
    } else {
        cat("***** Used cells number = 0, skip", i, '*****\n')
    } 
}

#file.copy("/SGRNJ03/pipline_test/chenkai/README/CNV_region-README.pdf",outdir)
# file.copy (from=readmePath, to=outdir)
print("------------------Done--------------------")