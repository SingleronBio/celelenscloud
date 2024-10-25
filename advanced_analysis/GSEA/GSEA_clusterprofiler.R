set.seed(10)
suppressMessages({
library(clusterProfiler)
library(plyr)
library(dplyr)
library(enrichplot)
library(argparser)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(ggupset)
library(ggtree)
library(limma)
library(Seurat)
library(msigdbr)
})

argv <- arg_parser('Pathway enrichment analysis with GSEA by clusterprofiler package')
argv <- add_argument(argv,'--prefix',help='the prefix of outfile',default='Enriched')
argv <- add_argument(argv,'--outdir',help='the output dir',default='.')
argv <- add_argument(argv,'--rds', help="the rds file")
argv <- add_argument(argv,'--diff', help="the diff file",default="F")
argv <- add_argument(argv,'--subcluster', help="subcluster list ,split by ,", default="all")
argv <- add_argument(argv,'--sample', help="sample list ,split by ,", default="all")
argv <- add_argument(argv,'--gname', help="sample list ,split by ,", default="F")
argv <- add_argument(argv,'--DefaultAssay', help="sub cluster seurat[integrated|]",default="RNA")
argv <- add_argument(argv,'--diffcluster', help="compare cluster  name:1:4vs2,3:4vs4,or T:NKvsB,NKvsB; split by ,",default="F")
argv <- add_argument(argv,'--group', help="compare group name:G1vsG2; split by ,",default="F")
argv <- add_argument(argv,'--minpct', help="only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations.",default="0.05")
argv <- add_argument(argv,'--gmt',help='Contains one or more gene sets. For each gene set, gives the gene set name and list of genes')
argv <- add_argument(argv,'--msigdb',help="use gene sets from msigdb database.",default="F")
argv <- add_argument(argv,'--species', help='the species:"Homo sapiens", "Mus musculus" etc',default="NULL")
argv <- add_argument(argv,'--category',help='the library: H,C1,C2,C3,C4,C5,C6,C7,C8',default=NULL)
argv <- add_argument(argv,'--subcategory', help='the sub library:BP,CC,MF. It is the next level of category',default=NULL)
argv <- add_argument(argv,'--subpathway',help='the sub pathway,It is the next level of subcategory',default="NULL")
argv <- add_argument(argv,'--cutoff',help='the cutoff value type to filter pathways, usable: pvalue,p.adjust, qvalues',default='p.adjust')
argv <- add_argument(argv,'--threshold',help='the cutoff values to filter pathways, usable: 0.05, 0.25',default='0.05')
argv <- add_argument(argv,'--pathway',help='compare pathways, name: "Cytokine-cytokine receptor interaction","Huntington disease","Oxidative phosphorylation"',default='F')
argv <- add_argument(argv,'--topshow',help='the number of top pathways to show',default=10,type='numeric')
argv <- add_argument(argv,'--minGSSize',help='the min gene number for each geneSet',default=10,type='numeric')
argv <- add_argument(argv,'--maxGSSize',help='the max gene number for each geneSet',default=500,type='numeric')
argv <- add_argument(argv,'--selectpathway',help='the select pathway file (only use in step2)')
argv <- add_argument(argv,'--saverds',help='save enrich rds or not',default=F, type='logit')
argv <- parse_args(argv)

gmt <- argv$gmt
msigdb <- argv$msigdb
species <- argv$species
category <- argv$category
subcategory <- argv$subcategory
subpathway <- argv$subpathway
prefix <- argv$prefix
rdsfile<-argv$rds
diff <- argv$diff 
do_Assay<-as.character(argv$DefaultAssay)
group <- unlist(strsplit(argv$group,split=","))
subcluster <- argv$subcluster
diffcluster <- argv$diffcluster
samplelist <- argv$sample
gnamelist <- argv$gname
outdir <- argv$outdir
minpct_u <- as.numeric(argv$minpct)
cutoff <- argv$cutoff
threshold <- as.numeric(argv$threshold)
pathway <- unlist(strsplit(argv$pathway,split=","))
topn <- argv$topshow
minGSSize  <- argv$minGSSize
maxGSSize  <- argv$maxGSSize
saverds   <- argv$saverds
selepath  <- argv$selectpathway
if ( !is.na(selepath) ) { selepath <- read.csv(selepath, header=F, sep='\t')$V1 }
if ( !is.na(selepath) ) { path <- paste0(outdir,'/',prefix,'/selected/') } else { path <- paste0(outdir,'/',prefix,'/') }
if(!dir.exists(path)){dir.create(path,recursive=T) ; print(paste0('create dir ',path) ) }else{ print('path existed !')}


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
 
scriptdir = dirnameScript() ## 将脚本所在目录赋予scriptdir变量

#function1
gsInfo <- function(object, geneSetID) {
    geneList <- object@geneList

    if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]

    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList

    df$Description <- object@result[geneSetID, "Description"]
    return(df)
}
gseaScores <- getFromNamespace("gseaScores", "DOSE")
ggtable <- function(d, p = NULL) {
    ggplotify::as.ggplot(tableGrob2(d, p))
}

gseareplot <- gseareplot<-function(x, geneSetID, title = "", color="green", base_size = 11,
                      rel_heights=c(1.5, .5, 1), subplots = 1:3,
                      pvalue_table = FALSE, ES_geom="line") {
    ES_geom <- match.arg(ES_geom, c("line", "dot"))

    geneList <- position <- NULL ## to satisfy codetool

    if (length(geneSetID) == 1) {
        gsdata <- gsInfo(x, geneSetID)
    } else {
        gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
    }
    p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
        theme_classic(base_size) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1)) +          #边框，背景线条
              geom_hline(yintercept = 0,color="grey" ) +                                 #y轴0值线
        scale_x_continuous(expand=c(0,0))

    if (ES_geom == "line") {
        es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description),
                              size=1)
    } else {
        es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description),
                               size=1, data = subset(gsdata, position == 1))
    }

    p.res <- p + es_layer +
        theme(legend.position = c(.8, .8), legend.title = element_blank(),
              legend.background = element_rect(fill = "transparent"))

    p.res <- p.res + ylab("Enrichment Score (ES)") +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.line.x=element_blank(),
              plot.margin=margin(t=.2, r =.2, b=0, l=.2, unit="cm"), #上边距、右边距、下边距、左边距
			  axis.title.y=element_text(size=14), #y轴标题大小
			  axis.text=element_text(size=10),  #坐标轴字体大小
			  plot.title = element_text(hjust = 0.5, size=15,face="bold")) #上图大标题居中 
    i <- 0
    for (term in unique(gsdata$Description)) {
        idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
        gsdata[idx, "ymin"] <- i
        gsdata[idx, "ymax"] <- i + 1
        i <- i + 1
    }
    p2 <- ggplot(gsdata, aes_(x = ~x)) +
        geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
        xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
        theme(legend.position = "none",
              plot.margin = margin(t=-.1, b=0,unit="cm"),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.line.x = element_blank(),
			  panel.border = element_rect(colour = "black", fill=NA, size=1)) +   #中间图边框颜色
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0))

    if (length(geneSetID) == 1) {
        v <- seq(1, sum(gsdata$position), length.out=9)
        inv <- findInterval(rev(cumsum(gsdata$position)), v)
        if (min(inv) == 0) inv <- inv + 1

        col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))

        ymin <- min(p2$data$ymin)
        yy <- max(p2$data$ymax - p2$data$ymin) * .3
        xmin <- which(!duplicated(inv))
        xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
        d <- data.frame(ymin = ymin, ymax = yy,
                        xmin = xmin,
                        xmax = xmax,
                        col = col[unique(inv)])
        p2 <- p2 + geom_rect(
                       aes_(xmin=~xmin,
                            xmax=~xmax,
                            ymin=~ymin,
                            ymax=~ymax,
                            fill=~I(col)),
                       data=d,
                       alpha=.9,
                       inherit.aes=FALSE)
    }

    df2 <- p$data #data.frame(x = which(p$data$position == 1))
    df2$y <- p$data$geneList[df2$x]
    new_df2 <- subset(df2, geneList > 0)
	zero_score_position <- length(rownames(new_df2))-1
    p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                              color="grey")+ geom_vline (xintercept = zero_score_position,color="grey",lty=2)       #下图中间颜色，下图加线条为虚线，
    p.pos <- p.pos + ylab("Ranked List Metric") +
        xlab("Rank in Ordered Dataset") +
        theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"),
		      axis.title.y=element_text(size=14), #y轴标题大小
			  axis.title.x=element_text(size=14),
			  axis.text=element_text(size=10))+  #坐标轴字体大小)
		annotate("text",x=zero_score_position/2,y=1,label="Pos (Positively Correlated)",color="red",size=5)+       #下图中间文字"正调控""负调控"
		annotate("text",x=length(rownames(df2))*0.75,y=-1,label="Neg (Negatively Correlated)",color="blue",size=5)+
		annotate("text",x=zero_score_position,y=0.3,label=paste("Zero score at ",zero_score_position),color="black",size=5)

    if (!is.null(title) && !is.na(title) && title != "")
        p.res <- p.res + ggtitle(title)

    if (length(color) == length(geneSetID)) {
        p.res <- p.res + scale_color_manual(values=color)
        if (length(color) == 1) {
            p.res <- p.res + theme(legend.position = "none")
            p2 <- p2 + scale_color_manual(values = "black")
        } else {
            p2 <- p2 + scale_color_manual(values = color)
        }
    }
#	print(x)
    if (pvalue_table) {
        pd <- x[geneSetID, c("NES", "pvalue", "p.adjust","enrichmentScore")]
		rownames(pd)=NULL
 
		tpd <- t(round(pd,3))#四舍五入到三位数
		NES=sprintf("%0.3f",tpd[1,1])
		Pvalue=sprintf("%0.3f",tpd[2,1])
		P.adjust=sprintf("%0.3f",tpd[3,1])
		enrichmentScore=as.numeric(format(tpd[4,1],digits =3,nsmall = 2))
		if (NES>=0){p.res <- p.res + theme(legend.position = "none") +
			 annotate("text",x=length(rownames(df2))*0.7,y=enrichmentScore*0.9,hjust=0,label=paste("NES: ",NES),color="black",size=6)+       #上图中间文字
			 annotate("text",x=length(rownames(df2))*0.7,y=enrichmentScore*0.8,hjust=0,label=paste("Pvalue: ",Pvalue),color="black",size=6)+
			 annotate("text",x=length(rownames(df2))*0.7,y=enrichmentScore*0.7,hjust=0,label=paste("P.adjust: ",P.adjust),color="black",size=6)
			 }else{ p.res <- p.res + theme(legend.position = "none") +
			 annotate("text",x=100,y=enrichmentScore*0.7,hjust=0,label=paste("NES: ",NES),color="black",size=6)+       #上图中间文字
			 annotate("text",x=100,y=enrichmentScore*0.8,hjust=0,label=paste("Pvalue: ",Pvalue),color="black",size=6)+
			 annotate("text",x=100,y=enrichmentScore*0.9,hjust=0,label=paste("P.adjust: ",P.adjust),color="black",size=6)}
}

    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] +
        theme(axis.line.x = element_line(),
              axis.ticks.x=element_line(),
              axis.text.x = element_text())

    if (length(subplots) == 1)
        return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2,
                                                        l=.2, unit="cm")))


    if (length(rel_heights) > length(subplots))
        rel_heights <- rel_heights[subplots]

    aplot::plot_list(gglist = plotlist, ncol=1, heights=rel_heights)
}
# function 2
# function 3
set_panel_size <- function(p=NULL, g=ggplot2::ggplotGrob(p), width=NULL, height=NULL) {
  panels <- grep("panel", g$layout$name)
  panel_index_w <- unique(g$layout$l[panels])
  panel_index_h <- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  if (!is.null(width))  {g$widths[panel_index_w] <- rep(width, nw) }
  if (!is.null(height)) {g$heights[panel_index_h] <- rep(height, nh)}
  return(g)
}
# function 4
get_term <- function(obj, topn, selterms){
	if (0<nrow(y_sig) & nrow(y_sig)<topn){topn<-nrow(y_sig)}else{topn<-topn}
	obj <- obj %>% arrange( get(cutoff) ) %>% slice(1:topn)
	if ( !is.na(selterms) ){ uterms = selterms}else{uterms = obj@result$Description}
	return(uterms)
}

# function 5
enrich_plots <- function(obj, uterms){
  if (length(uterms) < 2) { return('pathways too less to show , skip !!!') }
  obj = filter(obj, obj@result$Description %in% uterms)
  obj = enrichplot::pairwise_termsim(obj)
  len = ceiling(max(nchar(uterms),na.rm=T)/10)
  #
  netplt <- cnetplot(obj, showCategory= uterms, colorEdge=T, foldChange=geneList) + theme(aspect.ratio=1)
  st_plot(netplt, paste0(path,'/',prefix,'_GSEA_pathway_netplot'), height=12, width=12+len*ceiling(length(uterms)/16), res=150 )
  #
  hetplt <- heatplot(obj, showCategory= uterms, foldChange=geneList)
  wid <- len + length(unique(hetplt$data$Gene))/9
  hei <- length(uterms)*0.5
  hetplt = ggarrange(set_panel_size(hetplt, height=unit(length(uterms)*0.5,'cm') ))
  st_plot(hetplt, paste0(path,'/',prefix,'_GSEA_pathway_heatplot'), height=1.5+hei, width=min(90, wid)  )
  #
  if (length(uterms)>2){ 
    mapplt <- emapplot(obj, showCategory= uterms)
    st_plot(mapplt, paste0(path,'/',prefix,'_GSEA_pathway_mapplot'), height=10, width=12 )
  }else{print("mapplot failed")}
  #
  upsetplot <- upsetplot(obj, n=10)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),aspect.ratio=1,axis.text=element_text(size=9))
  st_plot(upsetplot, paste0(path,'/',prefix,'_GSEA_pathway_upsetplot'), height=6, width=6+len )
  #
  tryCatch({treeplot <- treeplot(obj,nWords=0,showCategory= uterms)}, warning = function(w){
    treeplot <- treeplot(obj,nWords=0,showCategory= uterms, nCluster = 3)
    st_plot(treeplot, paste0(path,'/',prefix,'_GSEA_pathway_treeplot'), height=6, width=6+len )
  }, error = function(e){print("treeplot failed")})
  #
  dotplt <- dotplot(obj,showCategory=uterms,color=cutoff,orderBy=cutoff,decreasing=F,label_format=100) + theme(legend.box='horizontal')
  if (nrow(obj) <= 3) { dotplt <- dotplt + scale_size_continuous(breaks=dotplt$data$Count, range=c(2,5)) + theme(legend.key.height=unit(4,'mm')) }
  dotplt = ggarrange(set_panel_size(dotplt, width=unit(3,'in')))
  st_plot(dotplt, paste0(path,'/',prefix,'_GSEA_pathway_dotplot'), height=1+0.23*length(uterms), width=5.5+len )
  #
  sourceScript("plotGseaTable.R")
  st_plot(plotGseaTable(egmt@geneSets[names(egmt@geneSets) %in% topshow],geneList,egmt@result,gseaParam = 0.5), paste0(path,'/',prefix,'_GSEA_pathway_tableplot'), height=1+0.23*length(uterms), width=5.5+len )
  png(paste0(path,'/',prefix,'_GSEA_pathway_tableplot.png'),height=100+23*length(uterms), width=100+200*len)
  print(plotGseaTable(egmt@geneSets[names(egmt@geneSets) %in% topshow],geneList,egmt@result,gseaParam = 0.5))
  dev.off()
}

#find diff expressed genes,logfc=0
if (diff == "F"){
# system(paste0("source /Public/Software/anaconda3/etc/profile.d/conda.sh; conda activate /SGRNJ/Public/Software/conda_env/Seuratv3; Rscript /SGRNJ03/PiplineTest01/Software_test/changwenjing/gsea/diff_exp_rds.R --prefix ", prefix," --outdir ",outdir," --rds ",rdsfile," --subcluster ",subcluster," --sample ",samplelist," --gname ",gnamelist," --diffcluster ",diffcluster," --group ",group," --minpct ",minpct_u," --logfcthreshold 0"))
system(paste0('Rscript ', scriptdir , '/',"diff_exp_rds.R --prefix ", prefix," --outdir ",outdir," --rds ",rdsfile," --subcluster ",subcluster," --sample ",samplelist," --gname ",gnamelist," --diffcluster ",diffcluster," --group ",group," --minpct ",minpct_u," --logfcthreshold 0"))
	if (length(list.files(outdir)[grepl("diffexpressed.xls", list.files(outdir))])==0){
		print('unrecognized files provided, please check your input')
		quit()
	}else{
		genedf <- read.table(paste0(outdir,"/",list.files(outdir)[grepl("diffexpressed.xls", list.files(outdir))]),sep='\t',header=T,row.names=1)
	}
}else{
	if(!file.exists(diff)){
		print('unrecognized files provided, please check your input')
		quit()
	}else{
		genedf <- read.table(diff,sep="\t",header=T,row.names=1)
		if ("pct.1" %in% colnames(genedf)){ #兼容bulk数据
			genedf <- genedf[(genedf["pct.1"] >= minpct_u) | (genedf["pct.2"] >= minpct_u),]
			loc <- grep("avg|log", colnames(genedf))
			genedf <- genedf[genedf[,loc] != 0,]
		}
	}
}
loc <- grep("avg|log", colnames(genedf))
geneList = genedf[, loc]
names(geneList)= rownames(genedf)
geneList = sort(geneList,decreasing = T)
if (msigdb != "F"){
    geneset <- msigdbr(species = species, category = category, subcategory = subcategory)  %>% 
    dplyr::select(gs_name,gene_symbol)
	if (subpathway != "NULL" ){
    subpathway <- read.table(subpathway)
    geneset <- geneset[geneset$gs_name %in% subpathway$V1,]
	}
}else{
    geneset <- read.gmt(gmt)
}
#start GSEA analysis
if (saverds && file.exists(paste0(outdir,"/",prefix,'_','GSEA_enrichment_result.rds'))){
    egmt = readRDS(paste0(outdir,"/",prefix,'_','GSEA_enrichment_result.rds'))
}else{
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F,pvalueCutoff = 1, minGSSize = minGSSize, maxGSSize = maxGSSize )
if (saverds) {saveRDS(egmt,file=paste0(outdir,"/",prefix,'_','GSEA_enrichment_result.rds'))}
y=data.frame(egmt)
write.table(y,paste0(path,'/',prefix,'_','GSEA_enrichment.xls',sep=''),sep="\t",quote=F,row.names = FALSE)
y_sig <- filter(y, get(cutoff) < threshold)
write.table(y_sig,paste0(path,'/',prefix,'_','GSEA_enrichment_sig.xls',sep=''),sep="\t",quote=F,row.names = FALSE)
}
y=data.frame(egmt)
if ( is.na(selepath) ) {y_sig <- filter(y, get(cutoff) < threshold) }else{ y_sig <- y}
topshow = get_term(egmt, topn, selepath)

if(nrow(y_sig)<1){ 
    print('no signifcant path to plot GSEA enrichment plot ') 
    }else {
	print(topshow)
	enrich_plots(egmt, topshow)
	print(y_sig$ID[topshow])
}
#start GSEA pathway compare plot
if (pathway != "F"){
	sourceScript('color_protocol.R')
	color<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
	colorall<- c(color_protocol, color)
	data<-as.data.frame(pathway)
	data$color<-colorall[1:length(pathway)]
	data<-data[order(data$pathway),]
	len = ceiling(max(nchar(pathway),na.rm=T)/10)
	if (len<3){len<-3}else{len<-len}
	print(len)
	}else{
		print("do not conduct pathway compare")
}
