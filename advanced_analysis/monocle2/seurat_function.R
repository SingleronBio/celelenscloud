#' @name getColorsMeta
#' @param rds seurat object
#' @param meta choose metadata column for get colors
#' @return Return a vector for colors,this vector of names is metadata column
#' @export
#'
#' @examples
#'
#' rds <- readRDS("/SGRNJ03/PiplineTest01/Pipline_test/huangwanxiang/scanpy_test/batch1/P21050401.diff_PRO.rds")
#' color_cluster <- getcolors_meta(rds,meta="cluster")
#'
getColorsMeta <- function(rds,meta){
    metadata <- rds@meta.data   
    meta_colors <- paste0(meta,"_colors")
    if (meta_colors %in% colnames(metadata)){
        data <- metadata[,c(meta,meta_colors)]
        data <- unique(data)
        data_colors <- data[,2]
        names(data_colors) <- data[,1]
        if(!is.null(levels(data[,1]))){
            data_colors <- data_colors[levels(data[,1])]
        }
        return(data_colors)
    }else {
        print(paste0("Error: ",meta_colors," is not in metadata,please check rds metadata column!"))
        return(NULL)
    }
}

#---------------------------------------------------------------------------------------------------------------
#' @name getColors
#' @param rds seurat object
#' @param meta choose metadata column for get colors,default:c("sample","group","cluster","raw_cluster")
#' @return Return a list for colors,this vector of names is metadata column
#' @export
#'
#' @examples
#'
#' rds <- readRDS("/SGRNJ03/PiplineTest01/Pipline_test/huangwanxiang/scanpy_test/batch1/P21050401.diff_PRO.rds")
#' colorlist <- getColors(rds,meta=c("sample","group","cluster","raw_cluster"))
#'
getColors <- function(rds,meta=c("sample","group","cluster","raw_cluster")){
    meta_colors <- list()
    for (i in meta){
        k <- paste0(i,"_colors")
        meta_colors[[k]] <- getColorsMeta(rds,meta=i)
    }
    return(meta_colors)
}

#---------------------------------------------------------------------------------------------------------------
#' @name setColorsMeta
#' @param rds seurat object
#' @param meta choose metadata column for set colors
#' @param value meta names 
#' @param color color for meta names 
#' @return Return seurat object
#' @export
#'
#' @examples
#'
#' rds <- readRDS("/SGRNJ03/PiplineTest01/Pipline_test/huangwanxiang/scanpy_test/batch1/P21050401.diff_PRO.rds")
#' group <- c("A","B","C")
#' color_group <- c("#59a4ce","#ff8ab6","#8a5626")
#' rds1 <- setColorsMeta(rds,meta="group",value=group,color=color_group)
#'

setColorsMeta <- function(rds,meta,value,color){
    metadata <- rds@meta.data 
    meta_colors <- paste0(meta,"_colors")
    names(color) <- value
    if (! meta_colors %in% colnames(metadata)){
        metadata[,meta_colors] <- ""
    }
    if (meta %in% colnames(metadata)){
        if (! meta_colors %in%  colnames(metadata)){
            a <- paste0(colnames(metadata,meta_colors))
            metadata$newdaa <- ""
            colnames(metadata) <- a
        }
    #    intercolor <-  intersect(color,unique(as.character(metadata[,meta_colors])))
        intercolor <-  intersect(color,unique(as.character(metadata[which(! metadata[,meta] %in% value),][,meta_colors])))
        if ( length(intercolor) >0    ){
            stop(paste0("function(setColorsMeta) Error: ",intercolor," is duplicated in " ,meta_colors,", please set new colors"))
        }
        for (i in value){
            if (i %in% unique(as.character(metadata[,meta]))){
                metadata[which(metadata[,meta] == i),][,meta_colors] <- color[i]
            } else {
                stop(paste0("function(setColorsMeta) Error:",i," is not in ",meta,",please check rds metadata column!"))
            }
        }
    metadata[,meta_colors] <- as.array(metadata[,meta_colors])
    rds <- AddMetaData(rds,metadata)
    print(paste0("change color of ",meta))
    print(getColors(rds,meta=meta))
    return(rds)
    }else {
        stop(paste0("function(setColorsMeta) Error: ",meta," is not in metadata,please check rds metadata column!"))
    }
}

#---------------------------------------------------------------------------------------------------------------
#' @name getPlotParams
#' @param rds seurat object
#' @param meta choose metadata column for get some plotparams,can choose "sample","group","cluster","raw_cluster","color","pointSize"
#' @return Return a list
#' @export
#'
#' @examples
#' rds <- readRDS("/SGRNJ03/PiplineTest01/Pipline_test/huangwanxiang/scanpy_test/batch1/P21050401.diff_PRO.rds")
#' colorlist <- getPlotParams(rds,meta="color")
#' pointsizelist <- getPlotParams(rds,meta="pointSize")
#' samplelist <- getPlotParams(rds,meta="sample")
#' grouplist <- getPlotParams(rds,meta="group")
#' clusterlist <- getPlotParams(rds,meta="cluster")
#'
getPlotParams <- function(rds, meta){
    params <- list()
    if (meta %in% colnames(rds@meta.data)){
        group <- as.character(unique(rds@meta.data[,meta]))
        gnumber<-length(group)
        params[[paste0(meta,"_mar_bottom")]] <- max(nchar(group))/2 + 2.5
        cluster <- levels(rds)
        if(!(NA %in% (as.numeric(cluster)))){
            params[[paste0(meta,"_cex")]] <- ifelse(gnumber>10, 1-(gnumber-10)*0.02, 1)
            params[[paste0(meta,"_mix")]] <- ifelse(gnumber>15, 1, 15/gnumber)
            params[[paste0(meta,"_ncol")]] <- 3
        }else{
            params[[paste0(meta,"_cex")]] <- ifelse(gnumber>20, 1-(gnumber-20)*0.02, 1)
            params[[paste0(meta,"_mix")]] <- ifelse(30/gnumber<1.5, 1.5, 30/gnumber)
            params[[paste0(meta,"_ncol")]] <- ceiling(gnumber/15)
        }
    } else if (meta == "color"){
        params <- getColors(rds,meta=c("sample","group","cluster","raw_cluster"))
    } else if (meta == 'pointSize'){
        cell_number<-length(colnames(rds))
        params[["pt_use"]] <- ifelse(cell_number > 6500, 0.1, ifelse(cell_number > 5500, 0.15, ifelse(cell_number > 4000, 0.2, ifelse(cell_number > 2500, 0.3, ifelse(cell_number > 1000, 0.4, 0.6)))))
    } else{
        stop(paste0("Error: ",meta," is not in metadata,please check rds metadata column!"))
    }
    return(params)
}

#---------------------------------------------------------------------------------------------------------------
#' @name statCompair_base
#' @param meta metadata dataframe
#' @param valuename Select the column name to be calculated
#' @param groupname Select column names for comparison
#' @param cmethod the type of test.you can chooese 't.test','anova','wilcox.test','kruskal.test'.Default is wilcox.test
#' @param p.adjust.method choose pvalue adjust method("holm","bonferroni","hochberg","hommel","BH","BY","fdr"),default is "bonferroni"
#' @param paired a logical indicating whether you want a paired test. Used only in 't.test' and in wilcox.test
#' @param alternative Select inspection type("two.sided","less","greater"),default is "two.sided"
#' @return Return a list
#' @export
#'
#' @examples
statCompair_base <- function(meta,valuename,groupname,cmethod="wilcox.test",p.adjust.method="bonferroni",paired = FALSE,alternative ="two.sided" ){
if (valuename %in% colnames(meta) | groupname %in% colnames(meta)){
        met1 <- meta[,c(valuename,groupname)]
        colnames(met1) <- c("value","compare")
        stat_test<- compare_means(value~compare,data = met1,method = cmethod, p.adjust.method = p.adjust.method,paired=paired,alternative = alternative)
        stat_test[,1] <- valuename
        stat_test$pair <- as.character(paired)
        stat_test$tail <- alternative
        stat_test <- as.data.frame(stat_test)
        df.mean=group_by(met1, compare) %>% summarize_each(funs(mean))
        df.med=group_by(met1,  compare) %>% summarize_each(funs(median))
        df.mean <- as.data.frame(df.mean)
        colnames(df.mean) <- c("compare","mean")
        df.med <- as.data.frame(df.med)
        colnames(df.med) <- c("compare","median")
        df.stat <- merge(df.mean,df.med,by="compare")
        df.stat$y <- valuename
    } else {
        stop(paste0(valuename," or ",groupname," is not in metadata, please check it !"))
    }
    return(list(compair=stat_test,stat=df.stat))
}


#---------------------------------------------------------------------------------------------------------------
#' @name statCompair
#' @param meta metadata dataframe
#' @param valuename Select the column name to be calculated
#' @param groupname Select column names for comparison
#' @param compare Select comparison relationship,default is "all" for test all comparison relationship,you can fill in this example("P-nM-4vsP-nM-5,P-nM-4vsP-nM-6")
#' @param cmethod the type of test.you can chooese 't.test','anova','wilcox.test','kruskal.test'.Default is wilcox.test
#' @param p.adjust.method choose pvalue adjust method("holm","bonferroni","hochberg","hommel","BH","BY","fdr"),default is "bonferroni"
#' @param paired a logical indicating whether you want a paired test. Used only in 't.test' and in wilcox.test
#' @param alternative Select inspection type("two.sided","less","greater"),default is "two.sided"
#' @return Return a list
#' @export
#'
#' @examples

statCompair <- function(meta,valuename,groupname,compare="all",cmethod="wilcox.test",p.adjust.method="bonferroni",paired = FALSE,alternative ="two.sided"){
    if(compare!="all"){
        compare <- unlist(strsplit(compare,split=","))
        for(i in 1:length(compare)){
            com <- compare[i]
            com <- unlist(strsplit(com,split="vs"))
            if (groupname %in% colnames(meta)){
                sub_meta <- meta[meta[,groupname] %in% com,]
                sub_stat <- statCompair_base(meta=sub_meta,valuename=valuename,groupname=groupname,cmethod=cmethod,p.adjust.method=p.adjust.method,paire=paired,alternative = alternative)
                if (i == 1){
                    stat <- sub_stat
                } else {
                    stat[["compair"]] <- rbind(stat[["compair"]],sub_stat[["compair"]])
                    stat[["stat"]] <- rbind(stat[["stat"]],sub_stat[["stat"]])
                    
                }
                stat[["compair"]] <- as.data.frame(stat[["compair"]])
                stat[["stat"]] <- as.data.frame(stat[["stat"]] %>% group_by(compare, mean, median) %>% slice_sample(n = 1))
            } else {
               stop(paste0(groupname," is not in metadata, please check it !")) 
            }
        }
    }else{
        stat <- statCompair_base(meta=meta,valuename=valuename,groupname=groupname,cmethod=cmethod,p.adjust.method=p.adjust.method,paired=paired,alternative = alternative)
    }
    return(stat)
}
