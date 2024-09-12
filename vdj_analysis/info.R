#
getInfo <- function(combined, split) {
  compare <- str_c(outdir,"/05.sc.Otherinfo/")
  
  #1 判断vdjType
  if (vdjType == "TCR") {
    selec_col <- c("CTgene", "TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ", "CTnt", "CTaa","CTstrict","CTstrictaa")
    
    col_name1 <- "VDJgeneA"
    col1 <- c("TRAV", "TRAJ")
    
    col_name2 <- "VDJgeneB"
    col2 <- c("TRBV", "TRBD", "TRBJ")
  }
  
  if (vdjType == "BCR") {
    selec_col <- c("CTgene", "IGLV", "IGLJ", "IGHV", "IGHD", "IGHJ", "CTnt", "CTaa","CTstrict","CTstrictaa")
    
    col_name1 <- "VDJgeneL"
    col1 <- c("IGLV", "IGLJ")
    
    col_name2 <- "VDJgeneH"
    col2 <- c("IGHV", "IGHD", "IGHJ")
  }
  
  #2 判断提取sanple/group统计频率
  if (split == "sample") {
    select_fre <- "Frequency"
  } else {
    select_fre <- "group_freq"
  }
  
  #3 判断cloneDefine
  if (cloneDefine == "VDJ") {
    define = "CTgene"
	sel_tmp <- c("CTaa","CTstrict")
  } else if(cloneDefine == "CDR3"){
    define = "CTaa"
	sel_tmp <- c("CTgene","CTstrict")
  } else if(cloneDefine == "VDJ_CDR3"){
	define = "CTstrict"
	sel_tmp <- c("CTgene","CTaa")
  } else if(cloneDefine == "VDJ_CDR3aa"){
	define = "CTstrictaa"
	sel_tmp <- c("CTgene")
  }
  
  #3
  if(length(combined) <= 7){
    venn_use <- list()
	plot_venn <- TRUE
  }else{
    plot_venn <- FALSE
  }
  for (i in 1:length(combined)) {
    a <- dplyr::select(combined[[i]], c(selec_col, get("select_fre"))) %>%
      dplyr::rename(clonetype = !!sym(ifelse(cloneDefine == "VDJ","CTgene",ifelse(cloneDefine == "CDR3","CTaa",ifelse(cloneDefine == "VDJ_CDR3","CTstrict","CTstrictaa")))), Freq = get("select_fre")) %>%
      unite(!!sym(col_name1), get("col1")) %>%
      unite(!!sym(col_name2), get("col2")) %>%
      dplyr::select(clonetype, Freq, get("col_name1"), get("col_name2"), CTnt, sel_tmp) %>%
      tidyr::drop_na() %>%
      arrange(desc(Freq))
#      dplyr::distinct(!!sym("clonetype"), .keep_all = TRUE)
#	   这里的distinct只保留了相同clonetype的第一列，有时候基于CTgene定义，CDR3会有不同的，故后边若有不同的CDR3序列，则paste到一起输出
	###将相同clonetype下，有不同CDR3序列的信息paste到一起输出
	clt <- unique(a$clonetype)
	clt_output <- NULL
	for(clonotype in clt){
		clt_use_tmp <- a[a$clonetype == clonotype,]
		if(length(unique(clt_use_tmp$CTnt)) == 1){
			clt_output <- rbind(clt_output,clt_use_tmp[1,])
		}else{
			clt_tmp <- data.frame(clonetype = clonotype,Freq = unique(clt_use_tmp$Freq))
			clt_tmp[,col_name1] <- unique(clt_use_tmp[,col_name1])
			clt_tmp[,col_name2] <- unique(clt_use_tmp[,col_name2])
			clt_tmp[,"CTnt"] <- paste(unique(clt_use_tmp$CTnt),collapse = ",")
			for(col_name in sel_tmp){
				clt_tmp[,col_name] <- paste(unique(clt_use_tmp[,col_name]),collapse = ",")
			}
			clt_output <- rbind(clt_output,clt_tmp)
		}
	}
	a <- clt_output
	####
    row_name <- seq(dim(a)[1])
    row_name <- str_c("Clonetype", row_name)
    rownames(a) <- row_name
    
    a <- rownames_to_column(a, var = "name")
    write.table(a, file=str_c(compare, "/", names(combined)[i], ".filter.cloneFreq.xls"), sep='\t', quote=F, row.names=F)
	if(plot_venn){
		venn_use[[names(combined)[i]]] <- a[,"clonetype"]
	}
  }
  if(plot_venn){
    pdf(paste0(compare, "/Venn.filter.",split,".pdf"),width = 8)
    venn(venn_use,ilabels = "counts",zcolor = color_protocol)
	dev.off()
	png(paste0(compare, "/Venn.filter.",split,".png"),width = 8/7*480)
    venn(venn_use,ilabels = "counts",zcolor = color_protocol)
    dev.off()
	inter <- attributes(venn(venn_use,ilabels = "counts",zcolor = color_protocol))
	inter <- inter$intersections
	sink(paste0(compare, "/Venn.filter.",split,".txt"))
	for(i in names(inter)){
		if(grepl("[:]",i)){
			cat(paste0(i,"\n"))
			cat(paste(inter[[i]],collapse = "\n"))
			cat("\n")
			cat("\n")
		}
	}
	sink()
  }
}


####提取rds中颜色(from fangsheng)
getColorsMeta <- function(rds,meta){
    metadata <- rds@meta.data
    meta_colors <- paste0(meta,"_colors")
    if (meta_colors %in% colnames(metadata)){
        data <- metadata[,c(meta,meta_colors)]
        data_colors <- data[,2]
        names(data_colors) <- data[,1]
        data_colors <- data_colors[!duplicated(data_colors)]
        if(!is.null(levels(data[,1]))){
            data_colors <- data_colors[levels(data[,1])]
        }
        return(data_colors)
    }else {
#        print(paste0("Error: ",meta_colors," is not in metadata,please check rds metadata column!"))
        return(NA)
    }
}

getColors <- function(rds,meta=c("sample","group","cluster","raw_cluster")){
    meta_colors <- list()
    for (i in meta){
        k <- paste0(i,"_colors")
        meta_colors[[k]] <- getColorsMeta(rds,meta=i)
    }
    return(meta_colors)
}




