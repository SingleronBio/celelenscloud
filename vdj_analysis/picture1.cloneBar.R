#1 按样本画图
cloneBarSample <- function(combined) {
  
  #1 数据处理
  data <- data.frame()
  dataFre <- data.frame()
  
  for (i in 1:length(combined)) {
    tem <- drop_na(combined[[i]], CTgene) %>%
      dplyr::select(sample, Frequency)
    
    tem$clonal_expansion <- as.character(cut(tem$Frequency, fresplit, fresplitname))
    
    tem <- group_by(tem, sample, clonal_expansion) %>%
      summarise(count = n())   
    
    #暂时不用
    data <- rbind(data, tem)
    
    tem <- dplyr::mutate(tem, fre = count/ sum(count))  #前面分组了，此处统计仍是按组统计
    dataFre <- rbind(dataFre, tem)
  }
  
  dataFre$clonal_expansion <- factor(dataFre$clonal_expansion, levels = fresplitname)
  dataFre$sample <- factor(dataFre$sample, levels = unique(samples))
  
  #2 画图
  color_ce <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B",
                "#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101",
                "#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B",
                "#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B",
                "#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060",
                "#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC",
                "#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8") 
  color_ce <- color_ce[1:length(fresplitname)]
  names(color_ce) <- fresplitname
  
  p <- ggplot(dataFre, aes(x = sample, y = fre)) +
    geom_col(aes(fill = clonal_expansion), position = 'stack', width = 0.7) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = color_ce) +
    theme_half_open() +
    labs(x = "", y = "Frequency") +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     size = 14),
          axis.text.y = element_text(size = 14)) +
    theme(legend.position = "right") +
    guides(fill=guide_legend(title = "")) +
    theme(aspect.ratio = 1)
  
  
  #a <- list()
  #a[[1]] <- dataFre
  #a[[2]] <- p
  #return(a)
  write.table(dataFre, file = str_c(outdir, "/01.sc.CloneBar/", prefix, ".sample_data.xls"), sep = '\t', quote = F, row.names = F, col.names = c("Sample", "Dscription", "count", "Frequency"))
  
  pdf(str_c(outdir, "/01.sc.CloneBar/", prefix, ".sample_data.pdf"), width = 8,height = 8)
  print(p)
  dev.off()
  
  png(str_c(outdir, "/01.sc.CloneBar/", prefix,".sample_data.png"), width = 800,height = 800)
  print(p)
  dev.off()
  
}

#2 按样本分组画图
cloneBarGroup <- function(combined) {
  
  #1 数据处理
  data <- data.frame()
  dataFre <- data.frame()
  
  for (i in 1:length(combined)) {
    tem <- drop_na(combined[[i]], CTgene) %>%
      dplyr::select(group, group_freq)
    
    tem$clonal_expansion <- as.character(cut(tem$group_freq, fresplit, fresplitname))
    
    tem <- group_by(tem, group, clonal_expansion) %>%
      summarise(count = n())
    
    #暂时不用
    data <- rbind(data, tem)
    
    tem <- dplyr::mutate(tem, fre = count/ sum(count))
    dataFre <- rbind(dataFre, tem)
  }
  
  dataFre$clonal_expansion <- factor(dataFre$clonal_expansion, levels = fresplitname)
  dataFre$group <- factor(dataFre$group, levels = unique(groups))
  
  #2 画图
  color_ce <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B",
                "#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101",
                "#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B",
                "#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B",
                "#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060",
                "#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC",
                "#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8") 
  color_ce <- color_ce[1:length(fresplitname)]
  names(color_ce) <- fresplitname
  
  p <- ggplot(dataFre, aes(x = group, y = fre)) +
    geom_col(aes(fill = clonal_expansion), position = 'stack', width = 0.7) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = color_ce) +
    theme_half_open() +
    labs(x = "", y = "Frequency") +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     size = 14),
          axis.text.y = element_text(size = 14)) +
    theme(legend.position = "right") +
    guides(fill=guide_legend(title = "")) +
    theme(aspect.ratio = 1)
  
  #a <- list()
  #a[[1]] <- dataFre
  #a[[2]] <- p
  #return(a)
  
  write.table(dataFre, file = str_c(outdir, "/01.sc.CloneBar/", prefix, ".group_data.xls"), sep = '\t', quote = F, row.names = F, col.names = c("Group", "Dscription", "count", "Frequency"))
  
  pdf(str_c(outdir, "/01.sc.CloneBar/", prefix, ".group_data.pdf"), width = 4,height = 8)
  print(p)
  dev.off()
  
  png(str_c(outdir, "/01.sc.CloneBar/", prefix,".group_data.png"), width = 400,height = 800)
  print(p)
  dev.off()
}

#####用来统计每个样本中匹配到VDJ信息的细胞占比
statClonetype_percentage <- function(combined){
	for(i in 1:length(names(combined))){
		sp_tmp <- names(combined)[i]
		tmp <- combined[[sp_tmp]]
		if(i == 1){
			out <- data.frame(sample = sp_tmp,ColonType_count = nrow(tmp[!is.na(tmp$CTgene),]),Cell_Number = nrow(tmp),ColonType_Freq = nrow(tmp[!is.na(tmp$CTgene),])/nrow(tmp))
		}else{
			out_tmp <- data.frame(sample = sp_tmp,ColonType_count = nrow(tmp[!is.na(tmp$CTgene),]),Cell_Number = nrow(tmp),ColonType_Freq = nrow(tmp[!is.na(tmp$CTgene),])/nrow(tmp))
			out <- rbind(out,out_tmp)
		}
	}
	write.table(out,paste0(outdir,"/00.QC/",prefix,".fraction_of_colontype.xls"),sep = "\t",quote = F,row.names = F)
}
