#1 按样本画图
haveOrNotsample <- function(combined) {
  
  #1 判断vdjType
  if (vdjType == "TCR") {
    p_lable = "No_TCR"
  } else {
    p_lable = "No_BCR"
  }
  
  #2 数据处理
  for (i in 1:length(combined)) {
    #克隆型频率那一列为integer，需要
    combined[[i]]$Frequency <- as.numeric(combined[[i]]$Frequency)
    
    data <- dplyr::select(combined[[i]], sample, cluster, Frequency)
    
    data$clonal_expansion <- as.character(cut(data$Frequency, fresplit, fresplitname))
    data$clonal_expansion[is.na(data$clonal_expansion)] <- get("p_lable")
    
    data <- group_by(data, cluster, clonal_expansion) %>%
      summarise(count = n()) %>%
      as.data.frame() %>%
      dplyr::mutate(fre = count/sum(count))
    
    data$clonal_expansion <- factor(data$clonal_expansion, levels = c(get("p_lable"), fresplitname))
    data$cluster <- factor(data$cluster, levels = cluster)
    
    #3 定义颜色
    color_provide <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B",
                       "#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101",
                       "#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B",
                       "#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B",
                       "#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060",
                       "#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC",
                       "#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8") 
    
    color_ce =  c("#CED0CD", color_provide[1:length(fresplitname)])   
    names(color_ce) = c(get("p_lable"), fresplitname)
    
    #4 画图
    p <- ggplot(data, aes(x = cluster, y = fre)) +
      geom_col(aes(fill = clonal_expansion), position = 'stack', width = 0.3) +
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
    
    write.table(data, file = str_c(outdir,"/04.sc.HaveOrNot/", names(combined)[i], ".sample_data.xls"), sep = '\t', quote = F, row.names = F, col.names = c("cluster", "Dscription", "count", "Frequency"))
    
    pdf(str_c(outdir,"/04.sc.HaveOrNot/", names(combined)[i], ".sample_data.pdf"), width = 8,height = 8)
    print(p)
    dev.off()
    
    png(str_c(outdir,"/04.sc.HaveOrNot/", names(combined)[i], ".sample_data.png"), width = 800,height = 800)
    print(p)
    dev.off()
  }
}


#2 按样本分组画图
haveOrNotGroup <- function(combined) {
  
  #1 判断vdjType
  if (vdjType == "TCR") {
    p_lable = "No_TCR"
  } else {
    p_lable = "No_BCR"
  }
  
  #2 数据处理
  for (i in 1:length(combined)) {
    
    #将样本分组相关数据列名替换成样本相关列名
    combined[[i]]$sample <- combined[[i]]$group
    combined[[i]]$Frequency <- combined[[i]]$group_freq
    
    #克隆型频率那一列为integer，需要改为numeric
    combined[[i]]$Frequency <- as.numeric(combined[[i]]$Frequency)
    
    data <- dplyr::select(combined[[i]], sample, cluster, Frequency)
    
    data$clonal_expansion <- as.character(cut(data$Frequency, fresplit, fresplitname))
    data$clonal_expansion[is.na(data$clonal_expansion)] <- get("p_lable")
    
    data <- group_by(data, cluster, clonal_expansion) %>%
      summarise(count = n()) %>%
      as.data.frame() %>%
      dplyr::mutate(fre = count/sum(count))
    
    data$clonal_expansion <- factor(data$clonal_expansion, levels = c(get("p_lable"), fresplitname))
    data$cluster <- factor(data$cluster, levels = cluster)
    
    #3 定义颜色
    color_provide <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B",
                       "#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101",
                       "#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B",
                       "#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B",
                       "#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060",
                       "#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC",
                       "#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8")
    
    color_ce =  c("#CED0CD", color_provide[1:length(fresplitname)])   
    names(color_ce) = c(get("p_lable"), fresplitname)
    
    #4 画图
    p <- ggplot(data, aes(x = cluster, y = fre)) +
      geom_col(aes(fill = clonal_expansion), position = 'stack', width = 0.3) +
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
    
    write.table(data, file = str_c(outdir,"/04.sc.HaveOrNot/", names(combined)[i], ".group_data.xls"), sep = '\t', quote = F, row.names = F, col.names = c("cluster", "Dscription", "count", "Frequency"))
    
    pdf(str_c(outdir,"/04.sc.HaveOrNot/", names(combined)[i], ".group_data.pdf"), width = 8,height = 8)
    print(p)
    dev.off()
    
    png(str_c(outdir,"/04.sc.HaveOrNot/", names(combined)[i], ".group_data.png"), width = 800,height = 800)
    print(p)
    dev.off()
  }
}

