library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
options( stringsAsFactors = F)

#1 新色板，样本数超过24个后，后面的样本使用旧色板
clustcol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B",
              "#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101",
              "#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B",
              "#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B",
              "#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060",
              "#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC",
              "#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8")

#2 自定义函数——对list进行筛选
selectTCR <- function(contig){
  
  #筛选指定列
  contigTRA <- contig[,c("TRAV", "TRAJ", "chain1", "sample", "group")]
  contigTRB <- contig[,c("TRBV", "TRBJ", "chain2", "sample", "group")]
  
  #对列重命名
  colnames(contigTRA) <- c("v_gene", "j_gene", "chain", "sample", "group")
  colnames(contigTRB) <- c("v_gene", "j_gene", "chain", "sample", "group")
  
  #按行合并
  contig <- rbind(contigTRA, contigTRB)
}

selectBCR <- function(contig){
  contigIGH <- contig[,c("IGHV", "IGHJ", "chain1", "sample", "group")]
  contigIGL <- contig[,c("IGLV", "IGLJ", "chain2", "sample", "group")]
  colnames(contigIGH) <- c("v_gene", "j_gene", "chain", "sample", "group")
  colnames(contigIGL) <- c("v_gene", "j_gene", "chain", "sample", "group")
  contig <- rbind(contigIGH, contigIGL)
}

dataProcess <- function(combined){
  
  if(vdjType == "TCR") {
    newcombined <- list()
    
    #使用lapply批量筛选sc.combinedSample（输出认为list）
    newcombined <- lapply(combined, selectTCR)
    
    names(newcombined) <- names(combined)
  }else{
    newcombined <- list()
    newcombined <- lapply(combined, selectBCR)
    names(newcombined) <- names(combined)
  }
  return(newcombined)
}

theme_bar <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.text.x = element_text(size = 12,angle=45,hjust=1,vjust=1,face="bold"),
          axis.title.x = element_blank(),
          axis.title.y=element_text(face = "bold",size = 14),
          axis.text = element_text(face = "bold",size = 12),
          axis.ticks.x = element_blank(),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          #legend.justification=c(1,0),
          legend.position = "right", 
          aspect.ratio = 0.6,         #指定横纵坐标轴比例
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=20)),
    )
}

#3 
plotFreqPerSample <- function(combined, top = 20,color = NA) {
  
  ##3.1 筛选指定列——按样本对sc.combinedSample筛选
  plotcombined <- dataProcess(combined)

  ##3.2 
  data <- data.frame()
  for(i  in 1:length(plotcombined)) {
    data <- rbind(data, plotcombined[[i]])
  }
  
  ##3.3 
  if (vdjType == "TCR"){
    type <- c("TRAV","TRAJ","TRBV","TRBJ")
    
    for (t in type){
      
      # 筛选A/B链
      if (grepl('TRA', t)){
        data.use <- data[data$chain == 'TRA',]
        gene <- ifelse(t == "TRAV", 'v_gene', 'j_gene')
      }else{
        data.use <- data[data$chain == 'TRB',]
        gene <- ifelse(t == "TRBV", 'v_gene', 'j_gene')
      }
      
      # 统计每个克隆型在每个样本中的频数
      table_s <- table(data.use[[gene]], data.use$sample)
    
      # 统计每个克隆型在所有样本中的频数
      sum <- rowSums(table_s)
      
      # 按每个克隆型在所有样本中的频数对table_s排序——###########也是根据克隆型在所有样本中频数最多的来筛选top克隆###########
      table_s <- table_s[names(sum[order(sum,decreasing = T)]),,drop = F]
      # 计算每个样本中克隆型的频率——###########展示时频率是按样本统计###########
      table.prop_s <- prop.table(x = table_s, margin=2)
      
      write.table(table_s, paste0(outdir, "/10.sc.TopNgenexpTop/", prefix, ".", subprefix, "." , t, '.FreqPerSample.xls'), sep = '\t',quote = F, row.names = T, col.names =  NA)
      write.table(table.prop_s, paste0(outdir, "/10.sc.TopNgenexpTop/", prefix, ".", subprefix, "." , t, '.PropPerSample.xls'), sep = '\t',quote = F, row.names = T, col.names =  NA)
      
      # 筛选作图数据——默认绘制top20，不足20全画
      data.top <- table.prop_s[1:min(top,nrow(table.prop_s)),,drop = F]
      
      # 宽变长
      data_m <- melt(data.top)
      colnames(data_m) <- c("Gene","sample","Percentage")
      
     	# 指定legend顺序（与输入相同）
      data_m$sample <- factor(data_m$sample, levels = levels_sp)

      #画图
      if(!is.na(color[1])){
		clustcol <- color
      }
      pb <- ggplot(data_m, mapping = aes(x = Gene, y = Percentage, fill = sample)) +
        geom_col(position = "dodge", width = 0.5) +
        scale_fill_manual(values = clustcol)+
        scale_y_continuous(expand = c(0, 0)) +
        #ylim(c(0,max(data_m$Percentage) + 0.05)) +
        theme_bar() +
        theme(legend.key.size = unit(0.5, "cm")) +  #图例大小
        guides(fill = guide_legend(ncol = 1))     #图例列数
      
      pdf(paste0(outdir, "/10.sc.TopNgenexpTop/", prefix, "." , subprefix, "." , t, '.SamplePercentage.barplot.pdf'), width = 9, height = 9)
      print(pb)
      dev.off()
      png(paste0(outdir, "/10.sc.TopNgenexpTop/", prefix, "." , subprefix, "." , t, '.SamplePercentage.barplot.png'), width = 900, height = 900)
      print(pb)
      dev.off()
    }
    
  } else if (vdjType == "BCR"){
    
    type <- c("IGHV","IGHJ","IGL_IGKV","IGL_IGKJ")
    
    for (t in type){
      
      # 筛选H/L/K链
      if (grepl('IGH', t)){
        data.use <- data[data$chain == 'IGH',]
        gene <- ifelse(t == "IGHV", 'v_gene', 'j_gene')
      }else{
        data.use <- data[data$chain != 'IGH',]
        gene <- ifelse(t == "IGL_IGKV", 'v_gene', 'j_gene')
      }
      
      #
      if(length(combined) != 1){
	   table_s <- table(data.use[[gene]], data.use$sample)
   	   
   	   #
   	   sum <- rowSums(table_s)
   	   
   	   #
   	   table_s <- table_s[names(sum[order(sum,decreasing = T)]),,drop = F]
   	   
   	   #
       table.prop_s <- prop.table(x = table_s,margin=2)
      
       write.table(table_s, paste0(outdir, "/10.sc.TopNgenexpTop/", prefix, ".", subprefix, "." , t, '.FreqPerSample.xls'), sep = '\t',quote = F, row.names = T, col.names =  NA)
       write.table(table.prop_s, paste0(outdir, "/10.sc.TopNgenexpTop/", prefix, ".", subprefix, "." , t, '.PropPerSample.xls'), sep = '\t',quote = F, row.names = T, col.names =  NA)
      
      #
       data.top <- table.prop_s[1:min(top,nrow(table.prop_s)),]
      
      #
       data_m <- melt(data.top)
       colnames(data_m) <- c("Gene","sample","Percentage")
	  }else{
		table_s <- table(data.use[[gene]], data.use$sample)
		sum <- as.data.frame(rowSums(table_s))
		table_s <- table_s[order(sum, decreasing = T),]
		data_m <- data.frame(table_s)
		data_m <- data.frame(Gene = rownames(data_m),sample = rep(names(combined),nrow(data_m)),Percentage = data_m$table_s/sum(data_m$table_s))
		data_m <- data_m[1:min(top,nrow(data_m)),]
		data_m <- data_m[order(data_m$Percentage,decreasing = T),]
		print(data_m)
		data_m$Gene <- factor(data_m$Gene,levels = data_m$Gene)
      }
      #指定legend顺序（与输入相同）	
      data_m$sample<-factor(data_m$sample, levels = levels_sp)
      #head(data_m)
      if(!is.na(color[1])){
		clustcol <- color	
	  }
      pb <- ggplot(data_m, mapping = aes(x = Gene, y = Percentage, fill = sample)) +
        geom_bar(position = position_dodge(0.7), width = 0.5, stat = "identity", color = "transparent", size = 0) +
        scale_fill_manual(values = clustcol)+
        scale_y_continuous(expand = c(0, 0)) +
        #ylim(c(0, max(data_m$Percentage) + 0.05)) +
        theme_bar() +
        theme(legend.key.size = unit(0.5, "cm")) +  #图例大小
        guides(fill = guide_legend(ncol = 1))     #图例列数
      
      pdf(paste0(outdir, "/10.sc.TopNgenexpTop/", prefix, ".", subprefix, "." , t, '.SamplePercentage.barplot.pdf'), width = 9 , height = 9)
      print(pb)
      dev.off()
      png(paste0(outdir, "/10.sc.TopNgenexpTop/", prefix, ".", subprefix, "." , t, '.SamplePercentage.barplot.png'), width = 900 , height = 900)
      print(pb)
      dev.off()
    }
  }
}

