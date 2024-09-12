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

#3 画图
plotVJHeatmap <- function(combined) {
  
  ##3.1 筛选指定列——按样本对sc.combinedSample筛选
  plotcombined <- dataProcess(combined)
  
  ##
  if (vdjType == "TCR"){
   
    for (i in 1:length(plotcombined)){
      
      # 剃掉v_gene/j_gene含NA的行
      data.use <- drop_na(plotcombined[[i]], v_gene, j_gene)
      data.use <- dplyr::filter(plotcombined[[i]], v_gene != "None" & j_gene != "None")

      # 统计克隆型之间的组装频数
      mtx <- table(data.use$v_gene, data.use$j_gene)
      
      # 统计克隆型之间的组装频率——每对克隆型频数/总克隆型频数
      mtx <- mtx/nrow(data.use)*100
      mtx <- matrix(mtx,ncol=ncol(mtx),dimnames = dimnames(mtx))
      print(class(mtx))
      # 画图
      write.table(t(mtx), file=paste0(outdir, "/11.sc.VJpairs/", prefix, ".", names(plotcombined)[i], ".VJ_pairs.heatmap.csv"), sep=",", quote = FALSE, row.names = TRUE)
      ph <- draw_ComplexHeatmap(t(mtx),scale_or_not = FALSE, heat_col = colorRampPalette(c("#FFFAFA","navy"))(50),cellwidth = 10,cellheight = 10, cluster_rows = FALSE, cluster_columns = FALSE)
      save_heatmap(ph,pdf = TRUE, png = TRUE, prefix = paste0(prefix, ".", names(plotcombined)[i], '.VJ_pairs'),outdir = paste0(outdir,"/11.sc.VJpairs/"))
#      pdf(paste0(outdir, "/11.sc.VJpairs/", prefix, ".", names(plotcombined)[i], '.VJ_pairs.heatmap.pdf'), width = 15 , height = 10)
#      ph <- pheatmap(t(mtx),
#                     border_color = "white",
#                     breaks=seq(0, max(mtx), length.out = 50),
#                     col = colorRampPalette(c("#FFFAFA","navy"))(50),
#                     cluster_rows = F,
#                     cluster_cols = F,
#                     main = names(plotcombined)[i])
#      dev.off()
      
#      png(paste0(outdir, "/11.sc.VJpairs/", prefix, ".", names(plotcombined)[i], '.VJ_pairs.heatmap.png'), width = 1500 , height = 1000)
#      ph <- pheatmap(t(mtx),
#                     border_color = "white",
#                     breaks=seq(0, max(mtx), length.out = 50),
#                     col = colorRampPalette(c("#FFFAFA","navy"))(50),
#                     cluster_rows = F, 
#                     cluster_cols = F,
#                     main = names(plotcombined)[i])
#      dev.off()

    }
  } else {
   
    for (i in 1:length(plotcombined)){

      # 剃掉v_gene/j_gene含NA的行
      data.use <- drop_na(plotcombined[[i]], v_gene, j_gene)
      data.use <- dplyr::filter(data.use, v_gene != "None" & j_gene != "None")

      #
      mtx <- as.matrix(table(data.use$v_gene, data.use$j_gene))
      #
      mtx <- mtx/nrow(data.use)*100
      #
      mtx <- as.matrix(t(mtx))
      mtx <- matrix(mtx,ncol=ncol(mtx),dimnames = dimnames(mtx))
	  print(class(mtx))
      write.table(mtx, file=paste0(outdir, "/11.sc.VJpairs/", prefix, ".", names(plotcombined)[i], ".VJ_pairs.heatmap.csv"), sep=",", quote = FALSE, row.names = TRUE)
      ph <- draw_ComplexHeatmap(mtx,scale_or_not = FALSE,heat_col = colorRampPalette(c("#FFFAFA","navy"))(50),cellwidth = 10,cellheight = 10, cluster_rows = FALSE, cluster_columns = FALSE)
	 # colorRampPalette(c("#FFFAFA","navy"))(50)
      save_heatmap(ph,pdf = TRUE, png = TRUE, prefix = paste0(prefix, ".", names(plotcombined)[i], '.VJ_pairs'),outdir = paste0(outdir,"/11.sc.VJpairs/"))
      
#      pdf(paste0(outdir, "/11.sc.VJpairs/", prefix, ".", names(plotcombined)[i], '.VJ_pairs.heatmap.pdf'), width = 15 , height = 10)
#      ph <- pheatmap(t(mtx),
#                     border_color = "white",
#                     breaks=seq(0, max(mtx), length.out = 50),
#                     col = colorRampPalette(c("#FFFAFA","navy"))(50),
#                     cluster_rows = F,
#                     cluster_cols = F,
#                     main = names(plotcombined)[i])
#      dev.off()
#      
#      png(paste0(outdir, "/11.sc.VJpairs/", prefix, ".", names(plotcombined)[i], '.VJ_pairs.heatmap.png'), width = 1500 , height = 1000)
#      ph <- pheatmap(t(mtx),
#                     border_color = "white",
#                     breaks=seq(0, max(mtx), length.out = 50),
#                     col = colorRampPalette(c("#FFFAFA","navy"))(50),
#                     cluster_rows = F, 
#                     cluster_cols = F,
#                     main = names(plotcombined)[i])
#      dev.off()
    }
  }
}
