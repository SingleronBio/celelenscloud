library(SeuratObject)
library(tidyverse)
library(scatterpie)
library(RColorBrewer)
library(scales)
library(cowplot)

#绘制代码
scClonetype <- function(rds, sample_p, group_p){
  
  #1自定义函数，绘制legend时调用
  ##1.1 自定义画圈函数
  geom_circle <- function(r, xc, yc, color="black", fill=NA, ...) {
    x <- xc + r*cos(seq(0, pi, length.out=100))
    ymax <- yc + r*sin(seq(0, pi, length.out=100))
    ymin <- yc + r*sin(seq(0, -pi, length.out=100))
    annotate("ribbon", x=x, ymin=ymin, ymax=ymax, color=color, fill=fill, ...)
  }
  
  ##1.2 自定义函数通过半径回溯原始值
  labeller = function(x) {
    return(ceiling((2^(x*dot_size)-2)))
  }
  
  ##1.3 自定义函数，将lable向量中偶数位置转变成负数
  lableReverse <- function(x){
    for (i in 1:length(x)) {
      if (i %% 2 == 0) {
        x[i] <- -x[i]
      } else {
        x[i] <- x[i]
      }
    }
    
    return(x)
  }
  
  #2 读取rds
  pro <- rds
  
  #3 转换sample/cluster为字符型
  pro@meta.data$sample <- as.character(pro@meta.data$sample)
  pro@meta.data$cluster <- as.character(pro@meta.data$cluster)
  
  #4 从meta.data中提取克隆型频率
  df3 = pro@meta.data %>%
    rownames_to_column(var = "cell") %>%
    dplyr::select(sample, cluster, Frequency) %>%
    dplyr::filter(Frequency != "NA") %>%
    dplyr::mutate(celltype = cluster)
  
  df3$clonal_expansion <- as.character(cut(df3$Frequency, fresplit, fresplitname))
  
  df3 <- dplyr::select(df3, -Frequency, -cluster) %>%
    group_by(sample, celltype, clonal_expansion) %>%
    summarise(count = n()) %>%
    spread(clonal_expansion, count, fill = 0) %>%
    as.data.frame()
  
  #5 有些数据不好，可能到这步df3中缺少频数分组的某一列
  for (i in 1:length(fresplitname)) {
    if(any(str_detect(colnames(df3), fresplitname[i]))){
      print(str_c("the column of ", fresplitname[i], " in data frame"))
    } else {
      df3[[fresplitname[i]]] <- "0"
    }
  }
  
  # 将频数分组列改成数值型
  for(col_n in 3:length(colnames(df3))){
    df3[[col_n]] = as.numeric(df3[[col_n]])
  }
  
  #6 分别按sample及celltype提取统计数据
  ##6.1 分样本不分细胞类型统计
  df_all_celltype = data.frame()
  for (samplei in unique(df3$sample)) {
    tmp1 = df3 %>% filter(sample == samplei)
    b = c(sample = samplei,celltype = "alltype",colSums(tmp1[,3:dim(tmp1)[2]])) %>% as.data.frame() %>% t() %>% as.data.frame(stringsAsFactors=FALSE)
    df_all_celltype = rbind(df_all_celltype,b) #rbind向数据框添加行
  }
  
  ##6.2 分细胞类型不分样本统计
  df_all_sample = data.frame()
  for (typei in unique(df3$celltype)) {
    tmp2 = df3 %>% filter(celltype == typei)
    c = c(sample = "allsample",celltype = typei,colSums(tmp2[,3:dim(tmp2)[2]])) %>% as.data.frame() %>% t() %>% as.data.frame(stringsAsFactors=FALSE)
    df_all_sample=rbind(df_all_sample,c)
  }
  
  ##6.3 样本与细胞类型都不分统计
  df_all_sample_celltype = data.frame()
  d = c(sample="allsample",celltype="alltype",colSums(df3[,3:dim(df3)[2]])) %>% as.data.frame() %>% t() %>% as.data.frame(stringsAsFactors=FALSE)
  df_all_sample_celltype = d
  
  ##2.4 对数据框df_all_celltype/df_all_sample列进行命名
  #colnames(df_all_celltype) <- colnames(df3)
  #colnames(df_all_sample) <- colnames(df3)
  
  #6.4 合并提取的数据框
  df4 = rbind(df_all_celltype,df_all_sample) %>% rbind(df_all_sample_celltype)
  
  ##6.5 上述代码操作得到的df4s所有列的数据类型均是字符型，需要将除前两列以外的列转换成数值型，前两列为字符型
  for(col_n in 3:length(colnames(df4))){
    df4[[col_n]]=as.numeric(df4[[col_n]]) #此处不能使用$及数据框列名对数据框列进行指定，只能使用[[]]及数据框列数来指定数据框的列（防止需要转
    #换的列列名不规范报错）
  }
  
  df4$sample = as.character(df4$sample)
  df4$celltype = as.character(df4$celltype)
  
  #7 将df4追加在df3上,并将前两列转换成指定数据类型
  df5 = df3 %>% rbind(df4) %>%
    rownames_to_column(var = "tmp") %>%  #重命名行名
    dplyr::select(-tmp)
  
  df5$sample=as.character(df5$sample)
  df5$celltype=as.character(df5$celltype)
  
  #8 输出数据
  ##8.1 按组展示数据
  if(group_p == "T") {
    df5_save <- dplyr::rename(df5, "group" = "sample")
    write.table(df5_save,file=str_c(outdir,"/06.sc.CloneTypePie/",prefix,".sc.CloneType.gorup.Percent.xls"),sep="\t",quote=F,row.names=T,col.names = NA)
  }
  
  ##8.2 按样本展示数据
  if(sample_p == "T") {
    write.table(df5,file=str_c(outdir,"/06.sc.CloneTypePie/",prefix,".sc.CloneType.sample.Percent.xls"),sep="\t",quote=F,row.names=T,col.names = NA)
  }
  
  #9 转换数据列数据类型
  ##9.1 使样本/group顺序与输入一致
  
  ###9.1.1 按样本展示
  #sample <- unique(df5$sample)
  if (sample_p == "T") {
    #sample <- c(samples[!(samples %in% dplyr::setdiff(samples, unique(pro@meta.data$sample)))], "allsample")
    sample <- c(samples, "allsample")
  } else {
    print("picture group")
  }
  
  ###9.1.2 按样本分组展示
  if (group_p == "T"){
    sample <- c(unique(groups), "allsample")
  } else {
    print("picture sample")
  }
  
  ##9.2 将数据框df5的sample列转换成数字（负数）——match函数根据元素，输出元素位置
  for(sample_id in sample){
    df5$sample[df5$sample==sample_id] <- -match(sample_id,sample)
  }
  
  ##9.2 将数据框df5的celltype列转换成数字（正数）
  celltype <- unique(df5$celltype)
  
  for(celltype_id in celltype){
    df5$celltype[df5$celltype==celltype_id] <- match(celltype_id,celltype)
  }
  
  ##9.3 经过上述处理df5前两列虽为数字，但其仍是为字符型，需要转换成数值型
  for(col_n in 1:length(colnames(df5))){
    df5[[col_n]]=as.numeric(df5[[col_n]])
  }
  
  #10 给数据框新增一列，用于展示扇形图大小
  ##10.1 自动选取合适的dot_size，确保饼图大小合适
  dot_size <- ceiling(log2(max(rowSums(df5[,3:length(colnames(df5))])) + 2) / 0.5) + 2
  
  ##10.2 确定饼图半径（<0.5）
  df5$radius = log2(rowSums(df5[,3:length(colnames(df5))]) + 2) / dot_size
  
  #11 从色板中获取三种颜色，并命名
  ##11.1 旧色板
  #color_ce=brewer.pal(9, "YlGnBu")[c(3,5,8)]   
  #names(color_ce)=c("Single","Medium","Large")
  
  ##11.2 新色板
  color_provide <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B",
                "#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101",
                "#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B",
                "#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B",
                "#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060",
                "#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC",
                "#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8")  
  color_ce <- color_provide[1:length(fresplitname)]
  names(color_ce) <- fresplitname
  
  #12 绘图
  ##12.1 不带圆圈图形
  p1 <- ggplot() +
    geom_scatterpie(data=df5,
                    aes(x=celltype, y=sample, r=radius),
                    cols=colnames(df5)[3:(length(colnames(df5))-1)],
                    color=NA) +
    coord_equal()+
    geom_segment(aes(x=0.5,
                     xend=(length(unique(df5$celltype))+0.5),
                     y=-(length(unique(df5$sample))-0.5),
                     yend=-(length(unique(df5$sample))-0.5)),
                 linetype = 5) +
    
    geom_segment(aes(x=(length(unique(df5$celltype))-0.5),
                     xend=(length(unique(df5$celltype))-0.5),
                     y=-(length(unique(df5$sample))+0.5),
                     yend=-0.5),
                 linetype = 5) +
    geom_segment(aes(x=0.5, xend=(length(unique(df5$celltype))+0.5),
                     y=-0.5, yend=-0.5)) +
    geom_segment(aes(x=0.5, xend=(length(unique(df5$celltype))+0.5),
                     y=-(length(unique(df5$sample))+0.5), yend=-(length(unique(df5$sample))+0.5))) +
    geom_segment(aes(x=0.5, xend=0.5,
                     y=-0.5, yend=-(length(unique(df5$sample))+0.5))) +
    geom_segment(aes(x=(length(unique(df5$celltype))+0.5), xend=(length(unique(df5$celltype))+0.5),
                     y=-0.5, yend=-(length(unique(df5$sample))+0.5))) +
    geom_segment(aes(x=(length(unique(df5$celltype))+3),
                     xend=(length(unique(df5$celltype))+3),
                     y=-0.5,
                     yend=-(length(unique(df5$sample))+0.5)),
                 color ="white") +
    scale_fill_manual(values = color_ce) +
    scale_x_continuous("",expand = c(0,0),breaks = (1:length(unique(df5$celltype))),labels = celltype) +
    scale_y_continuous("",expand = c(0,0),breaks = (-length(unique(df5$sample))):-1,labels = rev(sample)) +
    theme_minimal_grid(color = "black") +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     size = 14),
          axis.text.y = element_text(size = 14)) +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.9, 0.9))
  
  
  
  ##12.2 添加圆圈legend
  ###12.2.1 准备画圈数据
  #diff <- c(lable[1],lable[1:(length(lable) - 1)])
  #lable <- sort(summary(df5$radius), decreasing = TRUE)
  #真实半径的统计值df5$radius列
  radius <- summary(df5$radius)[c("Min.", "Mean", "Max.")]
  lable_s <- labeller(radius)
  
  #绘制图例圆圈时的半径
  picture_radius <- signif(radius, digits = 2)
  
  #加标签时位置（半径）
  lable_radius <- lableReverse(picture_radius)
  
  ###12.2.3 添加圆圈legend
  if (length(unique(rowSums(df5[, 3:dim(df5)[2]]))) == 1) {
    p <- p1 +
      geom_circle(r = picture_radius[1], 
                  xc = length(unique(df5$celltype)) + 1.5, 
                  yc = -((length(unique(df5$sample))))) +
      geom_segment(aes(x =length(unique(df5$celltype)) + 1.5, 
                       xend = length(unique(df5$celltype)) + 1.5 + max(picture_radius) + 0.4*max(picture_radius),
                       y = -((length(unique(df5$sample)))) + lable_radius[1], 
                       yend = -((length(unique(df5$sample)))) + lable_radius[1])) +
      annotate("text",
               x = length(unique(df5$celltype)) + 1.5 + max(picture_radius) + 0.4*max(picture_radius),
               y = -((length(unique(df5$sample)))) + lable_radius[1],
               label = min(rowSums(df5[3:5])),
               hjust = 0)
  }
  
  if (length(unique(rowSums(df5[, 3:dim(df5)[2]]))) == 2) {
    p <- p1 +
      geom_circle(r = picture_radius[1], 
                  xc = length(unique(df5$celltype)) + 1.5, 
                  yc = -((length(unique(df5$sample))))) +
      geom_segment(aes(x =length(unique(df5$celltype)) + 1.5, 
                       xend = length(unique(df5$celltype)) + 1.5 + max(picture_radius) + 0.4*max(picture_radius),
                       y = -((length(unique(df5$sample)))) + lable_radius[1], 
                       yend = -((length(unique(df5$sample)))) + lable_radius[1])) +
      annotate("text",
               x = length(unique(df5$celltype)) + 1.5 + max(picture_radius) + 0.4*max(picture_radius),
               y = -((length(unique(df5$sample)))) + lable_radius[1],
               label = min(rowSums(df5[3:5])),
               hjust = 0) +
      
      geom_circle(r = picture_radius[3], 
                  xc = length(unique(df5$celltype)) + 1.5, 
                  yc = -((length(unique(df5$sample))))) +
      geom_segment(aes(x =length(unique(df5$celltype)) + 1.5, 
                       xend = length(unique(df5$celltype)) + 1.5 + max(picture_radius) + 0.4*max(picture_radius),
                       y = -((length(unique(df5$sample)))) + lable_radius[3], 
                       yend = -((length(unique(df5$sample)))) + lable_radius[3])) +
      annotate("text",
               x = length(unique(df5$celltype)) + 1.5 + max(picture_radius) + 0.4*max(picture_radius),
               y = -((length(unique(df5$sample)))) + lable_radius[3],
               label = max(rowSums(df5[3:5])),
               hjust = 0)
  }
  
  if (length(unique(rowSums(df5[, 3:dim(df5)[2]]))) > 2) {
    p <- p1 +
      geom_circle(r = picture_radius[1], 
                  xc = length(unique(df5$celltype)) + 1.5, 
                  yc = -((length(unique(df5$sample))))) +
      geom_segment(aes(x =length(unique(df5$celltype)) + 1.5, 
                       xend = length(unique(df5$celltype)) + 1.5 + max(picture_radius) + 0.4*max(picture_radius),
                       y = -((length(unique(df5$sample)))) + lable_radius[1], 
                       yend = -((length(unique(df5$sample)))) + lable_radius[1])) +
      annotate("text",
               x = length(unique(df5$celltype)) + 1.5 + max(picture_radius) + 0.4*max(picture_radius),
               y = -((length(unique(df5$sample)))) + lable_radius[1],
               label = min(rowSums(df5[3:5])),
               hjust = 0) +
      
      geom_circle(r = picture_radius[2], 
                  xc = length(unique(df5$celltype)) + 1.5, 
                  yc = -((length(unique(df5$sample))))) +
      geom_segment(aes(x =length(unique(df5$celltype)) + 1.5, 
                       xend = length(unique(df5$celltype)) + 1.5 + max(picture_radius) + 0.4*max(picture_radius),
                       y = -((length(unique(df5$sample)))) + lable_radius[2], 
                       yend = -((length(unique(df5$sample)))) + lable_radius[2])) +
      annotate("text",
               x = length(unique(df5$celltype)) + 1.5 + max(picture_radius) + 0.4*max(picture_radius),
               y = -((length(unique(df5$sample)))) + lable_radius[2],
               label = labeller(radius)[2],
               hjust = 0) +
      
      geom_circle(r = picture_radius[3], 
                  xc = length(unique(df5$celltype)) + 1.5, 
                  yc = -((length(unique(df5$sample))))) +
      geom_segment(aes(x =length(unique(df5$celltype)) + 1.5, 
                       xend = length(unique(df5$celltype)) + 1.5 + max(picture_radius) + 0.4*max(picture_radius),
                       y = -((length(unique(df5$sample)))) + lable_radius[3], 
                       yend = -((length(unique(df5$sample)))) + lable_radius[3])) +
      annotate("text",
               x = length(unique(df5$celltype)) + 1.5 + max(picture_radius) + 0.4*max(picture_radius),
               y = -((length(unique(df5$sample)))) + lable_radius[3],
               label = max(rowSums(df5[3:5])),
               hjust = 0)
  }
  
  if (sample_p == "T") {
    suffix = ".sc.CloneType.sample.Percent"
  }
  
  if (group_p == "T") {
    suffix = ".sc.CloneType.group.Percent"
  }
  
  pdf(str_c(outdir,"/06.sc.CloneTypePie/",prefix, suffix,".pdf"), width = 15,height = 10)
  print(p)
  dev.off()
  
  png(str_c(outdir,"/06.sc.CloneTypePie/",prefix, suffix,".png"), width = 1500,height = 1000)
  print(p)
  dev.off()
}
