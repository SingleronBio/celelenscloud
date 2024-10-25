##3.1 按cluster及sample构建子集
subsetSc <- function(pro, subcluster, subsample) {
  
  #根据cluster构建子集
  if (subcluster != "all") {
    if (all(subcluster %in% levels(pro))){
      pro <- subset(pro, idents = subcluster)
      print(str_c("the ident after Process is ", str_c(sort(as.character(levels(pro))), collapse = "/ "), sep = ": "))
    } else {
      print("error: some cluster which provided not in the pro@meta.data, plase cheak the input cluster")
      quit()
    }
  } else {
    print("all cluster used to analysis")
  }
  
  #根据sample构建子集
  if (subsample != "all") {
    if (all(subsample %in% unique(pro@meta.data$sample))){
      Idents(pro) <- "sample"
      pro <- subset(pro, idents = subsample)
      Idents(pro) <- "cluster"
      print(str_c("the sample after Process is ", str_c(sort(as.character(unique(pro@meta.data$sample))), collapse = "/ "), sep = ": "))
    } else {
      print("error: some samples which provided not in the pro@meta.data, plase cheak the input samples")
      quit()
    }
  } else {
    print("all sample used to analysis")
  }
  return(pro)
}

##3.2 对cluster进行提取细胞
subsetClustercount <- function(pro) {
  if (subset == "F") {
    
    #不按细胞类型subset
    pro <- pro
  } else if (subset == "default") {
    
    #使用默认标准按细胞类型subset
    if (dim(pro@meta.data)[1] > 10000){
      pro <- subset(pro, downsample = round(10000 / length(levels(pro))))
    } else {
      pro <- pro
    }
    
  } else {
    
    #自定义数值，按细胞类型subset
    subset <- as.numeric(subset)
    if (!(is.numeric(subset))) {
      print("error: parameter of subset is wrong, please cheak 'subset'")
      quit()
    } else {
      pro <- subset(pro, downsample = subset)
    }
  }
  
  return(pro)
}

##3.3 修改cluster因子水平（多个细分亚群合并后可能造成cluster因子水平混乱）
clusterLevels <- function(pro, clusterlevels) {
  if (clusterlevels == "F") {
    pro <- pro
  } else if (clusterlevels == "default") {
    pro@meta.data$cluster <- factor(pro@meta.data$cluster, levels = sort(unique(as.character(pro@meta.data$cluster))))
    Idents(pro) <- "cluster"
  } else {
    pro@meta.data$cluster <- factor(pro@meta.data$cluster, levels = clusterlevels)
    Idents(pro) <- "cluster"
  }
  return(pro)
}

##3.4 对meta.data增加列(group/platform)
addColumn <- function(pro, addgorup, addplatform) {
  #添加group列
  if (addgorup == "T") {
    if (length(sample) == length(group)) {
      pro@meta.data$group <- pro@meta.data$sample
      
      for (i in 1:length(sample)) {
        pro@meta.data$group <- str_replace_all(pro@meta.data$group, str_c("^", sample[i],"$"), group[i])
      }
    } else {
      print("error: the length of sample is not match with the length of group, please cheak the input parameter 'sample', 'group'")
      quit()
    }
  }
  
  #添加platform列
  if (addplatform == "T") {
    if (length(sample) == length(platform)) {
      pro@meta.data$platform <- pro@meta.data$sample
      
      for (i in 1:length(sample)) {
        pro@meta.data$platform <- str_replace_all(pro@meta.data$platform, str_c("^",sample[i],"$"), platform[i])
      }
    } else {
      print("error: the length of sample is not match with the length of platform, please cheak the input parameter 'sample', 'platform'")
      quit()
    }
  }
  
  return(pro)
}

##3.5 指定rootcluster
pointRoot <- function(my_cds, rootcluster) {
  if (length(unique(pData(my_cds)$State)) > 1) {
    rootState <- as.numeric(as.data.frame(table(pData(my_cds)$State, pData(my_cds)$cluster)) %>%
                              dplyr::filter(Var2 == rootcluster) %>%
                              dplyr::filter(Freq == max(Freq)) %>%
                              pull(Var1) %>%
                              as.character())
    if (length(unique(rootState)) == 1) {
      return(rootState)
    } else {
      print(str_c("error: when set ", rootcluster, "as root, the unmber of root state was not equal unique"))
      quit()
    }
    
  } else {
    return(1)
  }
}

##3.6 设置配色版
colorSet <- function(newcol) {
  clustcol <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4",
                "#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00",
                "#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#E6E6FA","#FFDAB9")
  
  #判断是否使用新色板
  if (newcol == 'T') {
    color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
    clustcol <- c(color_protocol, clustcol)
  } else {
    clustcol<-c("red","blue","orange","green","Purple","black","yellow","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF",
                "#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#6A5ACD","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355",
                "#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#E6E6FA","#FFDAB9")
  }
  
  return(clustcol)
}
