#1 按cluster及sample构建子集
subsetSc <- function() {
  
  #判断sample输入参数是否有误
  if (all(samples %in% as.character(unique(pro@meta.data$sample)))) {
    
  } else {
    print("error: some samples which provided not in the meta.data, plase cheak the input sample")
    quit()
  }
  
  #根据cluster构建子集
  if (identical(cluster ,"all")) {
    print("all cluster used to analysis")
  } else {
    pro <- subset(pro, idents = cluster)
    print(str_c("the ident after Process is ", str_c(sort(as.character(levels(pro))), collapse = "/ "), sep = ": "))
  }
  
  #根据sample构建子集
  if (identical(sort(samples) , sort(as.character(unique(pro@meta.data$sample))))) {
    print("all sample used to analysis")
  } else {
    Idents(pro) <- "sample"
    pro <- subset(pro, idents = samples)
    Idents(pro) <- "cluster"
    print(str_c("the sample after Process is ", str_c(sort(as.character(unique(pro@meta.data$sample))), collapse = "/ "), sep = ": "))
  }
  return(pro)
}

#2 向meta.data中加入group列
addColumngroup <- function() {
  pro@meta.data$group <- pro@meta.data$sample
  
  for (i in 1:length(samples)) {
    pro@meta.data$group <- str_replace_all(pro@meta.data$group, str_c("^", samples[i],"$"), groups[i])
  }
  
  return(pro)
}













