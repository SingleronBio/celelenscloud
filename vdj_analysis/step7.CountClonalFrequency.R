#1 分组构建list
groupBy <- function(data, by) {
  
  #获取meta.data列名
  col_name <- colnames(data)
  
  #根绝by提取meta.data中by指定列的元素
  a <- as.character(unique(data[[match(get("by"), col_name)]]))
  
  #创建与a等长的list
  tem <- vector("list", length(a))
  names(tem) <- a
  
  #根据by指定列对meta.data进行分组
  for (i in 1:length(a)) {
    assign(a[i], data[data[[match(get("by"), col_name)]] == a[i] ,])
    tem[[i]] <- get(a[i])
  }
  
  return(tem)
}

#2 克隆型统计
CountClonalFrequency <- function(t, by, cloneDefine) {
  pro <- t
  col_name <- colnames(pro@meta.data)
  
  #1 根据指定键分组
  tem <- groupBy(data = pro@meta.data, by = by)
  
  #2 cloneDefine参数确定根据哪一列统计克隆型频率
  if (cloneDefine == "VDJ") {
    define = "CTgene"
  } else if(cloneDefine == "CDR3"){
    define = "CTaa"
  } else if(cloneDefine == "VDJ_CDR3"){
	define = "CTstrict"
  } else if(cloneDefine == "VDJ_CDR3aa"){
 	define = "CTstrictaa"
  }
  
  #3 根据by参数确定克隆型频率的列名
  if (by == "sample") {
    freColname = "Frequency"
  } else {
    freColname = "group_freq"
  }
  
  #4 统计克隆型频率
  df <- data.frame()
  for (i in 1:length(tem)) {
    #使用table()函数对克隆型频率进行统计，并修改列名
    b <- as.data.frame(table(tem[[i]][[match(get("define"), col_name)]]), stringsAsFactors = F) %>%
      rename_with(~ c(get("define"),get("freColname")), everything())
    
    #将频率统计结果加入数据框
    #b$value =as.numeric(b$value)
    tem[[i]] <- rownames_to_column(tem[[i]], var = "cell_id") %>% 
      left_join(b, by = get("define"))
    
    #合并数据框
    df <- rbind(df, tem[[i]])
  }
  
  #5 将统计的克隆型频率合并到meta.data
  pro@meta.data <- rownames_to_column(pro@meta.data, var = "cell_id") %>%
    left_join(., dplyr::select(df, c("cell_id", get("freColname"))), by = "cell_id") %>%
    column_to_rownames(var = "cell_id")
  
  #6 结果返回
  combindSc <- list()
  combindSc[[1]] <- tem
  combindSc[[2]] <- pro
  
  return(combindSc)
}



