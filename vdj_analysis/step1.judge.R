#1 获取输入rds的barcode前缀名称，确定barcode的标准：字符长度大于等于16；字符只包含A/T/C/G或A/T/C/G/-/1
cheackSample <- function(data){
  
  #获取单细胞数据barcode，每个样本提取一个barcode
  a <-  unique(data$sample)
  sc_rownames <- c() 
  for (i in 1:length(a)) {
    assign(paste("df",a[i],sep="_"), filter(data, sample == a[i]))
    df <- get(paste("df",a[i],sep="_"))
    sc_rownames <- c(sc_rownames, rownames(df)[1])
  }
  
  #指定barcode碱基模式，用于识别barcode，从而从barcode中分隔出barcode前缀
  pattern <- c("A","C","G","T","-","1","_")
  
  
  sc_spname <- c()
  
  #使用for循环遍历sc_rownames，从中提取barcode前缀
  for (s_barcode in sc_rownames) {
    ##这里由于barocde格式调整，需要修改对b的处理方法
	if(grepl("[ATCG]{8}_[ATCG]{8}_[ATCG]{8}.*",s_barcode)){
        tmp1 <- gsub("_[ATCG]{8}_[ATCG]{8}_[ATCG]{8}.*","",s_barcode)
        tmp2 <- gsub(paste0(tmp1,"_"),"",s_bracode)
        b <- c(tmp1,tmp2)
    }else{
        tmp1 <- gsub("_[ATCG]{8,}.*$","",s_barcode)
        tmp2 <- gsub(paste0(tmp1,"_"),"",s_barcode)
        b <- c(tmp1,tmp2)
    }
    #字符长度大于等于16
    c <- b[str_length(b)>=16]  
    
    #字符中包含“1”或“—”的数目小于等于1
    c2 <- c()
    for (i in c) {
      
      #barcode中要么没有“1”和“-”，有的话它俩数目相等，且<=1
      i1 <- str_count(i, "1")
      i2 <- str_count(i, "-")
      if ((i1 == i2) & (i1 <= 1)) {
        c2 <- c(c2, i)
      }
    }
    
    e <- str_split(c2,"")
    
    #判断barcode碱基位置
    for (i in 1:length(e)) {
      
      #只含有pattern
      if(all(unique(e[[i]]) %in% pattern)){
        site_tmp <- i
        break
      }
    }
    
    #根据barcode碱基临时位置，提取barcode碱基,并确定barcode碱基
    site <- match(c2[site_tmp],b)
    
    #根据barcode碱基碱基位置对barcode前缀进行连接（之前可能被分隔）
    barcode_spname <- str_c(b[1:site-1],collapse = "_")
    
    sc_spname <- c(sc_spname, barcode_spname)
  }
  
  if (all(spname %in% unique(sc_spname))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#2 判断sample传入参数是否被meta.data$sample包含，防止手误输错样本（此处在小本本上写上方晟的名字，并敲一下他的小脑壳）
cheackSampleCount <- function(data) {
  if (length(vdjlist) != length(spname)){
    print("the length of spname/ vdjlist/ is unequal, plase cheak the input spname/ vdjlist/")
    quit()
  }
  
  if (length(vdjlist) != length(samples)){
    print("the length of sample/ vdjlist/ is unequal, plase cheak the input sample/ vdjlist/")
    quit()
  }
  
  if (length(vdjlist) != length(groups)){
    print("the length of groups/ vdjlist/ is unequal, plase cheak the input groups/ vdjlist/")
    quit()
  }
}


