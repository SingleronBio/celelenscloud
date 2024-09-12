#1 定义TCR/BCR所属链
TCR <- function(contig) {
  contig <- contig[contig$chain %in% c("TRA","TRB"),]
}

BCR <- function(contig) {
  contig <- contig[contig$chain %in% c("IGL","IGH","IGK"),]
}

#2 productive筛选
keepProductive <- function(contiglist) {
  
  newcontiglist <- list()
  
  for(x in 1:length(contiglist)) {
    contigi <- contiglist[[x]]
    newcontiglist[[x]] <- contigi[contigi$productive=="True"|contigi$productive== "true",]
  }
  
  return(newcontiglist)
}

#3 过滤掉只包含一条链的细胞
filteringSingle <- function(contig) {
  
  table <- subset(as.data.frame(table(contig$barcode)), Freq == 1 )
  barcodes <- as.character(unique(table$Var1))
  `%!in%` = Negate(`%in%`)
  contig <- subset(contig, barcode %!in% barcodes)
  
  return(contig)
}

#4 多链过滤——单个细胞中每个链可能有多个，此处保留两链各一个
##4.1 TCR过滤
filteringMulti <- function(contig) {
  
  table <- subset(as.data.frame(table(contig$barcode,contig$chain)), Freq > 1)
  barcodes <- as.character(unique(table$Var1))
  
  multichain <- NULL
  for (j in seq_along(barcodes)) {
    chain <- contig[contig$barcode == barcodes[j],] %>%
      group_by(barcode, chain) %>% top_n(n = 1, wt = reads) %>% group_by(barcode, chain) %>% top_n(n = 1, wt = umis)
    if(max(table(chain$chain)) < 2) {
      multichain <- rbind(multichain, chain) }
  }
  
  `%!in%` = Negate(`%in%`)
  contig <- subset(contig, barcode %!in% barcodes)
  contig <- rbind(contig, multichain)
  
  return(contig)
}

##4.2 BCR过滤
filteringMultiBCR <- function(contig) {
  tmp <- contig
  tmp$chain <- gsub("IGK","IGL",contig$chain)
  table <- subset(as.data.frame(table(tmp$barcode,tmp$chain)), Freq > 1)
  barcodes <- as.character(unique(table$Var1))
  
  multichain <- NULL
  for (j in seq_along(barcodes)) {
    chain <- contig[contig$barcode == barcodes[j],] %>%
      group_by(barcode, chain) %>% top_n(n = 1, wt = reads) %>% group_by(barcode, chain) %>% top_n(n = 1, wt = umis)
    if( sum(unique(chain$chain)  %in% c("IGL","IGK")) == 2) {
      if("IGH" %in% unique(chain$chain)) {
        IGH <- subset(chain, chain=="IGH")
        IGL <- subset(chain, chain %in% c("IGL","IGK")) %>% group_by(barcode) %>% top_n(n = 1, wt = reads) %>% top_n(n = 1, wt = umis)
        chain <- rbind(IGH,IGL)
      }else{
        chain <- chain %>% group_by(barcode) %>% top_n(n = 1, wt = reads) %>% top_n(n = 1, wt = umis)
      }
    }
    if(max(table(gsub("IGK","IGL",chain$chain))) < 2) {
      multichain <- rbind(multichain, chain) }
  }
  
  `%!in%` = Negate(`%in%`)
  contig <- subset(contig, barcode %!in% barcodes)
  contig <- rbind(contig, multichain)
  
  return(contig)
}

#5 main
filterVDJ <- function(contiglist, vdjType=vdjType,chairPair=chairPair) {

  newcontiglist <- NULL
  
  #1 Productive列过滤
  contiglist <- keepProductive(contiglist)
  
  #2 Productive列过滤后QC绘图
  plist <- dataQC(contiglist, vdjType)
  
  
  if(vdjType=="TCR") {
    #3 链过滤
    contiglist <- lapply(contiglist, TCR)
    
    #4 多链只保留一对
    for(x in 1:length(contiglist)) {
      newcontiglist[[x]] <- filteringMulti(contiglist[[x]])
    }
    
    #5 过滤掉只包含单链的细胞
    if(chairPair==TRUE) {
      for(x in 1:length(newcontiglist)) {
        newcontiglist[[x]] <- filteringSingle(newcontiglist[[x]])
      }
    }
  }
  else {
    contiglist <-lapply(contiglist,BCR)
    
    for(x in 1:length(contiglist)) {
      newcontiglist[[x]] <- filteringMultiBCR(contiglist[[x]])
    }
    
    if(chairPair==TRUE) {
      for(x in 1:length(newcontiglist)) {
        newcontiglist[[x]] <- filteringSingle(newcontiglist[[x]])
      }
    }
  }
  
  return(list(newcontiglist, plist))
}































