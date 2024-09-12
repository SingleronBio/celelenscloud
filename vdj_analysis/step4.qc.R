#1 dataQC——展示各个克隆型频率对应的细胞数
dataQC <- function(newcontiglist, vdjType) {
  
  plist <- list()
  replist <- list()
  
  #1.1 TCR
  if(vdjType=="TCR") {
    for(i in 1:length(newcontiglist)){
      
      #统计A链频率
      TRA <- newcontiglist[[i]][newcontiglist[[i]]$chain %in% c("TRA"),]
      TRA_freq <- as.data.frame(table(TRA$barcode))
      TRA_freq$Freq <- as.factor(TRA_freq$Freq)
      
      #A链画图
      plist[[1]] <- ggplot(TRA_freq, aes(x=Freq)) +
        geom_bar( colour="black", fill="grey") + 
        labs(x="frequency",y="count",title="TRA chain counts of per cell") +
        theme_classic()
      
      #统计B链频率
      TRB <- newcontiglist[[i]][newcontiglist[[i]]$chain %in% c("TRB"),]
      TRB_freq <- as.data.frame(table(TRB$barcode))
      TRB_freq$Freq <- as.factor(TRB_freq$Freq)
      
      #B链画图
      plist[[2]] <- ggplot(TRB_freq, aes(x=Freq)) +
        geom_bar(colour="black", fill="grey") +
        labs(x="frequency",y="count",title="TRB chain counts of per cell") + 
        theme_classic()
      
      #AB链拼图
      p <- ggarrange(plotlist=plist, ncol=2, nrow=1, align="v")
      
      replist[[i]] <- p
      
    }
    return(replist)
    
    #1.2 BCR
  } else {
    for(i in 1:length(newcontiglist)) {
      
      IGH <- newcontiglist[[i]][newcontiglist[[i]]$chain %in% c("IGH"),]
      IGH_freq <- as.data.frame(table(IGH$barcode))
      IGH_freq$Freq <- as.factor(IGH_freq$Freq)
      
      plist[[1]] <- ggplot(IGH_freq, aes(x=Freq)) +
        geom_bar(colour="black", fill="grey") + 
        theme_classic()+
        labs(x="frequency",y="count",title="IGH chain counts of per cell")
      
      IGL_K <- newcontiglist[[i]][newcontiglist[[i]]$chain %in% c("IGL","IGK"),]
      IGL_K_freq <- as.data.frame(table(IGL_K$barcode))
      IGL_K_freq$Freq <- as.factor(IGL_K_freq$Freq)
      
      plist[[2]] <- ggplot(IGL_K_freq, aes(x=Freq)) +
        geom_bar(colour="black", fill="grey") +
        theme_classic() +
        labs(x="frequency",y="count",title="IGL&IGK chain counts of per cell")
      
      p <- ggarrange(plotlist=plist,ncol=2, nrow=1, align="v")
      
      replist[[i]] <- p
      
    }
    return(replist)
  }
}

#2 展示过滤后每个样本中包含配对免疫受体或单链免疫受体细胞数
pairChainCells <- function(contiglist, samples){
  data <- data.frame()	
  for(i in seq_along(contiglist)) {
    pnum <- nrow(subset(contiglist[[i]],chain1!="None"&chain2!="None"))
    snum <- nrow(contiglist[[i]])-pnum
    data <- rbind(data,c(samples[i],pnum,snum))
  }
  colnames(data) <- c("sample","pnum","snum")
  return(data)
}




