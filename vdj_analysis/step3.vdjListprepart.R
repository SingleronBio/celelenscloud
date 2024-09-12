#1 读取vdjlist数据，并将其存入list
mergeSample <- function(vdjlist,deal_barcode) {
  contiglist <- list()
  
  for(i in 1:length(vdjlist)) {
	if(deal_barcode == "no"){
	    contiglist[[i]] <- read.csv(vdjlist[i], stringsAsFactors = F, sep = ",")
	}else if(deal_barcode == "yes"){
		tmp <- read.csv(vdjlist[i], stringsAsFactors = F, sep = ",")
		tmp$barcode <- gsub("[-][0-9]$","",tmp$barcode)
		contiglist[[i]] <- tmp
		print(head(tmp))
	}
  }

  return(contiglist)
}

#2 在免疫组标准分析文件baecode前加上文库前缀（需要与转录组barcode前缀一致，一般为文库名）
addBarcodeprefix <- function(contiglist, spname) {
  outcontiglist <- list()
  
  for(i in 1:length(contiglist)) {
    contiglist[[i]]$barcode <- str_c(spname[i], "_", contiglist[[i]]$barcode)
    outcontiglist[[i]] <- contiglist[[i]]
  }
  
  return(outcontiglist)
}

#3 添加VDJgene列，以“-”连接v_gene、d_gene、j_gene
addVDJgenecol <- function(contiglist) {
  
  newcontiglist <- list()
  
  for(x in 1:length(contiglist)) {
    contiglist[[x]]$VDJgene <- paste(contiglist[[x]]$v_gene, contiglist[[x]]$d_gene, contiglist[[x]]$j_gene, sep = "-")
  }
  return(contiglist)
}






