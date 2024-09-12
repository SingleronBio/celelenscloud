#1 免疫组信息合并——按样本
##1.1 定义字符向量
HEAVYLINES <- c("chain1", "cdr3_aa1", "cdr3_nt1","VDJgeneH","IGHV","IGHD","IGHJ")
LIGHTLINES <- c("chain2", "cdr3_aa2", "cdr3_nt2","VDJgeneL","IGLV","IGLD","IGLJ")
l_lines <- c("IGLct", "cdr3", "cdr3_nt", "v_gene")
k_lines <- c("IGKct", "cdr3", "cdr3_nt", "v_gene")
h_lines <- c("IGHct", "cdr3", "cdr3_nt", "v_gene")
VDJGENES <- c("v_gene","d_gene","j_gene")

TRALINES <- c("chain1", "cdr3_aa1", "cdr3_nt1","VDJgeneA","TRAV","TRAD","TRAJ")
TRBLINES <- c("chain2", "cdr3_aa2", "cdr3_nt2","VDJgeneB","TRBV","TRBD","TRBJ")
CONTIGLINES <- c("chain", "cdr3", "cdr3_nt","VDJgene","v_gene","d_gene","j_gene")
CTLINES <- c("CTgene", "CTnt", "CTaa", "CTstrict","CTstrictaa", "cellType")

#' Author shaoyanyan
#' @param cells : assign T-AB cells or T-GD or B cells 
#' This function celltype specfic chain and cell
#' @return c(chain1, chain2, cellType): cellType-T-AB/T-GD/B
cellT <- function(cells) {
  if (cells == "T-AB") { 
    chain1 <- "TRA"
    chain2 <- "TRB" 
    cellType <- "T-AB" 
  } else if (cells == "T-GD") {
    chain1 <- "TRD"
    chain2 <- "TRG"
    cellType <- "T-GD" 
  } else if (cells == "B") {
    chain1 <- "IGH"
    chain2 <- "IGL"
    cellType <- "B" 
  }
  return(list(chain1, chain2, cellType))
}


#Sorting the V/D/J/C gene sequences for T and B cells
makeGenes <- function(cellType, data2, chain1, chain2) {
  if(cellType %in% c("T-AB", "T-GD")) {
    data2 <- data2 %>% 
      mutate(TCR1 = ifelse(chain == chain1, paste(with(data2, 
                                                       interaction(v_gene,  j_gene, c_gene))), NA)) %>%
      mutate(TCR2 = ifelse(chain == chain2, paste(with(data2, 
                                                       interaction(v_gene,  j_gene, d_gene, c_gene))), NA))
  }
  else {
    data2 <- data2 %>% 
      mutate(IGKct = ifelse(chain == "IGK", paste(with(data2, 
                                                       interaction(v_gene,  j_gene, c_gene))), NA)) %>%
      mutate(IGLct = ifelse(chain == "IGL", paste(with(data2, 
                                                       interaction(v_gene,  j_gene, c_gene))), NA)) %>%
      mutate(IGHct = ifelse(chain == "IGH", paste(with(data2, 
                                                       interaction(v_gene, j_gene, d_gene, c_gene))), NA))
  }
  return(data2)
  
}


#General combination of nucleotide, aa, and gene sequences for T/B cells
assignCT <- function(celltype, condf) {
  if (celltype %in% c("T-AB", "T-GD")) {
    condf$CTgene <- paste(condf$VDJgeneA, condf$VDJgeneB, sep="_")
    condf$CTnt <- paste(condf$cdr3_nt1, condf$cdr3_nt2, sep="_")
    condf$CTaa <- paste(condf$cdr3_aa1, condf$cdr3_aa2, sep="_")
    condf$CTstrict <- paste(condf$VDJgeneA, condf$cdr3_nt1,
                            condf$VDJgeneB, condf$cdr3_nt2, sep="_")
	condf$CTstrictaa <- paste(condf$VDJgeneA, condf$cdr3_aa1,
                            condf$VDJgeneB, condf$cdr3_aa2, sep="_")
  } else {
    condf$CTgene <- paste(condf$VDJgeneH, condf$VDJgeneL, sep="_")
    condf$CTnt <- paste(condf$cdr3_nt1, condf$cdr3_nt2, sep="_")
    condf$CTaa <- paste(condf$cdr3_aa1, condf$cdr3_aa2, sep="_") 
    condf$CTstrict <- paste(condf$VDJgeneH, condf$cdr3_nt1,
                            condf$VDJgeneL, condf$cdr3_nt2, sep="_")
	condf$CTstrictaa <- paste(condf$VDJgeneH, condf$cdr3_aa1,
                            condf$VDJgeneL, condf$cdr3_aa2, sep="_")

  }
  return(condf)
}

#' @param contig: the contig file of specfic sample
#' @param condf:  the contig file cells`s data frame 
#'				  rows of condf means the num of cells
#'		colnames of condf : 
#'		chain1 cdr3_aa1 cdr3_nt1 VDJgeneA TRAV TRAD TRAJ chain2 cdr3_aa2 cdr3_nt2 VDJgeneB TRBV TRBD TRBJ
#'		TRA                    V11_D2_J13 V11   D2  J13  TRB                     V13_D2_J1.3 V13 D2  J1.3
#' @return condf
parseTCR <- function(condf,contig) {
  uniqbarcode <- unique(condf$barcode)
  for (y in seq_along(uniqbarcode)) {
    barcode.i <- condf$barcode[y]
    location.i <- which(barcode.i == contig$barcode)
    chains <- unique(contig[location.i,"chain"])
    #for two chain cells: eg,one α&one β
    if (length(unique(chains)) == 2) {
      if(contig[location.i[1],c("chain")] %in% c("TRA","TRD")) {
        condf[y,TRALINES]<-contig[location.i[1],CONTIGLINES]
        condf[y,TRBLINES]<-contig[location.i[2],CONTIGLINES]
      }
      if(contig[location.i[1],c("chain")] %in% c("TRB","TRG")) {
        condf[y,TRBLINES]<-contig[location.i[1],CONTIGLINES]
        condf[y,TRALINES]<-contig[location.i[2],CONTIGLINES]
      }
    }
    else {
      class(contig$reads)
      contig1 <- contig[location.i,c("chain","reads")] %>% top_n(n = 1, wt = reads)	
      if(contig1$chain %in% c("TRA","TRD")) {
        condf[y,TRALINES]<-contig[location.i[1],CONTIGLINES]
        condf[y,TRBLINES]<- "None"
      }
      if(contig1$chain %in% c("TRB","TRG")) {
        condf[y,TRALINES]<- "None"
        condf[y,TRBLINES]<-contig[location.i[1],CONTIGLINES]
      }
    }
  }
  condf <- assignCT("T-AB",condf)	
  return(condf)
}

parseBCR <- function(condf,contig) {
  uniqbarcode <- unique(condf$barcode)
  for (y in seq_along(uniqbarcode)) {
    barcode.i <- condf$barcode[y]
    location.i <- which(barcode.i == contig$barcode)
    chains <- unique(contig[location.i,"chain"])
    #判断一下是否同时包含IGK和IGL
    if(sum(c("IGK","IGL") %in% chains) > 1) {
      contig <- contig[location.i,]  %>% group_by(barcode, chain) %>% top_n(n = 1, wt = reads)
      condf[y,LIGHTLINES] <- contig[which(contig$chain=="IGL"),][,CONTIGLINES]
      if(c("IGH") %in% chains) {
        condf[y,HEAVYLINES]<- contig[which(contig$chain=="IGH"),][,CONTIGLINES]
      }else{
        condf[y,HEAVYLINES]<- "None"		
      }
    }else if(length(chains) == 2) {
      if(contig[location.i[1],c("chain")] %in% c("IGH")) {
        condf[y,HEAVYLINES]<-contig[location.i[1],CONTIGLINES]
        condf[y,LIGHTLINES]<-contig[location.i[2],CONTIGLINES]
      } else {
        condf[y,HEAVYLINES]<-contig[location.i[2],CONTIGLINES]
        condf[y,LIGHTLINES]<-contig[location.i[1],CONTIGLINES]
      }
    }else if(length(chains) == 1) {
      if(contig[location.i[1],c("chain")] %in% c("IGH")) {
        condf[y,HEAVYLINES]<- contig[location.i[1],CONTIGLINES]
        condf[y,LIGHTLINES]<- "None"
      }
      else {
        condf[y,HEAVYLINES]<- "None"
        condf[y,LIGHTLINES]<- contig[location.i[1],CONTIGLINES]
      }
    }
  }
  condf <- assignCT("B",condf)
  return(condf)
}

combinedTCR <- function(contiglist, samples, groups) {
  condflist <- NULL
  for (i in seq_along(contiglist)) {
    data <- contiglist[[i]]
    unique_df <- unique(data$barcode)
    condf <- data.frame(matrix(NA, length(unique_df), 15))
    colnames(condf) <- c("barcode",TRALINES, TRBLINES)
    condf$barcode <- unique_df
    condf <- parseTCR(condf,as.data.frame(data))
    condf$ID <- samples[i]
    condf$groups <- groups[i]
    condflist[[i]] <- condf
  }
  names(condflist) <- samples
  return(condflist)
}

combinedBCR <- function(contiglist, samples, groups) {
  condflist <- NULL
  for (i in seq_along(contiglist)) {
    data <- contiglist[[i]]
    unique_df <- unique(data$barcode)
    condf <- data.frame(matrix(NA, length(unique_df), 15))
    colnames(condf) <- c("barcode",HEAVYLINES, LIGHTLINES)
    condf$barcode <- unique_df
    condf <- parseBCR(condf,data)
    condf$ID <- samples[i]
    condf$groups <- groups[i]
    condflist[[i]] <- condf
  }
  names(condflist) <- samples
  return(condflist)
}

#2 将筛选后得到的combind合并到meta.data
combinedJoinTometa.data <- function(combined, pro) {
  a <- data.frame()
  for (i in 1:length(combined)) {
    a <- rbind(a, combined[[i]])
  }
  
  pro@meta.data <- rownames_to_column(pro@meta.data, var = "barcode") %>%
    left_join(., a, by = "barcode") %>%
    column_to_rownames(var = "barcode")
  
  return(pro)
}
