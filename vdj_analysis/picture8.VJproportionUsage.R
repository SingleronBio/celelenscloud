library(ComplexHeatmap)
library(tidyverse)

#1 画幅
#交给heatmap函数计算，不再多做计算

#2 按样本画图
scVDJUsage <- function(sc,datatype=dataType,outdir=outdir) {
  t <- sc
  t@meta.data$cluster <- t@active.ident
  if(datatype=="TCR") {
    
    #判断是否是测到D基因，是就添加D基因信息，否则不添加
    a <- unique(t@meta.data$TRBD)
    a[is.na(a)] <- "None"
    b <- c("None","None")
    
    if (all(a==b)) {
      print("this data had no D gene")
    } else {
      print("this data belong to overall length, and then add gene of D info")
      TRBD.prop  <- prop.table(x=table(t@meta.data$cluster,t@meta.data$TRBD),margin=1)
      
      #删除被测到但没有名称的基因（名称为"None"），下同
      TRBD.prop  <- TRBD.prop[, !(colnames(TRBD.prop) %in% c("")),drop=F]
    }
    
    TRAV.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$TRAV),margin=1)
    TRAV.prop  <- TRAV.prop[, !(colnames(TRAV.prop) %in% c("")),drop=F]
    
    TRAJ.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$TRAJ),margin=1)
    TRAJ.prop  <- TRAJ.prop[, !(colnames(TRAJ.prop) %in% c("")),drop=F]
    
    TRBV.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$TRBV),margin=1)
    TRBV.prop  <- TRBV.prop[, !(colnames(TRBV.prop) %in% c("")),drop=F]
    
    TRBJ.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$TRBJ),margin=1)
    TRBJ.prop  <- TRBJ.prop[, !(colnames(TRBJ.prop) %in% c("")),drop=F]
    
    #
    env <- ls()
    #
    if(any(str_detect(env, "^TRBD.prop"))) {
      TRB <- cbind(TRBV.prop,TRBD.prop,TRBJ.prop) %>% as.data.frame(stringsAsFactors=FALSE)
      TRB <- dplyr::select(TRB, colnames(TRB)[!(colnames(TRB) %in% c("None"))])
      col.split <- c(rep("TRBV",ncol(TRBV.prop)),rep("TRBD",ncol(TRBD.prop)),rep("TRBJ",ncol(TRBJ.prop)))
      print("this data belong to overall length, and then add gene of D info")
    } else {
      TRB <- cbind(TRBV.prop,TRBJ.prop) %>% as.data.frame(stringsAsFactors=FALSE)
      TRB <- dplyr::select(TRB, colnames(TRB)[!(colnames(TRB) %in% c("None"))])
      col.split <- c(rep("TRBV",ncol(TRBV.prop)),rep("TRBJ",ncol(TRBJ.prop)))
      print("this data had no D gene")
    }
    
    TRA <- cbind(TRAV.prop,TRAJ.prop) %>% as.data.frame(stringsAsFactors=FALSE)
    TRA <- dplyr::select(TRA, colnames(TRA)[!(colnames(TRA) %in% c("None"))])
    TRA.col.split <- c(rep("TRAV",ncol(TRAV.prop)),rep("TRAJ",ncol(TRAJ.prop)))
    
        
    #filter NA values and None col
    TRB <- TRB %>% drop_na() %>% as.matrix()
    TRA <- TRA %>% drop_na() %>% as.matrix()
    write.table(cbind(rownames(TRA),TRA), file=paste0(outdir, "/09.sc.VJproportionUsage/", prefix, ".", unique(t@meta.data$sample), ".sc.A.VJproportion.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(cbind(rownames(TRB),TRB), file=paste0(outdir, "/09.sc.VJproportionUsage/", prefix, ".", unique(t@meta.data$sample), ".sc.B.VJproportion.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)

	anno_TRB <- data.frame(Type = substr(colnames(TRB),1,4))
	rownames(anno_TRB) <- colnames(TRB)
	hp.TRB <- try(draw_ComplexHeatmap(TRB,scale_or_not = FALSE, cellwidth = 10, cellheight = 10, column_split = "colAnnotation",colAnnotation  = anno_TRB))

	anno_TRA <- data.frame(Type = substr(colnames(TRA),1,4))
    rownames(anno_TRA) <- colnames(TRA)    
	hp.TRA <- try(draw_ComplexHeatmap(TRA,scale_or_not = FALSE, cellwidth = 10, cellheight = 10, column_split = "colAnnotation",colAnnotation  = anno_TRA))

    #保存图片
	if(!"try-error" %in% class(hp.TRB)){
		save_heatmap(hp.TRB,pdf = TRUE, png = TRUE, prefix = paste0(prefix,".",unique(t@meta.data$sample),".sc.TRB.VJproportion"),outdir = paste0(outdir, "/09.sc.VJproportionUsage/"))
	}
	if(!"try-error" %in% class(hp.TRA)){
		save_heatmap(hp.TRA,pdf = TRUE, png = TRUE, prefix = paste0(prefix,".",unique(t@meta.data$sample),".sc.TRA.VJproportion"),outdir = paste0(outdir, "/09.sc.VJproportionUsage/"))
	}
    #return(c(hp.TRB,hp.TRA,colnames_len))
  }
  
  if(datatype=="BCR") {
    
    #判断是否是测到D基因，是就添加D基因信息，否则不添加
    a <- unique(t@meta.data$IGHD)
    a[is.na(a)] <- "None"
    b <- c("None","None")
    
    if (all(a==b)) {
      print("this data had no D gene")
    } else {
      print("this data belong to overall length, and then add gene of D info")
      IGHD.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$IGHD),margin=1)
      IGHD.prop  <- IGHD.prop[, !(colnames(IGHD.prop) %in% c("")),drop=F]
    }
    
    IGHV.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$IGHV),margin=1)
    IGHV.prop  <- IGHV.prop[, !(colnames(IGHV.prop) %in% c("")),drop=F]
    
    IGHJ.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$IGHJ),margin=1)
    IGHJ.prop  <- IGHJ.prop[, !(colnames(IGHJ.prop) %in% c("")),drop=F]
    
    IGLV.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$IGLV),margin=1)
    IGLV.prop  <- IGLV.prop[, !(colnames(IGLV.prop) %in% c("")),drop=F]
    
    IGLJ.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$IGLJ),margin=1)
    IGLJ.prop  <- IGLJ.prop[, !(colnames(IGLJ.prop) %in% c("")),drop=F]
    
    #
    env <- ls()
    
    if(any(str_detect(env, "^IGHD.prop"))) {
      IGH <- cbind(IGHV.prop,IGHD.prop,IGHJ.prop) %>% as.data.frame(stringsAsFactors=FALSE)
      IGH <- dplyr::select(IGH, colnames(IGH)[!(colnames(IGH) %in% c("None"))])
      col.split <- c(rep("IGHV",ncol(IGHV.prop)),rep("IGHD",ncol(IGHD.prop)),rep("IGHJ",ncol(IGHJ.prop)))
      print("this data belong to overall length, and then add gene of D info")
    }else {
      IGH <- cbind(IGHV.prop,IGHJ.prop) %>% as.data.frame(stringsAsFactors=FALSE)
      IGH <- dplyr::select(IGH, colnames(IGH)[!(colnames(IGH) %in% c("None"))])
      col.split <- c(rep("IGHV",ncol(IGHV.prop)),rep("IGHJ",ncol(IGHJ.prop)))
      print("this data had no D gene")
    }
    
    IGL <- cbind(IGLV.prop,IGLJ.prop) %>% as.data.frame(stringsAsFactors=FALSE)
    IGL <- dplyr::select(IGL, colnames(IGL)[!(colnames(IGL) %in% c("None"))])
    IGL.col.split <- c(rep("IGLV",ncol(IGLV.prop)),rep("IGLJ",ncol(IGLJ.prop)))
    
        
    #filter NA
    IGH <- IGH %>% drop_na() %>% as.matrix()
    IGL <- IGL %>% drop_na() %>% as.matrix()
    write.table(cbind(rownames(IGH),IGH),file=paste0(outdir,"/09.sc.VJproportionUsage/",prefix,".",unique(t@meta.data$sample),".sc.IGH.VJproportion.csv"),sep=",",quote = FALSE,row.names = FALSE,col.names = TRUE)
    write.table(cbind(rownames(IGL),IGL),file=paste0(outdir,"/09.sc.VJproportionUsage/",prefix,".",unique(t@meta.data$sample),".sc.IGL.VJproportion.csv"),sep=",",quote = FALSE,row.names = FALSE,col.names = TRUE)

    
    anno_IGH <- data.frame(Type = substr(colnames(IGH),1,4))
	rownames(anno_IGH) <- colnames(IGH)
	hp.IGH <- try(draw_ComplexHeatmap(IGH,scale_or_not = FALSE, cellwidth = 10, cellheight = 10, column_split = "colAnnotation",colAnnotation  = anno_IGH))

	anno_IGL <- data.frame(Type = substr(colnames(IGL),1,4))
    rownames(anno_IGL) <- colnames(IGL)    
	hp.IGL <- try(draw_ComplexHeatmap(IGL,scale_or_not = FALSE, cellwidth = 10, cellheight = 10, column_split = "colAnnotation",colAnnotation  = anno_IGL))

    #保存图片
	if(!"try-error" %in% class(hp.IGH)){
		save_heatmap(hp.IGH,pdf = TRUE, png = TRUE, prefix = paste0(prefix,".",unique(t@meta.data$sample),".sc.IGH.VJproportion"),outdir = paste0(outdir, "/09.sc.VJproportionUsage/"))
	}
	if(!"try-error" %in% class(hp.IGL)){
		save_heatmap(hp.IGL,pdf = TRUE, png = TRUE, prefix = paste0(prefix,".",unique(t@meta.data$sample),".sc.IGL.VJproportion"),outdir = paste0(outdir, "/09.sc.VJproportionUsage/"))
	}
    #return(c(hp.IGH,hp.IGL))
  }
}

#3 按组画图
scVDJUsageGroup <- function(sc,datatype=dataType,outdir=outdir) {
  t <- sc
  t@meta.data$cluster <- t@active.ident
  if(datatype=="TCR") {

    #判断是否是测到D基因，是就添加D基因信息，否则不添加
    a <- unique(t@meta.data$TRBD)
    a[is.na(a)] <- "None"
    b <- c("None","None")

    if (all(a==b)) {
      print("this data had no D gene")
    } else {
      print("this data belong to overall length, and then add gene of D info")
      TRBD.prop  <- prop.table(x=table(t@meta.data$cluster,t@meta.data$TRBD),margin=1)

      #删除被测到但没有名称的基因（名称为"None"），下同
      TRBD.prop  <- TRBD.prop[, !(colnames(TRBD.prop) %in% c("")),drop=F]
    }

    TRAV.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$TRAV),margin=1)
    TRAV.prop  <- TRAV.prop[, !(colnames(TRAV.prop) %in% c("")),drop=F]

    TRAJ.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$TRAJ),margin=1)
    TRAJ.prop  <- TRAJ.prop[, !(colnames(TRAJ.prop) %in% c("")),drop=F]

    TRBV.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$TRBV),margin=1)
    TRBV.prop  <- TRBV.prop[, !(colnames(TRBV.prop) %in% c("")),drop=F]

    TRBJ.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$TRBJ),margin=1)
    TRBJ.prop  <- TRBJ.prop[, !(colnames(TRBJ.prop) %in% c("")),drop=F]

    #
    env <- ls()
	#
    if(any(str_detect(env, "^TRBD.prop"))) {
      TRB <- cbind(TRBV.prop,TRBD.prop,TRBJ.prop) %>% as.data.frame(stringsAsFactors=FALSE)
      TRB <- dplyr::select(TRB, colnames(TRB)[!(colnames(TRB) %in% c("None"))])
      col.split <- c(rep("TRBV",ncol(TRBV.prop)),rep("TRBD",ncol(TRBD.prop)),rep("TRBJ",ncol(TRBJ.prop)))
      print("this data belong to overall length, and then add gene of D info")
    } else {
      TRB <- cbind(TRBV.prop,TRBJ.prop) %>% as.data.frame(stringsAsFactors=FALSE)
      TRB <- dplyr::select(TRB, colnames(TRB)[!(colnames(TRB) %in% c("None"))])
      col.split <- c(rep("TRBV",ncol(TRBV.prop)),rep("TRBJ",ncol(TRBJ.prop)))
      print("this data had no D gene")
    }

    TRA <- cbind(TRAV.prop,TRAJ.prop) %>% as.data.frame(stringsAsFactors=FALSE)
    TRA <- dplyr::select(TRA, colnames(TRA)[!(colnames(TRA) %in% c("None"))])
    TRA.col.split <- c(rep("TRAV",ncol(TRAV.prop)),rep("TRAJ",ncol(TRAJ.prop)))

        
    #filter NA values and None col
    TRB <- TRB %>% drop_na() %>% as.matrix()
    TRA <- TRA %>% drop_na() %>% as.matrix()
	write.table(cbind(rownames(TRA),TRA),file=paste0(outdir,"/09.sc.VJproportionUsage/",prefix,".",unique(t@meta.data$group),".sc.A.VJproportion.csv"),sep=",",quote = FALSE,row.names = FALSE,col.names = TRUE)
    write.table(cbind(rownames(TRB),TRB),file=paste0(outdir,"/09.sc.VJproportionUsage/",prefix,".",unique(t@meta.data$group),".sc.B.VJproportion.csv"),sep=",",quote = FALSE,row.names = FALSE,col.names = TRUE)

	anno_TRB <- data.frame(Type = substr(colnames(TRB),1,4))
	rownames(anno_TRB) <- colnames(TRB)	
	hp.TRB <- try(draw_ComplexHeatmap(TRB,scale_or_not = FALSE, cellwidth = 10, cellheight = 10, column_split = "colAnnotation",colAnnotation  = anno_TRB))

	anno_TRA <- data.frame(Type = substr(colnames(TRA),1,4))
    rownames(anno_TRA) <- colnames(TRA)    
	hp.TRA <- try(draw_ComplexHeatmap(TRA,scale_or_not = FALSE, cellwidth = 10, cellheight = 10, column_split = "colAnnotation",colAnnotation  = anno_TRA))
    #保存图片
	if(!"try-error" %in% class(hp.TRB)){
		save_heatmap(hp.TRB,pdf = TRUE, png = TRUE, prefix = paste0(prefix,".",unique(t@meta.data$group),".sc.TRB.VJproportion"),outdir = paste0(outdir, "/09.sc.VJproportionUsage/"))
	}
	if(!"try-error" %in% class(hp.TRA)){
		save_heatmap(hp.TRA,pdf = TRUE, png = TRUE, prefix = paste0(prefix,".",unique(t@meta.data$group),".sc.TRA.VJproportion"),outdir = paste0(outdir, "/09.sc.VJproportionUsage/"))
	}



    #return(c(hp.TRB,hp.TRA,colnames_len))
  }
	if(datatype=="BCR") {

    #判断是否是测到D基因，是就添加D基因信息，否则不添加
    a <- unique(t@meta.data$IGHD)
    a[is.na(a)] <- "None"
    b <- c("None","None")

    if (all(a==b)) {
      print("this data had no D gene")
    } else {
      print("this data belong to overall length, and then add gene of D info")
      IGHD.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$IGHD),margin=1)
      IGHD.prop  <- IGHD.prop[, !(colnames(IGHD.prop) %in% c("")),drop=F]
    }

    IGHV.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$IGHV),margin=1)
    IGHV.prop  <- IGHV.prop[, !(colnames(IGHV.prop) %in% c("")),drop=F]

    IGHJ.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$IGHJ),margin=1)
    IGHJ.prop  <- IGHJ.prop[, !(colnames(IGHJ.prop) %in% c("")),drop=F]

    IGLV.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$IGLV),margin=1)
    IGLV.prop  <- IGLV.prop[, !(colnames(IGLV.prop) %in% c("")),drop=F]

    IGLJ.prop <- prop.table(x=table(t@meta.data$cluster,t@meta.data$IGLJ),margin=1)
    IGLJ.prop  <- IGLJ.prop[, !(colnames(IGLJ.prop) %in% c("")),drop=F]

    #
    env <- ls()
	if(any(str_detect(env, "^IGHD.prop"))) {
      IGH <- cbind(IGHV.prop,IGHD.prop,IGHJ.prop) %>% as.data.frame(stringsAsFactors=FALSE)
      IGH <- dplyr::select(IGH, colnames(IGH)[!(colnames(IGH) %in% c("None"))])
      col.split <- c(rep("IGHV",ncol(IGHV.prop)),rep("IGHD",ncol(IGHD.prop)),rep("IGHJ",ncol(IGHJ.prop)))
      print("this data belong to overall length, and then add gene of D info")
    }else {
      IGH <- cbind(IGHV.prop,IGHJ.prop) %>% as.data.frame(stringsAsFactors=FALSE)
      IGH <- dplyr::select(IGH, colnames(IGH)[!(colnames(IGH) %in% c("None"))])
      col.split <- c(rep("IGHV",ncol(IGHV.prop)),rep("IGHJ",ncol(IGHJ.prop)))
      print("this data had no D gene")
    }

    IGL <- cbind(IGLV.prop,IGLJ.prop) %>% as.data.frame(stringsAsFactors=FALSE)
    IGL <- dplyr::select(IGL, colnames(IGL)[!(colnames(IGL) %in% c("None"))])
    IGL.col.split <- c(rep("IGLV",ncol(IGLV.prop)),rep("IGLJ",ncol(IGLJ.prop)))

    
    #filter NA
    IGH <- IGH %>% drop_na() %>% as.matrix()
	print(IGH)
    IGL <- IGL %>% drop_na() %>% as.matrix()
	print(IGL)
	write.table(cbind(rownames(IGH),IGH),file=paste0(outdir,"/09.sc.VJproportionUsage/",prefix,".",unique(t@meta.data$group),".sc.IGH.VJproportion.csv"),sep=",",quote = FALSE,row.names = FALSE,col.names = TRUE)
    write.table(cbind(rownames(IGL),IGL),file=paste0(outdir,"/09.sc.VJproportionUsage/",prefix,".",unique(t@meta.data$group),".sc.IGL.VJproportion.csv"),sep=",",quote = FALSE,row.names = FALSE,col.names = TRUE)

	anno_IGH <- data.frame(Type = substr(colnames(IGH),1,4))
	rownames(anno_IGH) <- colnames(IGH)
	hp.IGH <- try(draw_ComplexHeatmap(IGH,scale_or_not = FALSE, cellwidth = 10, cellheight = 10, column_split = "colAnnotation",colAnnotation  = anno_IGH))

	anno_IGL <- data.frame(Type = substr(colnames(IGL),1,4))
    rownames(anno_IGL) <- colnames(IGL)    
	hp.IGL <- try(draw_ComplexHeatmap(IGL,scale_or_not = FALSE, cellwidth = 10, cellheight = 10, column_split = "colAnnotation",colAnnotation  = anno_IGL))

    #保存图片
	if(! "try-error" %in% class(hp.IGH)){
		save_heatmap(hp.IGH,pdf = TRUE, png = TRUE, prefix = paste0(prefix,".",unique(t@meta.data$group),".sc.IGH.VJproportion"),outdir = paste0(outdir, "/09.sc.VJproportionUsage/"))
	}
	if(! "try-error" %in% class(hp.IGL)){
		save_heatmap(hp.IGL,pdf = TRUE, png = TRUE, prefix = paste0(prefix,".",unique(t@meta.data$group),".sc.IGL.VJproportion"),outdir = paste0(outdir, "/09.sc.VJproportionUsage/"))
	}
    #return(c(hp.IGH,hp.IGL))
  }
}				
