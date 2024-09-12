#! /usr/bin/env Rscript
library(dplyr)
library(immunarch)
library(rlang)
##CTgene为VDJ基因型，CTaa为CDR3序列，根据cloneDefine确定用什么做统计绘图，用于绘图的一列会被命名为CDR3.aa
#source('/Public/Script/shouhou/SCRIPT/Seurat_Monocle_modify/color_protocol.R')
col_use <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B",
              "#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101",
              "#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B",
              "#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B",
              "#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060",
              "#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC",
              "#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8")

theCall <- function(x) {
    if (x == "gene") {
        x <- "CTgene"
    } else if(x == "nt") {
        x <- "CTnt"
    } else if (x == "aa") {
        x <- "CTaa"
    } else if (x == "gene+nt") {
        x <- "CTstrict"
    }
    return(x)
}

transform_data_for_immunarch <- function(sample_info_list,type,cloneDefine){
#	ct.define <- theCall("gene")
    sp <- names(sample_info_list)
    for(i in 1:length(sp)){
        if(i == 1){
            meta_data <- sample_info_list[[i]]
        }else{
            meta_data <- rbind(meta_data,sample_info_list[[i]])
        }
    }
    if(type == "TCR"){
#       meta_data<-meta_data[!is.na(meta_data$barcode),c("CTgene","CTnt","CTaa","CTstrict","TRAV","TRAD","TRAJ","TRBV","TRBD","TRBJ","sample","group")]
        meta_data<-meta_data[!is.na(meta_data$CTgene),c("CTgene","CTaa","CTstrict","CTstrictaa","sample","group")]
    }else{
#       meta_data<-meta_data[!is.na(meta_data$barcode),c("CTgene","CTnt","CTaa","CTstrict","IGHV","IGHD","IGHJ","IGLV","IGLD","IGLJ","sample","group")]
        meta_data<-meta_data[!is.na(meta_data$CTgene),c("CTgene","CTaa","CTstrict","CTstrictaa","sample","group")]
    }
#   colnames(meta_data)[which(colnames(meta_data)=="CTaa")]<-"CDR3"
#   colnames(meta_data)[which(colnames(meta_data)=="CTnt")]<-"VDJ"
    output<-list()
    output$data<-list()
    sp <- unique(as.character(meta_data$sample))
    for(i in sp){
        tmp<-meta_data[meta_data$sample==i,]
        if(cloneDefine == "VDJ"){
            ctstrict <- tmp$CTgene
        }else if(cloneDefine == "CDR3"){
            ctstrict <- tmp$CTaa
        }else if(cloneDefine == "VDJ_CDR3"){
			ctstrict <- tmp$CTstrict
		}else if(cloneDefine == "VDJ_CDR3aa"){
			 ctstrict <- tmp$CTstrictaa
		}
        if(length(ctstrict) != 0 ){
            freq <- c()
            for(n in 1:length(ctstrict)){
                freq_tmp <- length(ctstrict[ctstrict == ctstrict[n]])
                freq <- c(freq,freq_tmp)
                }
                tmp$Clones <- freq
                tmp$Proportion <- tmp$Clones/nrow(tmp)
                tmp <- tmp %>% dplyr::select("Clones","Proportion",everything())
                tmp <- data.frame(tmp)
                if(cloneDefine == "CDR3"){
                    colnames(tmp)[which(colnames(tmp) == "CTaa")] <- "CDR3.aa"
					tmp <- tmp[!duplicated(tmp[,"CDR3.aa"]),]
                    print("Clone type defined by CDR3")
                }else if(cloneDefine == "VDJ"){
                    colnames(tmp)[which(colnames(tmp) == "CTgene")] <- "CDR3.aa"
					tmp <- tmp[!duplicated(tmp[,"CDR3.aa"]),]
                    print("Clone type defined by VDJ")
                }else if(cloneDefine == "VDJ_CDR3"){
					colnames(tmp)[which(colnames(tmp) == "CTstrict")] <- "CDR3.aa"
					tmp <- tmp[!duplicated(tmp[,"CDR3.aa"]),]
					print("Clone type defined by VDJ_CDR3")
				}else if(cloneDefine == "VDJ_CDR3aa"){
					colnames(tmp)[which(colnames(tmp) == "CTstrictaa")] <- "CDR3.aa"
					tmp <- tmp[!duplicated(tmp[,"CDR3.aa"]),]
					print("Clone type defined by VDJ_CDR3aa")
				}
                tmp<-as_tibble(tmp)
                output$data[[i]]<-tmp
        }else{
                print(paste0("sample ",i," has no matched VDJ info."))
        }
    }
    sp_info <- data.frame(table(meta_data$sample,meta_data$group))
    sp_info <- sp_info[sp_info$Freq != 0,]
    sp_info <- sp_info[,c(1,2)]
    colnames(sp_info) <- c("Sample","Group")
    sp_info <- as_tibble(sp_info)
    output$meta <- sp_info
    return(output)
}




#####绘制所有immunarch
#####V2脚本中的多样性分析依据的可以是CTgene，CTstrict，CTaa
plot_diversity <- function(immunarch_data,method=NULL,outdir,plot_by = "Group",test = TRUE,color = NA){
    ######绘制的为单细胞转录组与免疫组库结合的免疫组库多样性结果
	if(!is.na(color[1])){
		col_use <- color
	}
    type <- method %||% c( "chao1","hill","inv.simp","d50")
    print(type)
	dir.create(paste0(outdir,"/02.sc.Diversity/",plot_by))
    for(i in type){
        p <- try(repDiversity(immunarch_data$data,.method = i) %>% vis(.by = plot_by,.meta = immunarch_data$meta,.test = test)+labs(x="")+scale_fill_manual(values = col_use),silent = T)
        if(!"try-error" %in% class(p)){
            pdf(paste0(outdir,"/02.sc.Diversity/",plot_by,"/",prefix,"_",i,"_diversity.pdf"))
            print(p)
            dev.off()
            png(paste0(outdir,"/02.sc.Diversity/",plot_by,"/",prefix,"_",i,"_diversity.png"))
            print(p)
            dev.off()
            ori_data <- data.frame(repDiversity(immunarch_data$data,.method = i))
            meta <- data.frame(immunarch_data$meta)
            gp <- as.character(meta$Group)
            names(gp) <- as.character(meta$Sample)
            if("Sample" %in% colnames(ori_data)){
                ori_data <- data.frame(group = gp[as.character(ori_data$Sample)],ori_data)
            }else{
                ori_data <- data.frame(group = gp[rownames(ori_data)],ori_data)
            }
            write.table(ori_data,paste0(outdir,"/02.sc.Diversity/",plot_by,"/",i,"_diversity.xls"),sep = "\t",quote= F,col.names = NA)
        }else{
#           print(as.character(p))
            print(paste0("Warning: Diversity plot of ",i," can't be ploted because some reasons."))
        }
    }
}

##############绘制hclust图，根据cloneDefine确定用哪部分(CTgene或者CTaa或者CTstrict)

plot_hclust <- function(immunarch_data,prefix,outdir){
	sp <- names(immunarch_data$data)
	for(i in sp){
		tmp <- data.frame(table(immunarch_data$data[[i]][,"CDR3.aa"]))
		print(head(tmp))
		rownames(tmp) <- tmp[,"CDR3.aa"]
		tmp <- tmp[,-which(colnames(tmp) == "CDR3.aa"),drop = F]
		tmp$Freq <- as.numeric(as.character(tmp$Freq))
		colnames(tmp)[1] <- i
		tmp[,i] <- tmp[,i]/sum(tmp[,i])
		if(which(sp == i) == 1){
			output <- tmp
		}else{
			tmp$use <- rownames(tmp)
			output$use <- rownames(output)
			output <- merge(output,tmp,by = "use",all = T)
			rownames(output) <- output$use
			output <- output[,-1,drop = F]
		}
	}
	output$Names <- rownames(output)
	output <- output %>% dplyr::select("Names",everything())
	output <- as_tibble(output)
	p <- vis(geneUsageAnalysis(output,"cosine+hclust"))
	k <- ggplot_build(p[[2]])$data[[3]]$xintercept
	p <- vis(geneUsageAnalysis(output,"cosine+hclust",.k=1))
	p <- p[[1]]+coord_cartesian(clip="off")+theme(plot.margin=unit(c(0,0,1,1),"cm"))
	pdf(paste0(outdir,"/03.sc.Clustertree/",prefix,"_hclust.pdf"),width = 13)
	print(p)
	dev.off()
	png(paste0(outdir,"/03.sc.Clustertree/",prefix,"_hclust.png"),width = 13/7*480)
	print(p)
	dev.off()
}
