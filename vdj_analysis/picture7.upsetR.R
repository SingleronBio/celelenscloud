suppressMessages({
library (UpSetR)
library (reshape2)
library (dplyr)
library (tidyr)
library (ggplot2)
#library (Seurat)
})

clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
source('/Public/Script/shouhou/SCRIPT/Seurat_Monocle_modify/color_protocol.R')
clustcol <- c(color_protocol, clustcol)

#rds <- readRDS ("/Public/Script/shouhou/SCRIPT/SigVDJ/SigVDJ/result/test.rds")
upsetR <- function (rds,outdir,num=10,width=13,height = 7,cloneDefine=cloneDefine,color = NA,vdjtype) {
    meta <- rds@meta.data
    for (sp in unique (meta$sample)) {
		print(paste0("---------",sp))
        sub <- meta[meta$sample %in% sp,]
		if(cloneDefine=="CDR3"){
			sub.ctp <- dplyr::select (sub,c("CTaa","cluster"))
			sub.ctp <- sub.ctp[!is.na (sub.ctp$CTaa),]
	        sub.ctp$id <- paste0(sub.ctp$cluster,".",sub.ctp$CTaa)
	        freq <- dplyr::group_by (sub.ctp,id) %>% summarise (freq = n()) %>% as.data.frame ()
    	    sub.ctp.unique <- sub.ctp[order(sub.ctp$CTaa),] %>% unique() %>% inner_join (freq,by = "id") %>% dplyr::select (-id)
		}else if(cloneDefine == "VDJ"){
			sub.ctp <- dplyr::select (sub,c("CTgene","cluster"))
			sub.ctp <- sub.ctp[!is.na (sub.ctp$CTgene),]
    	    sub.ctp$id <- paste0(sub.ctp$cluster,".",sub.ctp$CTgene)
        	freq <- dplyr::group_by (sub.ctp,id) %>% summarise (freq = n()) %>% as.data.frame ()
	        sub.ctp.unique <- sub.ctp[order(sub.ctp$CTgene),] %>% unique() %>% inner_join (freq,by = "id") %>% dplyr::select (-id)
		}else if(cloneDefine == "VDJ_CDR3"){
			sub.ctp <- dplyr::select (sub,c("CTstrict","cluster"))
            sub.ctp <- sub.ctp[!is.na (sub.ctp$CTstrict),]
            sub.ctp$id <- paste0(sub.ctp$cluster,".",sub.ctp$CTstrict)
            freq <- dplyr::group_by (sub.ctp,id) %>% summarise (freq = n()) %>% as.data.frame ()
            sub.ctp.unique <- sub.ctp[order(sub.ctp$CTstrict),] %>% unique() %>% inner_join (freq,by = "id") %>% dplyr::select (-id)
		}else if(cloneDefine == "VDJ_CDR3aa"){
			sub.ctp <- dplyr::select (sub,c("CTstrictaa","cluster"))
            sub.ctp <- sub.ctp[!is.na (sub.ctp$CTstrict),]
            sub.ctp$id <- paste0(sub.ctp$cluster,".",sub.ctp$CTstrict)
            freq <- dplyr::group_by (sub.ctp,id) %>% summarise (freq = n()) %>% as.data.frame ()
			sub.ctp.unique <- sub.ctp[order(sub.ctp$CTstrictaa),] %>% unique() %>% inner_join (freq,by = "id") %>% dplyr::select (-id)
		}
        colnames (sub.ctp.unique) <- c("clonotype","celltype","freq")
        write.table (sub.ctp.unique, file = str_c(outdir,"/08.sc.UpsetR/", prefix, "_", sp, ".all.subcluster.upsetR.xls"), sep='\t', quote = F, col.names = T, row.names = F)
        ## plot data
        if (length(unique (sub$cluster)) > 1) {
            pd <- sub.ctp.unique  %>% dplyr::select (c("celltype","clonotype"))
            pd$value = 1
            tr <- spread (pd,celltype,value)
            tr[is.na(tr)] <- 0
            ##plot
			print(paste0("---------",sp))
			if(is.na(color[1])){
				col_use <- clustcol[1:(length(colnames(tr))-1)]
			}else{
				if(is.null(names(color))){
					col_use <- color
				}else{
					#col_use <- color[colnames(tr)[2:ncol(tr)]]
					pd_use <- data.frame(table(sub.ctp.unique$celltype))
					pd_use <- pd_use[order(pd_use$Freq,decreasing = T),]
					col_use <- color[as.character(pd_use$Var1)]
				}
			}
            p <- try((upset (tr, nsets = length (colnames (tr)), nintersects = num, mb.ratio = c(0.5, 0.5),
                    order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE),
                    mainbar.y.label = paste0("Number of shared","\n",vdjtype," clonotypes"),sets.x.label = "Cluster size",
                    main.bar.color='#502E91',
                    sets.bar.color = col_use,
                    matrix.color="black",line.size =0.75,point.size = 2.5,
                    text.scale = 1.8)),silent = T)
			if(! "try-error" %in% class(p)){
	            pdf (str_c(outdir,"/08.sc.UpsetR/", prefix, "_", sp, ".sc.upsetR.pdf"), width = width, height = height)
	            print (p)
	            dev.off()
	            
	            png (str_c(outdir,"/08.sc.UpsetR/", prefix, "_", sp, ".sc.upsetR.png"), units = 'in', width = width, height = height, res=300)
	            print (p)
	            dev.off()
	            
	            write.table (tr, file = str_c(outdir, "/08.sc.UpsetR/", prefix, "_", sp, "sc.upsetR.xls"), sep='\t', quote=F, row.names = F, col.names = T)
			}else{
				print("================upsetR=================")
				print(as.character(p))
				print("=======================================")
			}
        }
    }
}

