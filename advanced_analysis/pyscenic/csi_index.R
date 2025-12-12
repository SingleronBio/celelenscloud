#! /usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(argparser)
library(pheatmap)

argv <- arg.parser("")
argv <- add_argument(argv,"--AUC",help = "use for person correltation coeffcient")
argv <- add_argument(argv,"--regulon_list",help = "regulons use in cluster")
argv <- add_argument(argv,"--prefix",help = "prefix")
argv <- add_argument(argv,"--outdir",help = "outdir")
argv <- add_argument(argv,"--k",help = "cluster number",default = 6)
argv <- parse.args(argv)

#####传参
if(grepl(".csv",argv$AUC,)){
        auc <- read.csv(argv$AUC,header = T,row.names = 1)
}else{
        auc <- read.table(argv$AUC,header= T,sep = "\t",row.names = 1)
}
print("AUC matrix has read")
regulon_list <- read.table(argv$regulon_list,header = F,sep = "\t")
regulon_use <- as.character(regulon_list[,1])
if( length(regulon_use[!regulon_use %in% rownames(auc)]) != 0 ){
        print(paste0("These regulons are not in auc matrix: ",paste(regulon_use[!regulon_use %in% rownames(auc)],collapse = ",")))
        regulon_use <- regulon_use[regulon_use %in% rownames(auc)]
}
print(regulon_use)
prefix <- argv$prefix
outdir <- argv$outdir
if(!file.exists(outdir)){
        dir.create(outdir)
}
#计算PCC矩阵
PCC_mat <- matrix(nrow = length(regulon_use),ncol = length(regulon_use))
colnames(PCC_mat) <- regulon_use
rownames(PCC_mat) <- regulon_use
#改进计算方法，省略冗余计算
for(i in 1:ncol(PCC_mat)){
#按列选
        com_1 <- colnames(PCC_mat)[i]
        for(n in 1:i){
#固定行后，按列选
                com_2 <- rownames(PCC_mat)[n]
                cor_num <- cor(as.numeric(as.character(auc[com_1,])),as.numeric(as.character(auc[com_2,])),method = "pearson")
                PCC_mat[com_1,com_2] <- cor_num
                PCC_mat[com_2,com_1] <- cor_num
                }
                width <- options()$width
                cat("[",paste(rep("-",(i/ncol(PCC_mat)*width)/4),collapse=""),paste(rep(" ",(width - i/ncol(PCC_mat)*width)/4),collapse = ""),"]",round(i/ncol(PCC_mat)*100),"%")
                Sys.sleep(0.05)
                if(i == ncol(PCC_mat)){
                        cat("\nDone!\n")
                }else{
                        cat("\r")
                }
}
print(PCC_mat)
CSI_mat <- PCC_mat
write.table(PCC_mat,paste0(outdir,"/",prefix,"_PCC_matrix.xls"),sep = "\t",quote = F,col.names = NA)
###计算CSI,改进计算方法，省去冗余计算
for(i in 1:ncol(PCC_mat)){
        #选定比较对象1,按列
        com_1 <- colnames(PCC_mat)[i]
        for(n in 1:i){
                #选定比较对象2，按列
                com_2 <- rownames(PCC_mat)[n]
                ref <- as.numeric(as.character(PCC_mat[com_1,com_2])) - 0.05
                #node_num 用于记录符合条件的node数(即同时小于ref)
                node_num <- 0
                for(m in regulon_use){
                        print(m)
                        if(m != com_1 & m != com_2){
                                tmp1 <- as.numeric(as.character(PCC_mat[com_1,m]))
                                tmp2 <- as.numeric(as.character(PCC_mat[com_2,m]))
                                if(tmp1 < ref & tmp2 < ref){
                                        node_num <- node_num + 1
                                }
                        }
                }
                csi <- node_num/length(regulon_use)
                CSI_mat[com_1,com_2] <- csi
                CSI_mat[com_2,com_1] <- csi
        }
}
print(CSI_mat)
write.table(CSI_mat,paste0(outdir,"/",prefix,"_CSI_mat.tsv"),sep = "\t",quote= F,col.names =NA)

p <- pheatmap(CSI_mat,cutree_rows = argv$k,cutree_cols = argv$k)
pdf(paste0(outdir,"/",prefix,"_CSI_heatmap.pdf"))
print(p)
dev.off()
png(paste0(outdir,"/",prefix,"_CSI_heatmap.png"))
print(p)
dev.off()
ggsave(file=paste0(outdir,"/",prefix,"_CSI_heatmap.svg"), plot=p)

reg_cluster <- cutree(p$tree_row,k=argv$k)
reg_cluster <- data.frame(reg_cluster)
colnames(reg_cluster) <- "Cluster"
reg_cluster$Regulon <- rownames(reg_cluster)
reg_cluster <- reg_cluster[order(reg_cluster$Cluster),]
write.table(reg_cluster,paste0(outdir,"/",prefix,"_cluster.xls"),sep = "\t",row.names = F,quote = F)
