library(pheatmap)
library (dplyr)
library (grid)
library (reshape2)
library (tidyverse)
library(plyr) 

source("/opt/lims_module_script/Heterogeneity_inferCNV/heatmap_function.R")
args<-commandArgs(T)
input=args[1]
dio=args[2]
cnvfilter=as.numeric(args[3])
armfilter=as.numeric(args[4])

print (paste0("cnvfilter: ",cnvfilter))
print (paste0("armfilter: ",armfilter))

#heatmap_col <- c("white","#0099CC","#CC0033")
heatmap_col <- c("#FEF0D9","#FDD49E","#FDBB84","#FC8D59","#E34A33","#B30000")
col_fun <- colorRampPalette(heatmap_col)(100)

plotFun <- function (data,outdir,filter=FALSE) {
    dir.create (outdir)
    arm <- read.table (data,header = T, sep='\t')
    if (filter) {
        arm$rowsum <- rowSums(arm[,2:ncol(arm)])
        num <- unname (quantile (arm$rowsum,cnvfilter))
        arm <- filter (arm,rowsum > num ) %>% select (-rowsum)
        rownames(arm) <- arm[,1]
        arm <- arm[,-1]
        ht<-draw_pheatmap(arm, scale  = "none",cluster_rows = T,cluster_cols = T,heat_col =col_fun)
        save_pheatmap(ht, prefix = paste0("filter.",cnvfilter,".cnv_state"), outdir = outdir)
    }else{
        rownames(arm) <- arm[,1]
        arm <- arm[,-1]
        ht<-draw_pheatmap(arm, scale  = "none",cluster_rows = T,cluster_cols = T,heat_col =col_fun)
        save_pheatmap(ht, prefix = "all.cnv_state", outdir = outdir)
    }

}


dat <- file.path(input,"cnv_stats.txt")

## all 
plotFun (dat,dio)
## filter  删除突变程度较弱长短臂（即删除小于突变总和50% 的长短臂）
plotFun (dat,dio,filter = TRUE)


## get table
#data  <- read.table ("/SGRNJ03/PiplineTest01/Software_test/zhengzhanye/infercnv.test/test.heatmap/result/03.tree.new/cnv_stats.txt",header = T,stringsAsFactor = F)
print (dat)
data <- read.table (dat,header = T,stringsAsFactor = F, sep='\t')
print (head(data))
write.table(data, file.path(input,"heatmap/cnv_stat.xls"), sep='\t', row.names=F)

##all
all <- melt (data)
colnames (all) <- c("arm","type","value")
arm.value <- unname (quantile (all$value,armfilter))
filt.all <- filter (all,value>arm.value) %>% dcast (arm~type)



pp <- list()

for (i in c(2:ncol(filt.all))) {
    id=colnames (filt.all)[i]
    hcl <- filt.all[,c(1,i)] %>% drop_na() 
    hcl <- hcl[order (hcl[,2],decreasing = T),] 
    rownames (hcl) <- c(1:nrow(hcl))
    hcl <- t(hcl) %>% as.data.frame()
    select.arm <- hcl["arm",] %>%as.data.frame() 
    if (ncol (select.arm) ==1) {
        colnames (select.arm) = "X1"
    }
    select.arm <- data.frame(id,select.arm)
    pp[[id]] <- select.arm
}

rs <- do.call(rbind.fill,pp)
rs[is.na(rs)] <- ''
write.table (rs,file.path(dio,paste0("cluster.filter.",armfilter,".state.xls")),quote= F,row.names =F,col.names = F,sep="\t")
