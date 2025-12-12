# library(pheatmap)
# library(argparser)

# argv <- arg_parser('')
# argv <- add_argument(argv,"--inputdir",help = "path of cnv_stats.txt file")
# argv <- add_argument(argv,"--outdir",help = "output dir")
# argv <- parse_args(argv)

# dir <- argv$inputdir
# outdir <- argv$outdir
# dir.create (outdir,showWarnings=T)
args<-commandArgs(T)
inputdir=args[1]
# outdir=args[2]
# dir.create (outdir,showWarnings=T)

data <- read.table (file.path(inputdir,"infercnv.observations.txt"),stringsAsFactors=F)
dd <- t(data)
#write.table (dd,file.path(inputdir,"infercnv.observations.covert.txt"),sep=' ',quote=F,row.names=T,col.names=T)


ref <- read.table (file.path(inputdir,"infercnv.references.txt"),stringsAsFactors=F)
tref <- t(ref)
#write.table (tref,file.path(inputdir,"infercnv.references.covert.txt"),sep=' ',quote=F,row.names=T,col.names=T)

rs <- rbind (dd,tref)
write.table (rs,file.path(inputdir,"infercnv.observations.references.covert.txt"),sep=' ',quote=F,row.names=T,col.names=T)