suppressMessages({
    library(Seurat)
    library(celda)
    library(argparser)
    library(ggplot2)
    library(reticulate)
    library(SingleCellExperiment)
})

argv <- arg_parser('Contamination_removal')
argv <- add_argument(argv, "--h5", help = "h5 matrix file")
argv <- add_argument(argv, "--spname", help = "the samples name")
argv <- add_argument(argv, "--outdir", help = "outdir")
argv <- parse_args(argv)

seurat2anndata <- function(seurat, outFile) {
    seurat@meta.data$cluster <- Idents(seurat)

    sc <- import("scanpy")
    adata <- sc$AnnData(
        X = Matrix::t(GetAssayData(seurat, slot='counts')),
        obs = seurat[["contamination"]],
        var = GetAssay(seurat)[[]]
    )

    tryCatch(
    { adata$obsm$update(X_umap = Embeddings(seurat, "umap")) },
    error = function(e) { print("NO UMAP EMBEDDINGS") }
    )

    adata$write(outFile)
}

dir.create(argv$outdir)
RAW1 <- Read10X_h5(argv$h5)
PRO <- CreateSeuratObject(counts = RAW1, project = argv$spname, min.cells = 5)
mt.genes <- rownames(PRO)[grep("^MT-", rownames(PRO), ignore.case = TRUE)]
PRO[["percent.mt"]] <- PercentageFeatureSet(PRO, features = mt.genes)
PRO <- subset(
    PRO,
    subset = nFeature_RNA > 200 &
        percent.mt < 50 &
        nCount_RNA > 0
)
PRO <- NormalizeData(object = PRO)
PRO <- FindVariableFeatures(PRO, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(PRO)
PRO <- ScaleData(PRO, features = all.genes)
PRO <- RunPCA(PRO, features = VariableFeatures(object = PRO))
PRO <- RunUMAP(PRO, reduction = "pca", dims = 1:20)
PRO <- FindNeighbors(PRO, reduction = "pca", dims = 1:20)
PRO <- FindClusters(PRO, resolution = 0.8)
PRO_sce <- as.SingleCellExperiment(PRO)

PRO_sce <- decontX(PRO_sce, z=PRO@meta.data$seurat_clusters)
counts.dcon <- round(PRO_sce@assays@data$decontXcounts)

data.con <- data.frame(cellID = colnames(PRO), ContaminationValue = PRO_sce$decontX_contamination)
colnames(data.con) <- c("barcode", "ContaminationValue")
write.table(data.con, paste(argv$outdir, "/", argv$spname, "_contaminationValue_per_barcode.xls", sep = ""), sep = "\t", quote = F, row.names = F)
PRO@meta.data$contamination <-PRO_sce$decontX_contamination
PRO@assays$RNA@counts <- counts.dcon 
seurat2anndata(PRO, paste(argv$outdir, "/", argv$spname, "_contamination.h5ad", sep = ""))
