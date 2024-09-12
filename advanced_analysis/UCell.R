## UCell function
getAnnData <- function(h5ad_path){
    suppressMessages({
        library('reticulate')
    })

    sc <- import("scanpy")
    adata <- sc$read(h5ad_path, backed = 'r')
    return(adata)
}


getSeuratMatrix <- function(adata, layer){
    suppressMessages({
        library('Seurat')
        library('Matrix')
    })

    if (layer == 'normalised') {
        layer_data <- adata$layers['normalised']
    } else if (layer == 'filtered') {
        layer_data <- adata$layers['filtered']
    } else
        stop("Layer does not exist")

    mtx_data <- sparseMatrix(
        j = layer_data@j,
        p = layer_data@p,
        x = layer_data@x,
        dims = layer_data@Dim,
        index1 = F,
        repr = "C"
    )
    rownames(mtx_data) <- rownames(adata$obs)
    colnames(mtx_data) <- rownames(adata$var)

    PRO <- CreateSeuratObject(counts = t(mtx_data), project = "auc")
    return(PRO)
}


getUCellScore <- function(PRO, gset, method){
    suppressMessages({
        library('UCell')
        library('Seurat')
    })

    marker <- list()
    marker[['gset']] <- unlist(gset)

    if (method == 'UCell') {

        PRO <- AddModuleScore_UCell(PRO, features = marker)
        result_df <- PRO@meta.data[, 'gset_UCell', drop = FALSE]

    } else if (method == 'AddModuleScore') {

        PRO <- AddModuleScore(PRO, features = marker, name='gset_Score')
        result_df <- PRO@meta.data[, 'gset_Score1', drop = FALSE]
    }

    return(result_df)
}

ucell <- function(h5ad_path, gset_list, layer, method){
    suppressMessages({
        library('dplyr')
    })

    result_list <- list()
    adata <- getAnnData(h5ad_path)
    PRO_raw <- getSeuratMatrix(adata, layer)
    if (class(gset_list)[1] == "matrix"){
        n_list = ncol(gset_list)
        for (i in 1:n_list){
            result_df <- getUCellScore(PRO_raw, gset_list[i,], method)
            names(result_df)[1] <- paste0("gset_Score_", as.character(i))
            result_list[[i]] <- result_df
        }
    }else{
        n_list <- length(gset_list)
        for (i in 1:n_list){
            result_df <- getUCellScore(PRO_raw, gset_list[[i]], method)
            names(result_df)[1] <- paste0("gset_Score_", as.character(i))
            result_list[[i]] <- result_df
        }
    }

    result_df <- bind_cols(result_list)
    return(result_df)
}