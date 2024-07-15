import anndata as ad
import click
import json
import os
import sys
import shutil
import subprocess
import numpy as np
import pandas as pd
from anndata import AnnData
from pathlib import Path
from scipy.sparse import csr_matrix
import scanpy as sc


def get_h5ad_path(matrix_info: dict):
    path = matrix_info['path']
    path_type = matrix_info['path_type']

    destination = "sample.h5ad"

    if path_type == 'LOCAL':
        try:
            shutil.copy(path, destination)
        except:
            raise Exception("Sorry, LOCAL file copy failed.")
    elif path_type == 'S3':
        try:
            pass
        except:
            raise Exception("Sorry, S3 file copy failed.")
    elif path_type == 'TOS':
        try:
            pass
        except:
            raise Exception("Sorry, TOS file copy failed.")
    else:
        return None

    return destination


def fix_var_names(adata: AnnData):
    gex_rows = list(
        map(lambda x: x.replace('_', '-'), adata.var.index)
    )
    adata.var.index = gex_rows

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    return adata


def fix_obs_names(adata: AnnData, name: str):
    gex_cols = list(
        map(lambda x: name + "_" + x, adata.obs.index)
    )
    adata.obs.index = gex_cols

    return adata


def get_n_pcs(adata):
    pct = adata.uns['pca']['variance_ratio'] / sum(adata.uns['pca']['variance_ratio']) * 100
    pctsum = pct.cumsum()
    for i in range(len(pct)):
        if (pctsum[i] > 90) and (pct[i] < 5):
            break
    diff = np.array(pct[:-1]) - np.array(pct[1:])
    index = [list(diff).index(i) for i in diff[diff > 0.1]]
    k = sorted(index, reverse=True)[0] + 1
    n_pcs = min(k, i)

    return n_pcs


def get_min_dist(adata):
    n_obs = adata.n_obs
    if n_obs < 1000:
        min_dist = 0.5
    elif 1000 < n_obs <= 100000:
        min_dist = 0.3
    elif 100000 < n_obs <= 300000:
        min_dist = 0.2
    else:
        min_dist = 0.1

    return min_dist


def hvg(adata, n_top_genes=2000, layer='filtered', *args, **kwargs):
    """
    Wrapper function for sc.highly_variable_genes()
    """
    sc.pp.highly_variable_genes(
        adata,
        layer=layer,
        n_top_genes=n_top_genes,
        min_disp=0.5,
        max_disp=np.inf,
        min_mean=0.0125,
        max_mean=3,
        flavor='seurat_v3',
        **kwargs
    )

    return adata


def leiden(adata, resolution=1, *args, **kwargs):
    """
    Wrapper function for sc.tl.leiden, for supporting multiple resolutions.
    """
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added='cluster',
        neighbors_key=None,
        obsp=None,
        **kwargs
    )

    return adata


def louvain(adata, resolution=1, *args, **kwargs):
    """
    Wrapper function for sc.tl.louvain, for supporting multiple resolutions.
    """
    sc.tl.louvain(
        adata,
        resolution=resolution,
        key_added='cluster',
        neighbors_key=None,
        obsp=None,
        **kwargs
    )

    return adata


def neighbors(adata, use_rep='X_pca', *args, **kwargs):
    """
    Wrapper function for sc.pp.neighbors(), for supporting multiple n_neighbors
    """
    n_pcs = get_n_pcs(adata)

    sc.pp.neighbors(
        adata,
        n_neighbors=15,
        n_pcs=n_pcs,
        use_rep=use_rep,
        key_added=None,
        **kwargs
    )

    return adata


def normalize(adata, *args, **kwargs):
    """
    Wrapper function for sc.pp.normalize_per_cell() and sc.pp.log1p(), mainly
    for supporting different ways of saving raw data.
    """
    adata.layers['filtered'] = adata.X.copy()
    sc.pp.normalize_total(
        adata,
        target_sum=1e4,
        **kwargs
    )
    sc.pp.log1p(
        adata,
        **kwargs
    )
    adata.layers['normalised'] = adata.X.copy()

    return adata


def pca(adata, *args, **kwargs):
    """
    Wrapper function for sc.pp.pca
    """
    sc.pp.pca(
        adata,
        **kwargs
    )

    return adata


def harmony(adata, batch_name, *args, **kwargs):
    """
    Wrapper function for sc.external.pp.harmony_integrate
    """
    sc.external.pp.harmony_integrate(
        adata,
        key=batch_name,
        basis='X_pca',
        adjusted_basis='X_pca_harmony',
        **kwargs
    )

    return adata


def scanorama(adata, batch_name, *args, **kwargs):
    """
    Wrapper function for sc.external.pp.scanorama_integrate
    """
    sc.external.pp.scanorama_integrate(
        adata,
        key=batch_name,
        basis='X_pca',
        adjusted_basis='X_scanorama'
    )

    return adata


def scale(adata, *args, **kwargs):
    """
    Wrapper function for sc.pp.scale
    """
    sc.pp.scale(
        adata,
        max_value=10,
        **kwargs
    )

    return adata


def tsne(adata, use_rep='X_pca', *args, **kwargs):
    """
    Wrapper function for sc.tl.tsne, for supporting named slot of tsne embeddings
    """
    sc.tl.tsne(
        adata,
        use_rep=use_rep,
        learning_rate=10,
        random_state=0,
        **kwargs
    )

    return adata


def umap(adata, *args, **kwargs):
    """
    Wrapper function for sc.tl.umap, for supporting named slot of umap embeddings
    """
    min_dist = get_min_dist(adata)

    sc.tl.umap(
        adata,
        min_dist=min_dist,
        random_state=0,
        **kwargs
    )

    return adata


def rank_genes_groups(adata, *args, **kwargs):
    """
    Wrapper function for sc.tl.rank_genes_groups
    """
    sc.tl.rank_genes_groups(
        adata,
        use_raw=False,
        n_genes=None,
        rankby_abs=False,
        pts=True,
        copy=False,
        method='wilcoxon',
        corr_method='benjamini-hochberg',
        tie_correct=True,
        **kwargs
    )

    return adata


def get_mtfilter(adata):
    MITO_GENE_PERCENT_LIST = [5, 10, 15, 20, 30, 50]
    adata.var["mito"] = adata.var.index.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], percent_top=None, log1p=False, inplace=True)
    mt_df = adata.obs["pct_counts_mito"]
    mt_count = []
    for mt_gene_percent in MITO_GENE_PERCENT_LIST:
        mt_count.append(len(mt_df[mt_df > mt_gene_percent]) / len(mt_df))
    mt_loc = mt_count.index(min(mt_count, key=lambda x: abs(x - 0.05)))
    mtfilter = int(MITO_GENE_PERCENT_LIST[mt_loc])

    return mtfilter


def subset_by_sample(adata_combined: AnnData, name: str):
    adata = adata_combined[
        adata_combined.obs['Sample ID'] == name
        ].copy()
    adata.X = np.nan_to_num(adata.X)
    sc.pp.filter_genes(adata, min_cells=1, inplace=True, copy=False)

    return adata


def select_data_integration(matrix_info: dict):
    samples = matrix_info['samples']
    type = matrix_info['type']

    h5ad_path = get_h5ad_path(matrix_info)
    adata = sc.read(h5ad_path)

    if type == 'USER':
        adata = AnnData(
            X=adata.X,
            obs=adata.obs[[]],
            var=adata.var[[]],
            dtype='float32'
        )
        adata.obs['Sample ID'] = samples[0]
        adata = fix_var_names(adata)
        adata = fix_obs_names(
            adata=adata,
            name=samples[0]
        )

        return adata

    else:
        return None


def filter_cells_genes(
        adata,
        analysis: dict,
        min_genes: float = 200,
        min_counts: float = 0
):
    """
    Wrapper function for sc.pp.filter_cells() and sc.pp.filter_genes(), mainly
    for supporting arbitrary filtering
    """
    sc.pp.filter_cells(
        adata,
        min_genes=min_genes,
        inplace=True,
        copy=False
    )
    sc.pp.filter_cells(
        adata,
        min_counts=min_counts,
        inplace=True,
        copy=False
    )

    if analysis.get('gene_state') == 'recommended':
        max_genes = int(
            np.percentile(
                adata.obs["n_genes"], 100 * 0.98
            ) + 1
        )
    else:
        if analysis.get('gene_num') is not None:
            max_genes = analysis.get('gene_num')[1]
        else:
            max_genes = 5000

    if analysis.get('umi_state') == 'recommended':
        max_counts = int(
            np.percentile(
                adata.obs["n_counts"], 100 * 0.98
            ) + 1
        )
    else:
        if analysis.get('umi_num') is not None:
            max_counts = analysis.get('umi_num')[1]
        else:
            max_counts = 30000

    sc.pp.filter_cells(
        adata,
        max_genes=max_genes,
        inplace=True,
        copy=False
    )
    sc.pp.filter_cells(
        adata,
        max_counts=max_counts,
        inplace=True,
        copy=False
    )
    sc.pp.filter_genes(
        adata,
        min_cells=5,
        inplace=True,
        copy=False
    )

    return adata


@click.group()
def cli():
    pass


@cli.command("save_h5ad_pipeline")
@click.argument('analysis_config', type=click.File('r'))
def save_h5ad_pipeline(analysis_config):
    analysis_config = json.load(analysis_config)
    matrix_info = analysis_config.get('matrix_info')
    species = analysis_config.get('species')

    adatas = {}
    for i in matrix_info:
        adata = select_data_integration(i)
        if adata:
            adatas[i.get('name')] = adata
    adata = ad.concat(adatas, axis=0, join='outer')

    del adatas
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    adata.obs['Species'] = pd.Categorical(
        [species] * adata.n_obs
    )

    if not isinstance(adata.X, csr_matrix):
        adata.X = csr_matrix(adata.X, shape=adata.shape)
    result_file = 'result.h5ad'
    adata.write(result_file, compression='lzf')

    with open('n_obs.txt', 'w') as f:
        f.write(str(adata.n_obs))


@cli.command("filter_pipeline")
@click.argument('analysis_config', type=click.File('r'))
def filter_pipeline(analysis_config):
    analysis_config = json.load(analysis_config)
    analysis = analysis_config.get('sc_parameters')

    result_path = 'result.h5ad'
    adata_combined = sc.read(result_path)

    if analysis.get('gene_state') == 'recommended':
        min_genes = 200
    else:
        if analysis.get('gene_num') is not None:
            min_genes = analysis.get('gene_num')[0]
        else:
            min_genes = 200

    if analysis.get('umi_state') == 'recommended':
        min_counts = 0
    else:
        if analysis.get('umi_num') is not None:
            min_counts = analysis.get('umi_num')[0]
        else:
            min_counts = 0

    adatas = {}
    mtfilters = []
    for name in adata_combined.obs['Sample ID'].cat.categories:
        adata = subset_by_sample(
            adata_combined,
            name
        )
        adata = filter_cells_genes(
            adata=adata,
            analysis=analysis,
            min_genes=min_genes,
            min_counts=min_counts
        )
        adatas[name] = adata
        mtfilters.append(get_mtfilter(adata))
    adata = ad.concat(adatas, axis=0, join='outer')
    del adatas

    if analysis.get('mitochondrial_state') == 'recommended':
        mtfilters_df = pd.DataFrame(mtfilters)
        pct_counts_mito = mtfilters_df.mode().sort_values(by=0, ascending=False)[0].iloc[0]
    else:
        if analysis.get('mtfilter') is not None:
            pct_counts_mito = analysis.get('mtfilter')
        else:
            pct_counts_mito = 50

    adata.var["mito"] = adata.var_names.str.contains("^[Mm][Tt]-")
    sc.pp.calculate_qc_metrics(
        adata,
        expr_type='counts',
        var_type='genes',
        qc_vars=["mito"],
        percent_top=None,
        use_raw=False,
        inplace=True,
        log1p=False
    )
    adata = adata[adata.obs['pct_counts_mito'] < pct_counts_mito, :]
    del adata.obs['n_genes']
    del adata.obs['n_counts']

    result_file = 'result.h5ad'
    adata.write(result_file, compression='lzf')


@cli.command("scanpy_pipeline")
@click.argument('analysis_config', type=click.File('r'))
def scanpy_pipeline(analysis_config):
    analysis_config = json.load(analysis_config)
    species = analysis_config.get('species')

    analysis = analysis_config.get('sc_parameters')
    resolution = analysis.get('resolution')
    batch_remove = analysis.get('batch_effect_removal')
    batch_name = analysis.get('batch_name', 'Sample ID')
    batch_remove_method = analysis.get('batch_remove_method', 'harmony')
    n_top_genes = analysis.get('n_top_genes', 2000)
    hvg_subscale = analysis.get('scale_large_data')
    remove_cell_cycle_effects = analysis.get('remove_cell_cycle_effects')

    result_path = 'result.h5ad'
    adata = sc.read(result_path)

    if hvg_subscale:
        adata.raw = adata

    normalize(adata=adata)

    if hvg_subscale:
        data = adata.copy()

    hvg(adata=adata, n_top_genes=n_top_genes, layer='filtered')

    if species in ['Homo sapiens', 'Mus musculus']:
        speciesx = 'human' if species == 'Homo sapiens' else 'mouse'
        cc_genes_file = os.path.join(sys.path[0], f'data/{speciesx}_cellcycle_genes.txt')
        cc_genes = pd.read_csv(cc_genes_file, sep="\t", index_col=0, header=None)
        s_genes = cc_genes.loc['s.genes'].str.split(',').to_list()[0]
        g2m_genes = cc_genes.loc['g2m.genes'].str.split(',').to_list()[0]
        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
        if remove_cell_cycle_effects:
            sc.pp.regress_out(adata, ['S_score', 'G2M_score'], 4)

    if hvg_subscale:
        adata = adata[:, adata.var.highly_variable]

    scale(adata=adata)
    pca(adata=adata)
    rep = 'X_pca'
    if batch_remove:
        if batch_remove_method == 'harmony':
            harmony(adata=adata, batch_name=batch_name)
            rep = 'X_pca_harmony'
        if batch_remove_method == 'scanorama':
            scanorama(adata=adata, batch_name=batch_name)
            rep = 'X_scanorama'
    neighbors(adata=adata, use_rep=rep)
    louvain(adata=adata, resolution=resolution)
    umap(adata=adata)
    tsne(adata=adata, use_rep=rep)

    cluster_dict = {
        str(num): str(num + 1)
        for num in range(adata.obs['cluster'].astype('int').max() + 1)
    }
    adata.obs['cluster'] = adata.obs['cluster'].map(cluster_dict).astype('category')
    cluster_order = [i for i in cluster_dict.values()]
    adata.obs['cluster'].cat.set_categories(cluster_order, inplace=True)

    if hvg_subscale:
        adata = adata.raw.to_adata()
        adata.X = data.layers["normalised"]
        adata.layers['normalised'] = data.layers['normalised']
        adata.layers['filtered'] = data.layers['filtered']
    else:
        adata.X = adata.layers["normalised"]

    key_added = '_'.join(['rank_genes_groups', 'cluster'])
    rank_genes_groups(adata, groupby='cluster', layer='normalised', key_added=key_added)

    result_file = 'result.h5ad'
    adata.write(result_file, compression='lzf')