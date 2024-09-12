from anndata import AnnData
from typing import Literal, Union
import click
import json
import numpy as np
import pandas as pd
import scanpy as sc


def chunks(length, n):
    for i in range(0, length, n):
        if i + n < length:
            yield i, i + n
        else:
            yield i, length


# no used, for backup
def calculate_rank_genes_groups_mean(adata, groupby='cluster', layer='normalised', key='mean', *args, **kwargs):
    """
    Wrapper function for sc.tl.rank_genes_groups
    """
    genes = adata.var_names.tolist()
    df = pd.DataFrame(index=adata.obs[groupby].cat.categories)
    for i, j in chunks(adata.n_vars, 1000):
        df = df.merge(
            pd.concat(
                [
                    adata.obs[[groupby]],
                    sc.get.obs_df(
                        adata,
                        keys=genes[i:j],
                        obsm_keys=(),
                        layer=layer,
                        gene_symbols=None,
                        use_raw=False
                    )
                ],
                axis=1
            ).groupby(
                by=groupby
            ).mean(),
            left_index=True,
            right_index=True,
        )

    return df.T


def mean_expression(
        adata: AnnData,
        groups: Union[list[str], Literal["all"]],
        groupby: str,
        layer: str,
        pseudocount: float = 1.0
):
    from scanpy._utils import select_groups
    from scanpy.preprocessing._utils import _get_mean_var

    groups_order, groups_masks_obs = select_groups(
        adata, groups, groupby
    )

    X = adata.layers[layer]
    var_names = adata.var_names

    n_genes = X.shape[1]
    n_groups = groups_masks_obs.shape[0]

    means = np.zeros((n_groups, n_genes))
    X = np.expm1(X)

    for imask, mask in enumerate(groups_masks_obs):
        X_mask = X[mask]
        means[imask], _ = _get_mean_var(X_mask)

    means = pd.DataFrame(
        np.log2(means + pseudocount),
        index=groups_order,
        columns=var_names
    ).T

    return means


def get_means_df(adata, selected_group=None, type='Seurat', *args, **kwargs):
    means_res = {}
    for key, layer in zip(['mean'], ['normalised']):
        if layer in adata.layers:
            if type == 'Seurat':
                means_res[key] = mean_expression(adata, groups='all', groupby='selected_group', layer=layer)[
                    selected_group
                ]
            else:
                means_res[key] = calculate_rank_genes_groups_mean(adata, groupby='selected_group', layer=layer)[
                    selected_group
                ]

    means_df = pd.concat(means_res, axis=1)
    means_df.reset_index(inplace=True)
    means_df = means_df.rename(columns={'index': 'names'})

    return means_df


def calculate_differential_expression(h5ad_path, query, layer='normalised'):
    adata = sc.read(h5ad_path)
    query1, query2 = query

    adata.obs['selected_group'] = None
    adata.obs['selected_group'] = adata.obs['selected_group'].mask(
        adata.obs[query1.keys()].isin(query1).all(axis=1),
        'selected_group1'
    )
    adata.obs['selected_group'] = adata.obs['selected_group'].mask(
        adata.obs[query2.keys()].isin(query2).all(axis=1),
        'selected_group2'
    ).astype('category')

    adata = adata[
        adata.obs['selected_group'].isin(
            ['selected_group1', 'selected_group2']
        )
    ]

    sc.tl.rank_genes_groups(
        adata=adata,
        groupby='selected_group',
        use_raw=False,
        groups=['selected_group1'],
        reference='rest',
        n_genes=None,
        rankby_abs=False,
        pts=True,
        key_added=None,
        method='wilcoxon',
        corr_method='benjamini-hochberg',
        tie_correct=True,
        layer=layer
    )
    df_group1 = sc.get.rank_genes_groups_df(
        adata,
        group='selected_group1'
    )
    df_means1 = get_means_df(adata, selected_group='selected_group1', type='Seurat')
    df_group1 = pd.merge(df_group1, df_means1, how='left', on='names')
    df_group1 = df_group1.astype(
        {
            i: 'float64'
            for i in df_group1.columns
            if i in ['scores', 'logfoldchanges', 'mean']
        }
    )
    df_group1 = df_group1.fillna(value=0)

    sc.tl.rank_genes_groups(
        adata=adata,
        groupby='selected_group',
        use_raw=False,
        groups=['selected_group2'],
        reference='rest',
        n_genes=None,
        rankby_abs=False,
        pts=True,
        key_added=None,
        method='wilcoxon',
        corr_method='benjamini-hochberg',
        tie_correct=True,
        layer=layer
    )
    df_group2 = sc.get.rank_genes_groups_df(
        adata,
        group='selected_group2'
    )
    df_means2 = get_means_df(adata, selected_group='selected_group2', type='Seurat')
    df_group2 = pd.merge(df_group2, df_means2, how='left', on='names')
    df_group2 = df_group2.astype(
        {
            i: 'float64'
            for i in df_group2.columns
            if i in ['scores', 'logfoldchanges', 'mean']
        }
    )
    df_group2 = df_group2.fillna(value=0)

    return df_group1, df_group2


@click.group()
def cli():
    pass


@cli.command("calculate_DEGs_pipeline")
@click.argument('analysis_config', type=click.File('r'))
def calculate_DEGs_pipeline(analysis_config):
    analysis_config = json.load(analysis_config)

    h5ad_path = "./temp_file/result.h5ad"
    query1 = analysis_config["scanpy_parameters"]["query1"]
    query2 = analysis_config["scanpy_parameters"]["query2"]
    query = [query1, query2]
    layer = analysis_config["scanpy_parameters"]["layer"]

    result_group1, result_group2 = calculate_differential_expression(h5ad_path, query, layer)
    result_group1.to_csv('./DEGs_result_group1.tsv', sep='\t', index=False)
    result_group2.to_csv('./DEGs_result_group2.tsv', sep='\t', index=False)


@cli.command("filter_DEGs_pipeline")
@click.argument('analysis_config', type=click.File('r'))
def filter_DEGs_pipeline(analysis_config):
    analysis_config = json.load(analysis_config)

    scores = analysis_config["filter_parameters"].get("scores", None)
    logfoldchanges = analysis_config["filter_parameters"].get("logfoldchanges", None)
    pvals = analysis_config["filter_parameters"].get("pvals", None)
    pvals_adj = analysis_config["filter_parameters"].get("pvals_adj", None)
    pct_nz_group = analysis_config["filter_parameters"].get("pct_nz_group", None)

    result_path = "./DEGs_result.tsv"
    df_group = pd.read_csv(result_path, sep='\t')

    if scores:
        if scores[0]:
            df_group = df_group[df_group["scores"] > scores[0]]
        if scores[1]:
            df_group = df_group[df_group["scores"] < scores[1]]

    if logfoldchanges:
        df_group = df_group[df_group["logfoldchanges"] > logfoldchanges]

    if pvals:
        df_group = df_group[df_group["pvals"] < pvals]

    if pvals_adj:
        df_group = df_group[df_group["pvals_adj"] < pvals_adj]

    if pct_nz_group:
        df_group = df_group[df_group["pct_nz_group"] > pct_nz_group]

    df_group.to_csv('./DEGs_filter_result.tsv', sep='\t', index=False)


@cli.command("plot_DEGs_pipeline")
@click.argument('analysis_config', type=click.File('r'))
def plot_DEGs_pipeline(analysis_config):
    analysis_config = json.load(analysis_config)

    go_path = "./enrich_go_result.tsv"
    kg_path = "./enrich_kg_result.tsv"
    rt_path = "./enrich_rt_result.tsv"
    wk_path = "./enrich_wk_result.tsv"
    re_go = pd.read_csv(go_path, sep='\t')
    re_kg = pd.read_csv(kg_path, sep='\t')
    re_rt = pd.read_csv(rt_path, sep='\t')
    re_wk = pd.read_csv(wk_path, sep='\t')

    re_go['logp'] = -np.log10(re_go['p.adjust'])
    re_go['GeneRatio2'] = re_go['GeneRatio']
    re_go['GeneRatio'] = re_go['GeneRatio'].apply(pd.eval)
    # re_go = re_go[['ONTOLOGY', 'Description', 'GeneRatio', 'logp', 'Count']]
    re_go.to_csv('./plot_go_result.tsv', sep='\t', index=False)

    re_kg['logp'] = -np.log10(re_kg['p.adjust'])
    re_kg['GeneRatio2'] = re_kg['GeneRatio']
    re_kg['GeneRatio'] = re_kg['GeneRatio'].apply(pd.eval)
    # re_kg = re_kg[['Description', 'GeneRatio', 'logp', 'Count']]
    re_kg.to_csv('./plot_kg_result.tsv', sep='\t', index=False)

    re_rt['logp'] = -np.log10(re_rt['p.adjust'])
    re_rt['GeneRatio2'] = re_rt['GeneRatio']
    re_rt['GeneRatio'] = re_rt['GeneRatio'].apply(pd.eval)
    # re_rt = re_rt[['Description', 'GeneRatio', 'logp', 'Count']]
    re_rt.to_csv('./plot_rt_result.tsv', sep='\t', index=False)

    re_wk['logp'] = -np.log10(re_wk['p.adjust'])
    re_wk['GeneRatio2'] = re_wk['GeneRatio']
    re_wk['GeneRatio'] = re_wk['GeneRatio'].apply(pd.eval)
    # re_wk = re_wk[['Description', 'GeneRatio', 'logp', 'Count']]
    re_wk.to_csv('./plot_wk_result.tsv', sep='\t', index=False)


if __name__ == '__main__':
    cli()
