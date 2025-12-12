import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import click
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from anndata import AnnData
from collections import Counter
from scipy.sparse import csr_matrix


def scrublet(adata, *args, **kwargs):
    """
    Wrapper function for sce.pp.scrublet()
    """
    try:
        sce.pp.scrublet(
            adata,
            expected_doublet_rate=0.05,
            **kwargs
        )
    except:
        sce.pp.scrublet(
            adata,
            expected_doublet_rate=0.25,
            **kwargs
        )

    if 'threshold' not in adata.uns['scrublet']:
        sce.pp.scrublet(
            adata,
            expected_doublet_rate=0.25,
            **kwargs
        )

    return adata


def rm_cont(matrix, prefix, outdir):
    """
    run decontx and output processed matrix
    """
    decontx_cmd = f"Rscript {sys.path[0]}/decontx.R \
                  --h5 '{matrix}' \
                  --spname '{prefix}' \
                  --outdir '{outdir}'"
    subprocess.run(decontx_cmd, shell=True, check=True)


def dict2array(dict_info):
    result = dict_info.items()
    data = list(result)

    if all(isinstance(value, float) for value in dict_info.values()):
        dt = [('Sample ID', 'U64'), ('info', '<f8')]
    elif all(isinstance(value, bool) for value in dict_info.values()):
        dt = [('Sample ID', 'U64'), ('info', '?')]
    else:
        dt = [('Sample ID', 'U64'), ('info', 'U64')]
    data = np.array(data, dtype=dt)

    return data


def _collect_datasets(dsets: dict, group: h5py.Group):
    for k, v in group.items():
        if isinstance(v, h5py.Dataset):
            dsets[k] = v[()]
        else:
            _collect_datasets(dsets, v)


def _read_v3_10x_h5(filename, *, start=None):
    """
    Read hdf5 file from Cell Ranger v3 or later versions.
    """
    with h5py.File(str(filename), "r") as f:
        try:
            dsets = {}
            _collect_datasets(dsets, f["matrix"])
            from scipy.sparse import csr_matrix
            M, N = dsets["shape"]
            data = dsets["data"]
            if dsets["data"].dtype == np.dtype("int32"):
                data = dsets["data"].view("float32")
                data[:] = dsets["data"]
            matrix = csr_matrix(
                (data, dsets["indices"], dsets["indptr"]),
                shape=(N, M),
            )
            obs_dict = {"obs_names": dsets["barcodes"].astype(str)}
            var_dict = {"var_names": dsets["name"].astype(str)}
            if "gene_id" not in dsets:
                var_dict["gene_ids"] = dsets["id"].astype(str)
            else:
                var_dict.update(
                    {
                        "gene_ids": dsets["gene_id"].astype(str),
                        "probe_ids": dsets["id"].astype(str),
                    }
                )
            var_dict["feature_types"] = dsets["feature_type"].astype(str)
            if "filtered_barcodes" in f["matrix"]:
                obs_dict["filtered_barcodes"] = dsets["filtered_barcodes"].astype(bool)
            if "features" in f["matrix"]:
                var_dict.update(
                    (
                        feature_metadata_name,
                        dsets[feature_metadata_name].astype(
                            bool if feature_metadata_item.dtype.kind == "b" else str
                        ),
                    )
                    for feature_metadata_name, feature_metadata_item in f["matrix"][
                        "features"
                    ].items()
                    if isinstance(feature_metadata_item, h5py.Dataset)
                    and feature_metadata_name
                    not in [
                        "name",
                        "feature_type",
                        "id",
                        "gene_id",
                        "_all_tag_keys",
                    ]
                )
            else:
                raise ValueError("10x h5 has no features group")
            adata = AnnData(
                matrix,
                obs=obs_dict,
                var=var_dict,
            )
            return adata
        except KeyError:
            raise Exception("File is missing one or more required datasets.")


def _write_10x_h5(adata, output):
    """
    write hdf5 file from Cell Ranger v3 or later versions.
    """
    matrix = adata.X
    barcode = adata.obs_names
    genes = adata.var_names
    with h5py.File(output, 'w') as f:
        group = f.create_group("matrix")
        group.create_dataset('barcodes', data=np.array(barcode))
        group.create_dataset('data', data=matrix.data)
        group.create_dataset('indices', data=matrix.indices)
        group.create_dataset('indptr', data=matrix.indptr)
        group.create_dataset('shape', data=matrix.T.shape)
        features = group.create_group("features")
        features.create_dataset('_all_tag_keys', data=np.array(['genome'], dtype='|O'))
        features.create_dataset('feature_type', data=np.array(len(genes) * ['Gene Expression'], dtype='|O'))
        features.create_dataset('genome', data=np.array(len(genes) * ['Unknown'], dtype='|O'))
        features.create_dataset('id', data=np.array(genes))
        features.create_dataset('name', data=np.array(genes))


def subset_by_sample(adata_combined, name):
    adata = adata_combined[adata_combined.obs['Sample ID'] == name].copy()
    adata.X = np.nan_to_num(adata.X)
    sc.pp.filter_genes(adata, min_cells=1, inplace=True, copy=False)

    return adata


def get_all_celltype_tag_from_obs(obs_keys: list):
    celltype_tag = []
    if 'annot_sub2' in obs_keys:
        celltype_tag.append('annot_sub2')
    if 'annot_sub' in obs_keys:
        celltype_tag.append('annot_sub')
    if 'annot_full' in obs_keys:
        celltype_tag.append('annot_full')

    return celltype_tag


@click.group()
def cli():
    pass


@cli.command("doublet_removal")
@click.argument('analysis_config', type=click.File('r'))
def doublet_removal(analysis_config):
    analysis_config = json.load(analysis_config)
    analysis = analysis_config.get('sc_parameters')

    result_path = './temp_file/result.h5ad'
    adata_combined = sc.read(result_path)

    adatas = {}
    doublet_removal_info = {}
    threshold_info = {}
    for name in adata_combined.obs['Sample ID'].cat.categories:
        print(f'--- {name} doublet_detection task in progress ---')
        try:
            adata = subset_by_sample(adata_combined, name)
            adata = scrublet(adata)

            # fix
            if adata.obs['predicted_doublet'].isnull().any():
                adata = adata[adata.obs['predicted_doublet'].notnull()]
                adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype(bool)

            doublet_scores_sim = adata.uns['scrublet']['doublet_scores_sim']
            sim1, sim2, *_ = np.split(doublet_scores_sim, [adata.n_obs, 2 * adata.n_obs])
            adata.obs['sim1'], adata.obs['sim2'] = sim1, sim2
            adatas[name] = adata

            doublet_removal_info[name] = adata.obs['predicted_doublet'].sum() / adata.n_obs
            threshold_info[name] = adata.uns['scrublet']['threshold']
        except:
            raise Exception(f'--- {name} doublet_detection task failed ---')

    adata = ad.concat(adatas, axis=0, join='outer')
    del adatas

    # for detection
    adata_detection = adata.copy()
    adata_detection.uns['doublet_removal_info'] = dict2array(doublet_removal_info)
    adata_detection.uns['threshold_info'] = dict2array(threshold_info)

    result_file = 'doublet_detection.h5ad'
    adata_detection.write(result_file, compression='lzf')

    # for result
    if analysis.get('doublet_detection'):
        adata_removal = adata[np.invert(adata.obs["predicted_doublet"]), :].copy()
    else:
        adata_removal = adata_combined.copy()

    for key in ['doublet_score', 'predicted_doublet', 'sim1', 'sim2']:
        if key in adata_removal.obs.keys():
            del adata_removal.obs[key]

    result_file = 'result.h5ad'
    adata_removal.write(result_file, compression='lzf')


@cli.command("contamination_removal")
@click.argument('analysis_config', type=click.File('r'))
def contamination_removal(analysis_config):
    analysis_config = json.load(analysis_config)
    analysis = analysis_config.get('sc_parameters')

    result_path = './temp_file/result.h5ad'
    adata_combined = sc.read(result_path)

    doublet_dir = Path("before_removal")
    doublet_dir.mkdir(mode=0o777, parents=True, exist_ok=True)

    adatas = {}
    contamination_removal_info = {}

    clu_num = Counter(adata_combined.obs['Sample ID'])
    for name in adata_combined.obs['Sample ID'].cat.categories:
        print(f'--- {name} contamination_detection task in progress ---')
        try:
            adata = subset_by_sample(adata_combined, name)
            before_removal_h5 = f"./before_removal/{name}.h5"
            _write_10x_h5(adata, before_removal_h5)
            if clu_num[name] < 100:
                # pass if cell number lower than 100
                pass
            else:
                rm_cont(before_removal_h5, name, "./contamination_removal")
                h5ad_path = f"./contamination_removal/{name}_contamination.h5ad"
                adata = sc.read(h5ad_path)
                if not isinstance(adata.X, csr_matrix):
                    adata.X = csr_matrix(adata.X, shape=adata.shape)
            adatas[name] = adata

            cont_per_bc = Path(f"./contamination_removal/{name}_contaminationValue_per_barcode.xls")
            if cont_per_bc.exists():
                predict = pd.read_table(cont_per_bc)
                contamination_removal_info[name] = predict['ContaminationValue'].median()
            else:
                contamination_removal_info[name] = 0

        except:
            raise Exception(f'--- {name} contamination_detection task failed ---')

    # for detection
    adata = ad.concat(adatas, axis=0, join='outer', label="Sample ID", merge="same")
    del adatas
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.uns['contamination_removal_info'] = dict2array(contamination_removal_info)

    if not isinstance(adata.X, csr_matrix):
        adata.X = csr_matrix(adata.X, shape=adata.shape)
    result_file = 'contamination_detection.h5ad'
    adata.write(result_file, compression='lzf')

    # for result
    if analysis.get('contamination_removal'):
        matrix_names = adata_combined.obs['Sample ID'].cat.categories

        if analysis.get('contamination_state') == 'default' or analysis.get('contamination_state') is None:
            adata_removal = adata.copy()

        if analysis.get('contamination_state') == 'based_on_single_samples':
            adatas = {}
            for name in matrix_names:
                if name not in analysis.get('based_on_single_samples'):
                    tmp = _read_v3_10x_h5(f"./before_removal/{name}.h5")
                else:
                    tmp = subset_by_sample(adata, name)
                adatas[name] = tmp
            adata_removal = ad.concat(adatas, axis=0, join='outer', label="Sample ID")
            del adata_removal.obs['contamination']
            del adatas

        adata_removal.var_names_make_unique()
        adata_removal.obs_names_make_unique()

        if not isinstance(adata_removal.X, csr_matrix):
            adata_removal.X = csr_matrix(adata_removal.X, shape=adata_removal.shape)

        meta_type = analysis_config.get('meta_type')
        meta_data = analysis_config.get('meta_data')
        annot_type = analysis_config.get('annot_type')

        if all(
                [
                    meta_type == "JSON",
                    meta_data
                ]
        ):
            metadata = pd.DataFrame.from_dict(meta_data).rename(columns={'SampleID': 'Sample ID'}).replace('', np.nan)
            df = adata_removal.obs.merge(metadata, on='Sample ID', how='left').set_index(adata_removal.obs.index)

            # Standardized annotation types for public datasets only
            if annot_type == "standardized":
                annot_key = get_all_celltype_tag_from_obs(adata_combined.obs_keys())
                df = df.merge(
                    adata_combined.obs[annot_key],
                    left_index=True,
                    right_index=True
                )
            adata_removal.obs = df
    else:
        adata_removal = adata_combined.copy()

    result_file = 'result.h5ad'
    adata_removal.write(result_file, compression='lzf')


if __name__ == '__main__':
    cli()
