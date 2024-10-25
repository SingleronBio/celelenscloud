import os
import subprocess
import pandas as pd
from tools import createDir
from v4 import run_v4

def run_pipeline(kwargs):
    print(kwargs)
    outdir = kwargs['outdir']
    prefix = kwargs['prefix']
    _CURDIR = os.path.dirname(os.path.abspath(__file__))
    result_dir = f'{outdir}/result'
    plot_dir = f'{outdir}/plot'
    createDir(result_dir)
    createDir(plot_dir)
    matrix = get_matrix_deg(
        Rapp = f'{_CURDIR}/R/get_matrix_deg.R',
        rds_path = kwargs['rds'],
        prefix = prefix,
        outdir = result_dir,
        species = kwargs['species'],
        ncell = kwargs['ncell'],
        sample = kwargs['sample'],
        group = kwargs['group'],
        degs_file_path = kwargs['degs'],
        diffcluster = kwargs['diffcluster'],
        diffgroup = kwargs['diffgroup'],
        logfc = kwargs['logfc'],
        p_adj_cutoff = kwargs['p_adj_cutoff'],
        p_cutoff = kwargs['p_cutoff']
    )
    run_cellphonedb(
        cpdb_file_path = f'{_CURDIR}/database/cellphonedb.zip',
        meta_file_path = f'{result_dir}/{prefix}_meta.tsv',
        counts_file_path = f'{result_dir}/{prefix}_count.tsv',
        outdir = result_dir,
        prefix = prefix,
        microenvs_file_path = kwargs['microenvs_file_path'],
        degs_file_path = kwargs['degs'],
        method = kwargs['method'],
        version = kwargs['cpdb_version'],
    )
    plot_cellphonedb(
        Rapp_dir = f'{_CURDIR}/R',
        outdir = plot_dir,
        prefix = prefix,
        cpdb_result = result_dir,
        database_dir = f'{_CURDIR}/database',
        meta = f'{result_dir}/{prefix}_meta.tsv',
        rds = kwargs['rds'],
    )
    os.system(f'''awk -F "\t" "NR>1 && NR<32 {'{'}print $2{'}'}" {plot_dir}/significant_means.xls > {result_dir}/rank30_rows.xls''')
    # extract rank30 pairs
    rank30 = pd.read_csv(f'{result_dir}/rank30_rows.xls', sep='\t', header=None)
    rank30[1].to_csv(f'{result_dir}/rank30_rows.txt', index=None, header=None)
    print('work done !')


def get_matrix_deg(
    Rapp: str,
    rds_path: str,
    prefix: str='cellphonedb',
    outdir: str='./',
    species: str='homo_sapiens',
    version: str='V3',
    platform: str='RNA',
    ncell: int=500,
    subcluster: str='all',
    sample: str='all',
    group: str='all',
    degs_file_path: str='F',
    diffcluster: str='F',
    diffgroup: str='F',
    logfc: float=0.2,
    p_adj_cutoff: float=0.01,
    p_cutoff: float=1
):
    cmd = f"""
Rscript {Rapp} \
    --rds {rds_path} \
    --prefix {prefix} \
    --outdir {outdir} \
    --species {species} \
    --version {version} \
    --platform {platform} \
    --ncell {ncell} \
    --subcluster {subcluster} \
    --sample {sample} \
    --group {group} \
    --degs_file_path {degs_file_path} \
    --diffcluster {diffcluster} \
    --diffgroup {diffgroup} \
    --logfc {logfc} \
    --p_adj_cutoff {p_adj_cutoff} \
    --p_cutoff {p_cutoff} \
"""
    print(cmd)
    return_code = subprocess.run(cmd, shell=True)
    return return_code


def run_cellphonedb(
    cpdb_file_path: str,
    meta_file_path: str,
    counts_file_path: str,
    outdir: str,
    prefix: str,
    microenvs_file_path: str=None,
    degs_file_path: str=None,
    method: str='statistical',
    version: str='v4',
    threads: int=3,
    threshold: float=0.1,
    iterations: int=1000
):
    microenvs_file_path = None if microenvs_file_path == 'F' else microenvs_file_path
    degs_file_path = None if degs_file_path == 'F' else degs_file_path
    if version == 'v4':
        result = run_v4(
            cpdb_file_path = cpdb_file_path,
            meta_file_path = meta_file_path,
            counts_file_path = counts_file_path,
            counts_data = "ensembl",
            output_path = outdir,
            output_suffix = prefix,
            microenvs_file_path = microenvs_file_path,
            method = method,
            degs_file_path = degs_file_path
        )
        cmd = f"""
mv {outdir}/*analysis_deconvoluted_{prefix}.txt {outdir}/deconvoluted.xls
mv {outdir}/*_analysis_means_*.txt  {outdir}/means.txt
mv {outdir}/*_analysis_pvalues*.txt  {outdir}/pvalues.txt
mv {outdir}/*_analysis_significant_means_*.txt  {outdir}/significant_means.txt
mv {outdir}/*_analysis_relevant_interactions*.txt {outdir}/relevant.xls
        """
    else:
        cmd = f"""
cellphonedb method statistical_analysis {meta_file_path} {counts_file_path} \
    --output-path {output_path} \
    --threads {threads} \
    --threshold {threshold} \
    --iterations {iterations} \
    --quiet
mv {outdir}/deconvoluted.txt {outdir}/deconvoluted.xls
"""
    print(cmd)
    return_code = subprocess.run(cmd, shell=True)
    return return_code


def plot_cellphonedb(
    Rapp_dir: str,
    outdir: str,
    prefix: str,
    cpdb_result: str,
    database_dir: str,
    meta: str,
    rds: str,
    newcol: str='auto',
    display_numbers: str='F',
    clusterorder: str='F'
):

    cmd = f'''
Rscript {Rapp_dir}/change_L_R.R --input {cpdb_result} --outdir {outdir} --prefix {prefix}
rm {outdir}/*.txt

Rscript {Rapp_dir}/chord.R --input {outdir} --outdir {outdir} --prefix {prefix} --rds {rds} --newcol {newcol}
Rscript {Rapp_dir}/plot.net.R --networkfile {outdir}/order/{prefix}_count_network_order.xls --outdir {outdir} --newcol {newcol} --prefix {prefix} --rds {rds}
Rscript {Rapp_dir}/plot_heatmaps.R --meta {meta} --pvalue {outdir}/pvalues.xls --outdir {outdir} --countmax F --clusterorder {clusterorder} --display_numbers {display_numbers} --count_filename {prefix}_heatmap_count.pdf --log_filename {prefix}_heatmap_log.pdf --count_network_filename count_network.xls --relevant {outdir}/relevant.xls

Rscript {Rapp_dir}/get_rows.R {outdir}/order/{prefix}_pvalues_order.xls {outdir}/order/{prefix}_relevant_order.xls {outdir}/order/{prefix}_significant_means_order.xls {database_dir}/LR/Chemokines_LR.xls {database_dir}/LR/Checkpoint_LR.txt {database_dir}/LR/Cytokine_LR.txt {database_dir}/LR/Grow_factor_LR.txt {outdir}/order {outdir}/order/{prefix}_count_network_order.xls {prefix}

Rscript {Rapp_dir}/cellphoneDB_picture.R --input {outdir}/order --outdir {outdir}/order --iterations 1000 --pvalue {outdir}/order/{prefix}_pvalues_order.xls --relevant {outdir}/order/{prefix}_relevant_order.xls --mean {outdir}/order/{prefix}_means_order.xls --prefix {prefix}

rm {outdir}/order/*celltype.txt {outdir}/order/*_rows.txt

perl {Rapp_dir}/prepare_cellphoneDB.v1.0.pl -cellphone {outdir}/significant_means.xls -pvalues {outdir}/pvalues.xls -relevant {outdir}/relevant.xls -means {outdir}/means.xls -outdir {outdir} -sample {prefix}_celltype_LR_sigmean 

mkdir -p {outdir}/{prefix}/heatmap
mv {outdir}/{prefix}_celltype_LR_sigmean.cellphone.xls {outdir}/order/{prefix}_target_LR.xls {outdir}/{prefix}
mv {outdir}/count_network.xls {outdir}/{prefix}/{prefix}_count_network.xls
mv {outdir}/{prefix}_heatmap_count.pdf {outdir}/{prefix}_heatmap_log.pdf {outdir}/{prefix}/heatmap

convert -density 400 -background white -alpha remove {outdir}/{prefix}/heatmap/{prefix}_heatmap_count.pdf {outdir}/{prefix}/heatmap/{prefix}_heatmap_count.png
convert -density 400 -background white -alpha remove {outdir}/{prefix}/heatmap/{prefix}_heatmap_log.pdf {outdir}/{prefix}/heatmap/{prefix}_heatmap_log.png
'''

    print(cmd)
    return_code = subprocess.run(cmd, shell=True)
    return return_code
