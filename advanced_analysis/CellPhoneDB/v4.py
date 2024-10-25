import os
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

def run_v4(
    cpdb_file_path: str,
    meta_file_path: str,
    counts_file_path: str,
    output_path: str,
    output_suffix: str,
    counts_data: str="ensembl",
    microenvs_file_path: str=None,
    method: str='statistical',
    degs_file_path: str=None,
):
    if method == 'statistical':
        result = cpdb_statistical_analysis_method.call(
            cpdb_file_path = cpdb_file_path,
            meta_file_path = meta_file_path,
            counts_file_path = counts_file_path,
            counts_data = "ensembl",
            output_path = output_path,
            output_suffix = output_suffix,
            microenvs_file_path = microenvs_file_path
        )
    else:
        result = cpdb_degs_analysis_method.call(
            cpdb_file_path = cpdb_file_path,
            meta_file_path = meta_file_path,
            counts_file_path = counts_file_path,
            counts_data = "ensembl",
            output_path = output_path,
            output_suffix = output_suffix,
            microenvs_file_path = microenvs_file_path,
            degs_file_path = degs_file_path
        )
    return result