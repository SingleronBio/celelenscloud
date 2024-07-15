import json
import click
import pandas as pd
import scanpy as sc
from pathlib import Path
from anndata import AnnData
from scipy.sparse import csr_matrix


#######################Convert the uploaded matrix into h5ad format#######################
# For a 10X format matrix, three files are required: 'barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz'
# (for cellranger v3 matrix) or 'barcodes.tsv', 'genes.tsv', 'matrix.mtx' (for cellranger v2 matrix).
# The difference between the two lies in 'features.tsv.gz' and 'genes.tsv'.
# The former contains three columns: gene id, gene symbol, and feature_types (for RNA data, it is 'Gene Expression'),
# while the latter contains only gene symbol and gene id.
# The three files need to be placed in a folder and compressed, either as zip, tar.gz, or tar format.
# For a TSV matrix file, the column names should be barcode names, the row names should be gene names,
# and it should be tab-separated. Note that the barcode names cannot consist entirely of numbers.
# For a h5 matrix file, any standardized matrix from CellRanger or homemade can be accepted.
# It only needs to ensure that it contains barcode, gene, and expression information.
# For a loom matrix file, it is important to note that scanpy by default
# stores the barcode names in the 'CellID' column in 'obs',
# and the gene names in the 'Gene' column in 'var'.
# For a h5ad file, matrix information will be retrieved from .X,
# and it is necessary to ensure that raw counts are stored in .X.
# Additionally, the version of anndata must be lower than 0.8.0.
# For a rds file, any valid rds will be accepted


def decompress_10x_mtx(path: str):
    """extract 10x matrix from a compressed folder"""
    p = Path(path)
    if p.is_file():
        target_dir = Path(f'{p}_files')
        target_dir.mkdir(parents=True, exist_ok=True)
        suffix = p.suffix

        supported_files = [
            'barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz',
            'barcodes.tsv', 'genes.tsv', 'matrix.mtx'
        ]

        if suffix == '.zip':
            from zipfile import ZipFile
            zip_file = ZipFile(p)
            for member in zip_file.filelist:
                if not member.is_dir():
                    file = Path(member.filename)
                    if file.name in supported_files:
                        zip_file.extract(member, path=target_dir)
                        Path(f'{target_dir}/{file}').rename(f'{target_dir}/{file.name}')

        if suffix == '.gz':
            import tarfile
            tar_file = tarfile.open(p, mode='r:gz')
            for member in tar_file.getmembers():
                if member.isfile():
                    file = Path(member.name)
                    if file.name in supported_files:
                        tar_file.extract(member, path=target_dir)
                        Path(f'{target_dir}/{file}').rename(f'{target_dir}/{file.name}')

        if suffix == '.tar':
            import tarfile
            tar_file = tarfile.open(p, mode='r')
            for member in tar_file.getmembers():
                if member.isfile():
                    file = Path(member.name)
                    if file.name in supported_files:
                        tar_file.extract(member, path=target_dir)
                        Path(f'{target_dir}/{file}').rename(f'{target_dir}/{file.name}')

        return target_dir

    elif p.is_dir():
        return p


def check_10x_version(path):
    """Check the version of the 10x matrix. If possible, repair the matrix."""
    features_path = Path(path) / 'features.tsv.gz'
    genes_path = Path(path) / 'genes.tsv'

    if features_path.exists():
        print('This file appears to be a cellRanger v3 file.')
        try:
            data = pd.read_csv(features_path, header=None, sep='\t')
            if len(data.columns) == 1:
                data[1] = data[0]
                data[2] = 'Gene Expression'
                data.to_csv(
                    features_path,
                    header=False,
                    index=False,
                    sep='\t',
                    compression='gzip'
                )
            elif len(data.columns) != 3:
                raise ValueError(
                    f'An error occurred while reading features.tsv.gz, '
                    f'it only contains {len(data.columns)} columns.'
                )
        except Exception as e:
            print(f'An error occurred while processing features.tsv.gz: {e}')

    elif genes_path.exists():
        print('This file appears to be a cellRanger v2 file.')
        try:
            data = pd.read_csv(genes_path, header=None, sep='\t')
            if len(data.columns) == 1:
                data[1] = data[0]
                data.to_csv(
                    genes_path,
                    header=False,
                    index=False,
                    sep='\t'
                )
            elif len(data.columns) != 2:
                raise ValueError(
                    f'An error occurred while reading genes.tsv, '
                    f'it only contains {len(data.columns)} columns.'
                )
        except Exception as e:
            print(f'An error occurred while processing genes.tsv: {e}')

    else:
        raise ValueError('No matrix file was found in the provided directory.')


def fix_10x_barcode_dtype(path, name):
    """Compatible with BD's 10X matrix"""
    filenames = ['barcodes.tsv.gz', 'barcodes.tsv']

    for filename in filenames:
        file_path = path / filename
        if file_path.exists():
            bc_path = str(file_path)
            break
    else:
        raise ValueError('No barcodes file found.')

    bc = pd.read_csv(bc_path, header=None)
    if bc[0].dtype == 'int64':
        print('A BD format 10x matrix found')
        bc[0] = name + '_' + bc[0].astype('str')

        if bc_path.endswith('.gz'):
            bc.to_csv(bc_path, header=False, index=False, compression='gzip')
        else:
            bc.to_csv(bc_path, header=False, index=False)

    else:
        pass


def process_and_write_data(
    adata,
    output_path,
    sample_id
):
    """process anndata object and write"""
    if not isinstance(adata.X, csr_matrix):
        adata.X = csr_matrix(adata.X, shape=adata.shape)

    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs['Sample ID'] = sample_id
    adata.write(output_path, compression='lzf')

    # save data basic information
    info_records = {
        'n_obs': adata.n_obs,
        'n_vars': adata.n_vars
    }
    with open('info_records.json', 'w') as f:
        json.dump(info_records, f)


@click.group()
def cli():
    pass


@cli.command("trans_h5")
@click.option("-i", "--input", help="Input h5 file")
@click.option("-n", "--name", help="Specify the name of the matrix")
@click.option("-o", "--output", help="Output h5ad file")
def trans_h5(input, name, output):
    adata = sc.read_10x_h5(input)
    process_and_write_data(adata, output, name)


@cli.command("trans_loom")
@click.option("-i", "--input", help="Input loom file")
@click.option("-n", "--name", help="Specify the name of the matrix")
@click.option("-o", "--output", help="Output h5ad file")
def trans_loom(input, name, output):
    adata = sc.read_loom(input)
    process_and_write_data(adata, output, name)


@cli.command("trans_tsv")
@click.option("-i", "--input", help="Input tsv file")
@click.option("-n", "--name", help="Specify the name of the matrix")
@click.option("-o", "--output", help="Output h5ad file")
def trans_tsv(input, name, output):
    adata = sc.read_csv(input, delimiter='\t').T
    process_and_write_data(adata, output, name)


@cli.command("trans_10x_mtx")
@click.option("-i", "--input", help="Input 10x_mtx file")
@click.option("-n", "--name", help="Specify the name of the matrix")
@click.option("-o", "--output", help="Output h5ad file")
def trans_10x_mtx(input, name, output):
    p = decompress_10x_mtx(input)
    _ = check_10x_version(p)
    fix_10x_barcode_dtype(p, name)

    adata = sc.read_10x_mtx(p)
    process_and_write_data(adata, output, name)


@cli.command("trans_h5ad")
@click.option("-i", "--input", help="Input h5ad file")
@click.option("-n", "--name", help="Specify the name of the matrix")
@click.option("-o", "--output", help="Output h5ad file")
def trans_h5ad(input, name, output):
    adata = sc.read(input)
    if not isinstance(adata.X, csr_matrix):
        adata.X = csr_matrix(adata.X, shape=adata.shape)
    adata = AnnData(
        X=adata.X,
        obs=pd.DataFrame(index=adata.obs_names),
        var=pd.DataFrame(index=adata.var_names),
        dtype='float32'
    )
    process_and_write_data(adata, output, name)


@cli.command("trans_rds")
@click.option("-i", "--input", help="Input rds file")
@click.option("-n", "--name", help="Specify the name of the matrix")
@click.option("-o", "--output", help="Output h5ad file")
def trans_rds(input, name, output):
    import rpy2.robjects as ro

    r_load = ro.r['readRDS']
    seurat_object = r_load(input)

    data_matrix_r = ro.r('GetAssayData')(
        seurat_object,
        slot="counts"
    )
    x = ro.r('as')(data_matrix_r, 'dgCMatrix')

    values = ro.r('as.vector')(x.do_slot('x'))
    indices = ro.r('as.vector')(x.do_slot('i'))
    indptr = ro.r('as.vector')(x.do_slot('p'))
    shape = tuple(x.slots['Dim'])[::-1]

    sparse_matrix = csr_matrix(
        (values, indices, indptr),
        shape=shape
    ).astype('float32')
    barcodes = ro.r('colnames')(data_matrix_r)
    genes = ro.r('rownames')(data_matrix_r)

    adata = AnnData(X=sparse_matrix)
    adata.obs_names = barcodes
    adata.var_names = genes
    process_and_write_data(adata, output, name)


if __name__ == "__main__":
    cli()
