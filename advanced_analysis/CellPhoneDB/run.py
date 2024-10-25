#!/usr/bin/env python
# -*- coding: utf-8 -*-
import click
from cpdb import run_pipeline

def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo('Singleron cellphonedb version 1.0')
    ctx.exit()

@click.command()
@click.option('-v', '--version', is_flag=True, callback=print_version, expose_value=False, is_eager=True)
@click.option('-r', '--rds', help='rds file', required=True, type=str)
@click.option('--prefix', default="cellphonedb", help='prefix', type=str)
@click.option('--outdir', default=".", help='output dir', type=str)
@click.option('--species', default="homo_sapiens", help='the species, include homo_sapiens, mus_musculus, rattus_norvegicus, macaca_mulatta, sus_scrofa', type=str)
@click.option('--clusterorder', default="F", help='clusterorder', type=str)
@click.option('--subcluster', default="all", help='subcluster list ,split by ,', type=str)
@click.option('--sample', default="all", help='"sample list ,split by ,', type=str)
@click.option('--group', default="all", help='group list ,split by ,', type=str)
@click.option('--platform', default="RNA", help='the sequencing platform,RNA or Spatial', type=str)
@click.option('--ncell', default=500, help='the cell number extracted from each celltype', type=int)

@click.option('--cpdb_version', default="v4", help='the analytical version of cellphonedb', type=str)
@click.option('--threads', default=3, help='Max of threads to process the data', type=int)
@click.option('--threshold', default=0.1, help='percent of cells expressing a gene', type=float)
@click.option('--iterations', default=1000, help='Number of pvalues analysis iterations', type=int)
@click.option('--method', default="statistical", help='cellphonedb method: "statistical", "degs"', type=str)

@click.option('--degs', default="F", help='differential expression file', type=str)
@click.option('--diffcluster', default="F", help='compare cluster name: T/NK:B,B/MPs; split by ,', type=str)
@click.option('--diffgroup', default="F", help='compare group name: BCells=G1:G2/G3,TCells=G1:G2/G3; split by ,', type=str)
@click.option('--logfc', default=0.2, help='logfc', type=float)
@click.option('--p_adj_cutoff', default=0.01, help='differential expression file', type=float)
@click.option('--p_cutoff', default=1, help='differential expression file', type=float)

@click.option('--microenvs_file_path', default="F", help='micro-environment file, Its content is used to limit cluster interactions', type=str)
@click.option('--newcol', default="auto", help='celltype color, include auto, color_v1 and color_v2, color_v1 will use cluster_colors column of rds meta.data; color_v2 will use newcol', type=str)
@click.option('--display_numbers', default="F", help='display numbers of heatmap', type=str)
@click.option('--flip', default="F", help='whether flip the x and y', type=str)
def main(*args, **kwargs):
    run_pipeline(kwargs)

if __name__ == '__main__':
    main()