#!/usr/bin/env python3

# @author: Peter van Rijn
# @date: 2026/07/08

# Decompose shWGS insert size metrics into healthy/tumor fractions using 
# non-negative matrix factorization (NMF), as described by Renaud et al. (2022,
# https://doi.org/10.7554/elife.71569), and plot the original insert sizes
# over the tumor/healthy signals provided by the aforementioned paper
# (resources/Renaud_et_al_Fig2B_WGS_Sigs_data.tsv) and include the estimated
# tumor fraction resulting from NMF. 

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import scipy as sp
import sklearn as skl
import seaborn as sn

from collections import defaultdict
from typing import List


def nmf(dfisize: pl.DataFrame, dfrenaud: pl.DataFrame) -> pl.DataFrame:
    '''Perform non-negative matrix factorization (NMF) '''

    sig1 = dfrenaud.filter( pl.col('signature_id').eq('Signature1') ).select('count').to_series().to_numpy()
    sig2 = dfrenaud.filter( pl.col('signature_id').eq('Signature2') ).select('count').to_series().to_numpy()

    # H matrix is the precomputed healthy/tumor signals from Renaud et al. (2022),
    # but we need to convert it to fraction per insert size
    sigs = np.array([sig1, sig2])
    row_sums = sigs.sum(axis=1)
    fixed_H = sigs / row_sums[:, np.newaxis]  # Normalize
    
    # Do the same for the actual insert size data
    dfisize = dfisize.with_columns(
        (pl.col('count') / pl.col('count').sum()).alias('fraction'),
    )
    normdata = dfisize.select('fraction').to_numpy().T
    normdata = normdata / normdata.sum(axis=1)

    # Calculate the W matrix
    # Parameters were taken from the M & M section in Renaud et al. (2022)
    # Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC9363120/#s4
    W, H, n_iter = skl.decomposition.non_negative_factorization(
        # dfisize.to_pandas().T,
        normdata,
        n_components=2,
        init='custom',
        solver='mu',
        beta_loss='kullback-leibler',
        update_H=False,
        H=fixed_H,
        max_iter=200,
        random_state=0
    )

    return pl.DataFrame(W, schema=['Signature 1', 'Signature 2'])


def find_isize_peaks(df: pl.DataFrame, isize_column: str='insert_size',
                     min_isize: int=50, max_isize: int=200, step: int=150):
    '''Find insert size peaks'''
    
    # Tumor signature peaks sit between 50 and 250, remove everything else
    df = df.filter(
        pl.col(isize_column).ge(min_isize) &
        pl.col(isize_column).le(max_isize)
    )

    max_isize = df.select(isize_column).max().item()
    peaks = set()

    # For some reason, not all relevant peaks are found when applying find_peaks()
    # to the full dataset. Splitting it up in smaller chunks seems to do the trick
    for center in range(min_isize, max_isize, step):

        left  = max(0, center-step)
        right = min(center+step, max_isize)

        dftmp = df.filter(pl.col(isize_column).ge(left) &
                          pl.col(isize_column).le(right))
        
        if dftmp.is_empty():
            continue

        # SciPy is not compatible with polars, convert to numpy
        data = dftmp.select('count').to_series().to_numpy()
        tmppeaks, _ = sp.signal.find_peaks(data, distance=10, height=0.001)
        
        peaks.update([peak + center for peak in tmppeaks if peak + center <= max_isize])
        # peaks.update([peak + step for peak in tmppeaks if peak + step <= 200])
    
    return sorted(peaks)


def plot(dfisz: pl.DataFrame, dfnmf: pl.DataFrame, dfrnd: pl.DataFrame,
         sample: str, figure_path: str, peaks: List=[]) -> None:
    '''Plotterdeplotterdeplot!'''
    
    fig, ax = plt.subplots(figsize=(10, 5))

    # sample insert sizes
    sn.lineplot(dfisz, x='Insert Size', y='% of Reads', linewidth=2.0,
                color='blue', hue='cfDNA signature', zorder=3, ax=ax)
    
    # NMF-calculated tumor signal fraction 
    ts = int(dfnmf.select('Signature 2').item() * 100)
    text = f'NMF Tumor Signal: {ts}%'
    ax.text(.75, .9, text, transform=ax.transAxes, fontweight="bold")

    # Plot Renaud et al. signals in the background
    sn.lineplot(data=dfrnd, x='Insert Size', y='% of Reads', estimator=None,
                hue='cfDNA signature', palette=['orange', 'brown'], linewidth=1.0,
                zorder=5, linestyle='--', ax=ax)

    # Reduce the frequency of the x axis ticks
    max_isize = dfisz.select('Insert Size').max().item()
    ax.set_xticks(sorted([n for n in range(0, max_isize, 100)] + list(peaks)))
    ax.tick_params(axis='x', labelbottom=True)
    ax.tick_params(axis='y', labelleft=True)

    # Rotate labels and set x-axis offset of peaks
    for l in ax.get_xticklabels():
        tick_txt = l.get_text()
        if tick_txt in set(map(str, peaks)):
            l.set_rotation(45)
            l.set_color('brown')
            x, y = l.get_position()
            l.set_position((x, y-0.04))

    # 
    for peak in peaks:
        ax.axvline(peak, linestyle=':', color='brown', linewidth=1.0, alpha=0.5, zorder=2)

    handles, labels = ax.get_legend_handles_labels()
    ax.get_legend().remove()
    ax.legend(handles=handles, labels=labels)
    ax.set_title(f'Sample {sample}')
    ax.set_xlabel("Insert Size")
    fig.tight_layout()
    fig.savefig(figure_path)


def smooth(df: pl.DataFrame, smoothing_factor: int=1) -> pl.DataFrame:
    '''Smooth insert size data by a factor'''
    d = defaultdict(int)

    for row in df.iter_rows(named=True):
        new_isize = (row['insert_size'] // smoothing_factor) * smoothing_factor
        d[new_isize] = d.get(new_isize, 0) + row.get('count', 0)

    return pl.from_dict({'insert_size': d.keys(), 'count': d.values()})


def add_fraction(df: pl.DataFrame, count_column='count', output_column='fraction') -> pl.DataFrame:
    '''Add a column with fractional values of the insert size counts'''
    return df.with_columns(
        ( pl.col(count_column).truediv(pl.col(count_column).sum()) ).alias(output_column)
    )


def filter_isize(df, isize_col: str='insert_size', min_isize: int=31,
                 max_isize: int=700) -> pl.DataFrame:
    '''Filter for min/max insert sizes'''
    return df.filter(
        pl.col(isize_col).ge(min_isize) &
        pl.col(isize_col).le(max_isize)
    )


def load_tumor_healthy_signatures(path: str) -> pl.DataFrame:
    '''Load healthy/tumor signals from Renaud et al. (2022)'''

    return pl.read_csv(
        path,
        separator='\t',
        new_columns=['signature_id', 'frag_len', 'count', 'insert_size']
    )
    

def load_picard_stats(path: str) -> pl.DataFrame:
    '''Load picard insert size metrics into polars dataframe, add insert size fractions'''
    
    return pl.read_csv(
        path,
        separator='\t',
        skip_rows=11,
        new_columns=['insert_size', 'count']
    )


def load_samtools_stats(path: str) -> pl.DataFrame:
    '''Load samtools insert size metrics into polars dataframe, add insert size fractions'''
    
    return pl.read_csv(
        path,
        separator='\t',
        columns=[0, 1],
        new_columns=['insert_size', 'count']
    )


def get_sample_name(path: str) -> str:
    '''Assume sample name is first field after splitting on dots and underscores'''
    
    return os.path.basename(path).split('.')[0].split('_')[0]


def parse_cmd() -> argparse.Namespace:
    parser = argparse.ArgumentParser('Non-Negative Matrix Factorization')
    # parser.add_argument('--insert-size-metrics', '-i', type=str, required=True,
    #                     help='Insert size metrics file generated by `picardtools CollectInsertSizeMetrics`')
    parser.add_argument('--samtools-isize-metrics', '-s', type=str, required=True,
                        help='''Tab separated insert size metrics file generated by
                               `samtools stats` followed by `grep '^IS' | cut -f 2-3`''')
    parser.add_argument('--signature-file', '-f', type=str, required=True,
                        help='File containing healthy/tumor insert size signals (Renaud et al. 2022)')
    parser.add_argument('--output-figure', '-o', type=str, required=True,
                        help='Output figure path (PNG, PDF)')
    parser.add_argument('--output-matrix', type=str, required=False,
                        help='Output matrix path (optional)')
    parser.add_argument('--min-insert-size', type=int, default=31,
                        help='Minimum insert size')
    parser.add_argument('--max-insert-size', type=int, default=701,
                        help='Maximum insert size')
    parser.add_argument('--smoothing-factor', type=int, default=1,
                        help='Smooth/Compress insert size data')
    parser.add_argument('--sample-name', type=str, required=False,
                        help='''Sample name to use when plotting. If not set, the 
                                script will try to infer the sample name from the
                                samtools insert size metrics file.''')
    
    return parser.parse_args()
    

def main() -> None:
    
    args = parse_cmd()

    # Load healthy/tumor cfDNA signatures from Renaud et al. (2022)
    dfrenaud = load_tumor_healthy_signatures(args.signature_file)
    dfrenaud = filter_isize(
        dfrenaud,
        min_isize=args.min_insert_size,
        max_isize=args.max_insert_size
    )

    # Load data sample data, smooth if needed and calculate insert size fractions
    dfisize = load_samtools_stats(args.samtools_isize_metrics)
    dfisize = filter_isize(
        dfisize,
        min_isize=args.min_insert_size,
        max_isize=args.max_insert_size
    )
    dfisize = smooth(dfisize, smoothing_factor=args.smoothing_factor)
    dfisize = add_fraction(dfisize)

    # both dataframe dimensions (our's and Renaud's insert sizes) should match for NMF
    if len(dfisize) < len(dfrenaud):
        dfrenaud = dfrenaud.filter(
            pl.col('insert_size').is_in(dfisize.select('insert_size').to_series().implode())
        )
    elif len(dfisize) > len(dfrenaud):
        dfisize = dfisize.filter(
            pl.col('insert_size').is_in(dfrenaud.select('insert_size').to_series().implode())
        )

    # Split insert sizes in healthy/tumor fractions with NMF
    dfnmf = nmf(dfisize, dfrenaud)
    
    # Retrieve tumor signature and perform peak calculation
    dftumor = dfrenaud.filter(
        pl.col('signature_id').eq('Signature2')
    )
    tumor_peaks = find_isize_peaks(dftumor)

    # Prettify column names for plotting
    dfisize = (
        dfisize.with_columns(
            pl.lit('Sample signature').alias('cfDNA signature')
        )
        .rename({
                'fraction':    '% of Reads',
                'insert_size': 'Insert Size',
        })
    )

    dfrenaud = (
        dfrenaud.with_columns(  # make plotting friendly names
            pl.col('signature_id').str.replace_many(
                ('Signature1', 'Signature2'),
                ('Healthy signature', 'Tumor signature'),
            )
        )
        .rename({
            'count':        '% of Reads',
            'insert_size':  'Insert Size',
            'signature_id': 'cfDNA signature'
        })
    )

    sample = args.sample_name if args.sample_name else get_sample_name(args.samtools_isize_metrics)

    # PLOTTTTTT
    plot(dfisize, dfnmf, dfrenaud, sample, args.output_figure, tumor_peaks)


if __name__ == '__main__':
    main()