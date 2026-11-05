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

from collections import defaultdict
from typing import List, Union

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import polars.selectors as cs
import scipy as sp
import sklearn as skl
import seaborn as sn


def nmf(dfisize: pl.DataFrame, dfrenaud: pl.DataFrame) -> pl.DataFrame:
    '''Perform non-negative matrix factorization (NMF) '''

    isizes = dfisize.get_column('insert_size')
    dfrenaud = dfrenaud.filter(
        pl.col('insert_size').ge(isizes.min()) &
        pl.col('insert_size').le(isizes.max())
    )

    # Extract Healthy/Tumor signals
    sig1 = dfrenaud.filter( pl.col('signature_id').eq('Signature1') ).get_column('count').to_numpy()
    sig2 = dfrenaud.filter( pl.col('signature_id').eq('Signature2') ).get_column('count').to_numpy()

    # H matrix is the precomputed healthy/tumor signals from Renaud et al. (2022),
    # but we need to convert it to fraction per insert size
    sigs = np.array([sig1, sig2])
    row_sums = sigs.sum(axis=1)
    fixed_H = sigs / row_sums[:, np.newaxis]  # Normalize

    # Pivot to wide format: insert_size as columns and fractions as values, while
    # keeping sample_name as index column for rows
    dfisize = (
        dfisize.pivot(on='insert_size', index='sample_name', values='fraction')
            .fill_null(0)
    )
    normdata = dfisize.select(cs.numeric()).to_numpy()

    # Calculate the W matrix
    # Parameters were taken from the M & M section in Renaud et al. (2022)
    # Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC9363120/#s4
    W, H, n_iter = skl.decomposition.non_negative_factorization(
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

    return pl.DataFrame(W, schema=['Signature 1', 'Signature 2']).with_columns(
        dfisize.get_column('sample_name')
    )


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

    return sorted(peaks)


def plot(dfisz: pl.DataFrame, dfnmf: pl.DataFrame, dfrnd: pl.DataFrame,
         figure_path: str, peaks: List=[], sample: str='') -> None:
    '''Plotterdeplotterdeplot!'''

    fig, ax = plt.subplots(figsize=(10, 5))

    isizes = dfisz.get_column('Insert Size')
    dfrnd = dfrnd.filter(
        pl.col('Insert Size').ge(isizes.min()) &
        pl.col('Insert Size').le(isizes.max())
    )

    # Plot sample insert sizes. Seaborn does not support polars dataframes;
    # these need to be converted to pandas dataframes
    sn.lineplot(dfisz.to_pandas(), x='Insert Size', y='% of Reads', linewidth=2.0,
                color='blue', hue='cfDNA signature', zorder=3, ax=ax)

    # NMF-calculated tumor signal fraction
    ts = int(dfnmf.select('Signature 2').item() * 100)
    text = f'NMF Tumor Signal: {ts}%'
    ax.text(.75, .9, text, transform=ax.transAxes, fontweight="bold")

    # Plot Renaud et al. signals in the background. Seaborn does not support
    # polars dataframes; these need to be converted to pandas dataframes
    sn.lineplot(data=dfrnd.to_pandas(), x='Insert Size', y='% of Reads', estimator=None,
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

    # Plot a vertical line for each determined peak
    for peak in peaks:
        ax.axvline(peak, linestyle=':', color='brown', linewidth=1.0, alpha=0.5, zorder=2)

    handles, labels = ax.get_legend_handles_labels()
    ax.get_legend().remove()
    ax.legend(handles=handles, labels=labels)
    if sample:
        ax.set_title(f'Sample {sample}')
    ax.set_xlabel("Insert Size")
    fig.tight_layout()
    fig.savefig(figure_path, dpi=600)


def smooth(df: pl.DataFrame, smoothing_factor: int=1) -> pl.DataFrame:
    '''Smooth insert size data by a factor'''
    d = defaultdict(int)

    for row in df.iter_rows(named=True):
        new_isize = (row['insert_size'] // smoothing_factor) * smoothing_factor
        d[new_isize] = d.get(new_isize, 0) + row.get('count', 0)

    newdf = pl.from_dict({'insert_size': d.keys(), 'count': d.values()})

    sample = df.get_column('sample_name').unique()[0]
    newdf = newdf.with_columns(sample_name=pl.lit(sample))

    return newdf


def add_fraction(df: pl.DataFrame, count_column='count', output_column='fraction') -> pl.DataFrame:
    '''Add a column with fractional values of the insert size counts'''
    return df.with_columns(
        pl.col(count_column).truediv(pl.col(count_column).sum()).alias(output_column)
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


def load_samtools_stats(path: str, sample: str) -> pl.DataFrame:
    '''Load samtools insert size metrics into polars dataframe, add insert size fractions'''

    return (
        pl.read_csv(
            path,
            separator='\t',
            columns=[0, 1],
            new_columns=['insert_size', 'count']
        )
        .with_columns(
            pl.lit(sample).alias('sample_name')
        )
    )


def get_sample_name(path: str) -> str:
    '''Assume sample name is first field after splitting on dots and underscores'''

    return os.path.basename(path).split('.')[0].split('_')[0]


def parse_cmd() -> argparse.Namespace:

    parser = argparse.ArgumentParser('Non-Negative Matrix Factorization')

    parser.add_argument('--samtools-isize-metrics', '-s', type=str, nargs='+', required=True,
                        help='''Tab separated insert size metrics file(s) generated by
                               `samtools stats` followed by `grep '^IS' | cut -f 2-3`''')

    parser.add_argument('--signature-file', '-f', type=str, required=True,
                        help='File containing healthy/tumor insert size signals (Renaud et al. 2022)')

    parser.add_argument('--output-directory', '-o', type=str, required=True,
                        help='Directory where output figures will be written to (PNG or PDF)')

    parser.add_argument('--output-matrix', action='store_true',
                        help='Write NMF matrix to TSV file (optional)')

    parser.add_argument('--min-insert-size', type=int, default=31,
                        help='Minimum insert size')

    parser.add_argument('--max-insert-size', type=int, default=701,
                        help='Maximum insert size')

    parser.add_argument('--smoothing-factor', type=int, default=1,
                        help='Smooth/Compress insert size data')

    return parser.parse_args()


def main() -> None:

    args = parse_cmd()

    # Load healthy/tumor cfDNA signatures from Renaud et al. (2022)
    dfrenaud = load_tumor_healthy_signatures(args.signature_file)

    samples = [get_sample_name(f) for f in args.samtools_isize_metrics]

    # Load data sample data, filter & smooth if needed and calculate insert size fractions
    dfisize = [load_samtools_stats(f, sample)
               for f, sample in zip(args.samtools_isize_metrics, samples)]

    dfisize = [filter_isize(df, min_isize=args.min_insert_size, max_isize=args.max_insert_size)
               for df in dfisize]

    dfisize = [smooth(df, smoothing_factor=args.smoothing_factor)
               for df in dfisize]

    # Add fractions and concatenate to a single dataframe
    dfisize = pl.concat([add_fraction(df) for df in dfisize], how='vertical')

    # Split insert sizes in healthy/tumor fractions with NMF
    dfnmf = nmf(dfisize, dfrenaud)

    # Retrieve tumor signature and perform peak calculation
    dftumor = dfrenaud.filter( pl.col('signature_id').eq('Signature2') )
    tumor_peaks = find_isize_peaks(dftumor)

    # Prettify column names for plotting
    dfisize = dfisize.with_columns(
        pl.lit('Sample signature').alias('cfDNA signature')
    ).rename({
        'fraction':    '% of Reads',
        'insert_size': 'Insert Size',
    })

    dfrenaud = dfrenaud.with_columns(
        pl.col('signature_id').str.replace_many(
            ('Signature1',        'Signature2'),
            ('Healthy signature', 'Tumor signature')
        )
    ).rename({
        'count':        '% of Reads',
        'insert_size':  'Insert Size',
        'signature_id': 'cfDNA signature'
    })

    # Don't create output directory until we know for sure NMF succeeded
    os.makedirs(args.output_directory, exist_ok=True)

    for sample, dfn in dfnmf.partition_by('sample_name', as_dict=True).items():

        if isinstance(sample, tuple):
            sample = sample[0]

        # PLOTTTTTT
        dfisz = dfisize.filter(pl.col('sample_name').eq(sample))
        output_figure_path = f'{args.output_directory}/{sample}.isize-NMF.png'

        plot(dfisz, dfn, dfrenaud, output_figure_path, sample=sample,
             peaks=tumor_peaks)

        if args.output_matrix:
            output_matrix_path = f'{args.output_directory}/{sample}.isize-NMF.tsv'
            dfn.write_csv(output_matrix_path, separator='\t')


if __name__ == '__main__':
    main()