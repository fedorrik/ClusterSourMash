#!/usr/bin/env python3
import argparse
import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


def read_matrix(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep='\t')
    except Exception as e:
        raise ValueError(f"failed to read {path}: {e}")
    if 'sample' not in df.columns:
        raise ValueError(f"{path}: expected first column named 'sample'")
    df = df.set_index('sample')
    if df.shape[0] != df.shape[1]:
        raise ValueError(f"{path}: matrix must be square, got {df.shape}")
    if list(df.index) != list(df.columns):
        raise ValueError(f"{path}: row/column labels must match and be in same order")
    try:
        df = df.apply(pd.to_numeric)
    except Exception as e:
        raise ValueError(f"{path}: matrix contains non-numeric values: {e}")
    return df


def choose_figsize(n: int) -> float:
    return max(6, min(18, 0.35 * n))


def save_matrix(df: pd.DataFrame, path: str) -> None:
    out = df.copy()
    out.index.name = 'sample'
    out.to_csv(path, sep='\t')


def plot_clustermap(df: pd.DataFrame, linkage_matrix, out_png: str, title: str, cmap: str, dpi: int, mask_diagonal: bool) -> None:
    plot_df = df.copy()
    if mask_diagonal:
        np.fill_diagonal(plot_df.values, np.nan)

    g = sns.clustermap(
        plot_df,
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        cmap=cmap,
        figsize=(choose_figsize(df.shape[0]), choose_figsize(df.shape[0])),
        cbar_kws={'label': title},
    )
    g.ax_heatmap.set_title(title)
    plt.savefig(out_png, dpi=dpi, bbox_inches='tight')
    plt.close(g.fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Summarize support matrices: mean, SD, and coefficient of variation. Row/column order is taken from the reference matrix.'
    )
    parser.add_argument('--bootstrap-dir', required=True, help='Directory with support matrices (*.tsv).')
    parser.add_argument('--reference-matrix', required=True, help='Main/original matrix used for clustering order.')
    parser.add_argument('--glob', default='*.tsv', help='Glob for support matrices [default: %(default)s]')
    parser.add_argument('--out-dir', default='support_matrix_comparison', help='Output directory [default: %(default)s]')
    parser.add_argument('--prefix', default='support_matrix_variation', help='Output file prefix inside --out-dir [default: %(default)s]')
    parser.add_argument('--dpi', type=int, default=300, help='Image DPI [default: %(default)s]')
    parser.add_argument('--std-cmap', default='viridis', help='Colormap for SD clustermap [default: %(default)s]')
    parser.add_argument('--cv-cmap', default='magma', help='Colormap for CV clustermap [default: %(default)s]')
    parser.add_argument('--mask-diagonal', action='store_true', help='Mask diagonal in clustermaps.')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    ref = read_matrix(args.reference_matrix)
    pattern = os.path.join(args.bootstrap_dir, args.glob)
    paths = sorted(glob.glob(pattern))
    if not paths:
        sys.exit(f'No matrix files found for pattern: {pattern}')

    print(f'Found {len(paths)} support matrices')
    mats = []
    for path in paths:
        mat = read_matrix(path)
        if list(mat.index) != list(ref.index):
            sys.exit(f"Label mismatch between reference matrix and '{path}'")
        mats.append(mat.values)

    arr = np.stack(mats, axis=0)
    mean = arr.mean(axis=0)
    std = arr.std(axis=0)
    eps = 1e-12
    cv = std / (mean + eps)

    mean_df = pd.DataFrame(mean, index=ref.index, columns=ref.columns)
    std_df = pd.DataFrame(std, index=ref.index, columns=ref.columns)
    cv_df = pd.DataFrame(cv, index=ref.index, columns=ref.columns)

    z = linkage(squareform(ref.values), method='average')

    pfx = os.path.join(args.out_dir, args.prefix)
    save_matrix(mean_df, pfx + '.mean.tsv')
    save_matrix(std_df, pfx + '.std.tsv')
    save_matrix(cv_df, pfx + '.cv.tsv')

    with open(pfx + '.summary.txt', 'w') as out:
        out.write(f'n_support_matrices\t{len(paths)}\n')
        out.write(f'n_samples\t{ref.shape[0]}\n')
        out.write(f'mean_sd\t{np.nanmean(std):.8f}\n')
        out.write(f'max_sd\t{np.nanmax(std):.8f}\n')
        out.write(f'mean_cv\t{np.nanmean(cv):.8f}\n')
        out.write(f'max_cv\t{np.nanmax(cv):.8f}\n')

    plot_clustermap(std_df, z, pfx + '.std.clustermap.png', 'Distance SD', args.std_cmap, args.dpi, args.mask_diagonal)
    plot_clustermap(cv_df, z, pfx + '.cv.clustermap.png', 'Coefficient of Variation', args.cv_cmap, args.dpi, args.mask_diagonal)

    print('Saved outputs to', os.path.abspath(args.out_dir))
