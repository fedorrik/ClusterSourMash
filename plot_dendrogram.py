#!/usr/bin/env python3
import argparse
import glob
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform


def read_matrix(path: str) -> pd.DataFrame:
    try:
        m = pd.read_csv(path, sep='\t').set_index('sample')
    except Exception as e:
        sys.exit(f"Failed to read matrix '{path}': {e}")
    if m.shape[0] != m.shape[1]:
        sys.exit(f"Matrix '{path}' must be square.")
    return m.apply(pd.to_numeric)


def compute_linkage(m: pd.DataFrame):
    return linkage(squareform(m.values), method='average')


def clusters_from_linkage(z: np.ndarray, labels):
    n = len(labels)
    clusters = {i: frozenset([labels[i]]) for i in range(n)}
    for i, row in enumerate(z):
        a, b = int(row[0]), int(row[1])
        clusters[n + i] = clusters[a] | clusters[b]
    return clusters


def internal_nodes_from_linkage(z: np.ndarray, labels):
    n = len(labels)
    clusters = clusters_from_linkage(z, labels)
    nodes = []
    for i in range(n - 1):
        node_id = n + i
        nodes.append({
            'node_id': node_id,
            'leaves': clusters[node_id],
            'height': float(z[i, 2]),
            'size': int(z[i, 3]),
            'row_idx': i,
        })
    return nodes


def select_nodes(nodes, n_leaves):
    return [x for x in nodes if x['size'] < n_leaves]  # exclude root


def read_support_linkages(support_dir: str, labels):
    files = sorted(glob.glob(os.path.join(support_dir, '*.tsv')))
    if not files:
        sys.exit(f"No support matrices (*.tsv) found in '{support_dir}'.")
    z_list = []
    for path in files:
        m = read_matrix(path)
        if list(m.index) != list(labels):
            sys.exit(f"Label mismatch in support matrix '{path}'.")
        z_list.append(compute_linkage(m))
    return files, z_list


def compute_support(main_z, support_zs, labels):
    main_nodes = internal_nodes_from_linkage(main_z, labels)
    selected = select_nodes(main_nodes, len(labels))
    if not support_zs:
        return {}, selected

    counts = defaultdict(int)
    for z in support_zs:
        clades = {x['leaves'] for x in internal_nodes_from_linkage(z, labels) if x['size'] < len(labels)}
        for node in selected:
            if node['leaves'] in clades:
                counts[node['leaves']] += 1

    support_pct = {
        node['leaves']: 100.0 * counts[node['leaves']] / len(support_zs)
        for node in selected
    }
    return support_pct, selected


def annotate_support(ax, dendr, support_pct, nodes, min_height=0.0, fontsize=9):
    if not support_pct:
        return
    node_by_row = {node['row_idx']: node for node in nodes}
    for row_idx, (xs, ys) in enumerate(zip(dendr['icoord'], dendr['dcoord'])):
        y = ys[1]
        if y < min_height:
            continue
        node = node_by_row.get(row_idx)
        if node is None:
            continue
        support = support_pct.get(node['leaves'])
        if support is None:
            continue
        x = 0.5 * (xs[1] + xs[2])
        ax.text(x, y, f"{int(round(support))}", ha='center', va='bottom', fontsize=fontsize)


def save_clusters_file(path, dendr, m):
    ordered_samples = [m.index[i] for i in dendr['leaves']]
    ordered_clusters = dendr['color_list']
    with open(path, 'w') as out:
        for sample, cluster in zip(ordered_samples, ordered_clusters):
            out.write(f"{sample}\t{cluster}\n")


def plot_single_dendrogram(m, out_path, threshold=None, no_axis=False, width=None, height=None):
    n = m.shape[0]
    dw = width or max(int(round(0.15 * n)), 5)
    dh = height or int(round(dw / 3))
    z = compute_linkage(m)
    thr = threshold if threshold is not None else 0.7 * z[:, 2].max()
    fig, ax = plt.subplots(figsize=(dw, dh))
    dendr = dendrogram(z, labels=m.index.tolist(), color_threshold=thr, leaf_rotation=90, ax=ax)
    if no_axis:
        ax.axis('off')
    fig.savefig(out_path, bbox_inches='tight', dpi=200)
    plt.close(fig)
    return dendr, z


def plot_clustermap(m, out_path, size=None, no_axis=False):
    n = m.shape[0]
    cs = size or max(int(round(0.3 * n)), 5)
    g = sns.clustermap(m, cmap='coolwarm', figsize=(cs, cs))
    if no_axis:
        g.ax_heatmap.axis('off')
    g.fig.savefig(out_path, bbox_inches='tight', dpi=300)
    plt.close(g.fig)


def plot_support_dendrograms(support_dir, labels, threshold=None, no_axis=False, width=None, height=None):
    files = sorted(glob.glob(os.path.join(support_dir, '*.tsv')))
    if not files:
        sys.exit(f"No support matrices (*.tsv) found in '{support_dir}'.")
    print(f"  Plotting dendrograms for {len(files)} support matrices into {support_dir}")
    for path in files:
        m = read_matrix(path)
        if list(m.index) != list(labels):
            sys.exit(f"Label mismatch in support matrix '{path}'.")
        stem = os.path.splitext(os.path.basename(path))[0]
        out_path = os.path.join(support_dir, f'{stem}.dendrogram.png')
        plot_single_dendrogram(m, out_path, threshold=threshold, no_axis=no_axis, width=width, height=height)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw dendrogram and optional clustermap from one distance matrix, with optional branch support from support matrices.')
    parser.add_argument('matrix', help="TSV matrix with first column 'sample'.")
    parser.add_argument('--dendrogram-out', default='dendrogram.png', help='Dendrogram PNG [default: %(default)s]')
    parser.add_argument('--clustermap-out', default='clustermap.png', help='Clustermap PNG [default: %(default)s]')
    parser.add_argument('--support-dir', default=None, help='Directory with support matrices (*.tsv).')
    parser.add_argument('--support-min-height', type=float, default=0.0, help='Annotate support only for branches with y >= this value [default: %(default)s]')
    parser.add_argument('--support-fontsize', type=int, default=9, help='Font size for support labels.')
    parser.add_argument('--plot-support-dendrograms', action='store_true', help='Also draw a dendrogram for each support matrix and save it into that same directory.')
    parser.add_argument('--threshold', type=float, default=None, help='Dendrogram color threshold; default=0.7*max(height).')
    parser.add_argument('--skip-clustermap', action='store_true', help='Skip clustermap drawing.')
    parser.add_argument('--no-axis', action='store_true', help='Hide axes.')
    parser.add_argument('--dendrogram-width', type=int, default=None, help='Dendrogram width in inches.')
    parser.add_argument('--dendrogram-height', type=int, default=None, help='Dendrogram height in inches.')
    parser.add_argument('--clustermap-size', type=int, default=None, help='Clustermap size in inches.')
    args = parser.parse_args()

    m = read_matrix(args.matrix)
    print('  Matrix shape:', m.shape)
    print('  Computing linkage (average)')
    z = compute_linkage(m)
    thr = args.threshold if args.threshold is not None else 0.7 * z[:, 2].max()

    support_pct = {}
    selected_nodes = []
    if args.support_dir is not None:
        print(f'  Reading support matrices from {args.support_dir}')
        support_files, support_zs = read_support_linkages(args.support_dir, m.index.tolist())
        print(f'  Loaded {len(support_files)} support matrices')
        support_pct, selected_nodes = compute_support(z, support_zs, m.index.tolist())
        out_tsv = args.dendrogram_out + '.support.tsv'
        rows = []
        for node in sorted(selected_nodes, key=lambda x: x['height'], reverse=True):
            if node['height'] < args.support_min_height:
                continue
            rows.append({
                'height': node['height'],
                'size': node['size'],
                'support_percent': round(support_pct.get(node['leaves'], np.nan), 3),
                'leaves': ','.join(sorted(node['leaves'])),
            })
        pd.DataFrame(rows).to_csv(out_tsv, sep='\t', index=False)
        print('  Saved branch support table →', out_tsv)
        if args.plot_support_dendrograms:
            plot_support_dendrograms(args.support_dir, m.index.tolist(), threshold=args.threshold, no_axis=args.no_axis, width=args.dendrogram_width, height=args.dendrogram_height)

    n = m.shape[0]
    dw = args.dendrogram_width or max(int(round(0.15 * n)), 5)
    dh = args.dendrogram_height or int(round(dw / 3))
    print('  Plotting dendrogram →', args.dendrogram_out)
    fig, ax = plt.subplots(figsize=(dw, dh))
    dendr = dendrogram(z, labels=m.index.tolist(), color_threshold=thr, leaf_rotation=90, ax=ax)
    if support_pct:
        annotate_support(ax, dendr, support_pct, selected_nodes, min_height=args.support_min_height, fontsize=args.support_fontsize)
    if args.no_axis:
        ax.axis('off')
    fig.savefig(args.dendrogram_out, bbox_inches='tight', dpi=200)
    plt.close(fig)

    save_clusters_file(args.dendrogram_out + '.clusters', dendr, m)
    print('  Saved cluster order →', args.dendrogram_out + '.clusters')

    if not args.skip_clustermap:
        print('  Plotting clustermap →', args.clustermap_out)
        plot_clustermap(m, args.clustermap_out, size=args.clustermap_size, no_axis=args.no_axis)