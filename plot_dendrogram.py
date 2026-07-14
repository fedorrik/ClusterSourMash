#!/usr/bin/env python3
import argparse
import glob
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, MultipleLocator
from matplotlib.colors import Normalize, to_hex
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform



def read_matrix(path: str) -> pd.DataFrame:
    try:
        m = pd.read_csv(path, sep='\t').set_index('sample')
    except Exception as e:
        sys.exit(f"Failed to read matrix '{path}': {e}")
    if m.shape[0] != m.shape[1]:
        sys.exit(f"Matrix '{path}' must be square.")
    return m.apply(pd.to_numeric)


def read_sample_group_colors(path: str):
    colors = {}
    try:
        with open(path) as fh:
            for line_no, line in enumerate(fh, start=1):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    sys.exit(f"Invalid color TSV line {line_no} in '{path}': expected group<TAB>R,G,B.")
                group, rgb_text = parts
                try:
                    rgb = tuple(int(x) for x in rgb_text.split(','))
                except ValueError:
                    sys.exit(f"Invalid RGB value on line {line_no} in '{path}': '{rgb_text}'.")
                if len(rgb) != 3 or any(x < 0 or x > 255 for x in rgb):
                    sys.exit(f"Invalid RGB value on line {line_no} in '{path}': '{rgb_text}'.")
                colors[group] = tuple(x / 255.0 for x in rgb)
    except OSError as e:
        sys.exit(f"Failed to read sample group colors '{path}': {e}")
    if not colors:
        sys.exit(f"No sample group colors found in '{path}'.")
    return colors


def validate_sample_group_colors(labels, group_colors, path):
    missing = sorted({label.split('_')[-1] for label in labels} - set(group_colors))
    if missing:
        sys.exit(
            f"Sample group colors file '{path}' is missing groups matching sample name suffixes: "
            + ','.join(missing)
        )


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


def quote_newick_label(label):
    label = str(label)
    if label == '':
        return "''"
    if any(c in label for c in "()[]':;, \t\n\r"):
        return "'" + label.replace("'", "''") + "'"
    return label


def format_newick_length(length):
    return f"{max(float(length), 0.0):.10g}"


def linkage_to_newick(z, labels, support_pct=None):
    support_pct = support_pct or {}
    n = len(labels)
    clusters = {i: frozenset([labels[i]]) for i in range(n)}
    children = {}
    heights = {i: 0.0 for i in range(n)}

    for i, row in enumerate(z):
        node_id = n + i
        left, right = int(row[0]), int(row[1])
        children[node_id] = (left, right)
        heights[node_id] = float(row[2])
        clusters[node_id] = clusters[left] | clusters[right]

    root = n + len(z) - 1

    def serialize(node_id, parent_height=None):
        node_height = heights[node_id]
        branch = '' if parent_height is None else ':' + format_newick_length(parent_height - node_height)
        if node_id < n:
            return quote_newick_label(labels[node_id]) + branch

        left, right = children[node_id]
        label = ''
        if node_id != root and clusters[node_id] in support_pct:
            label = str(int(round(support_pct[clusters[node_id]])))
        subtree = f"({serialize(left, node_height)},{serialize(right, node_height)}){label}"
        return subtree + branch

    return serialize(root) + ';'


def save_newick(path, z, labels, support_pct=None):
    with open(path, 'w') as out:
        out.write(linkage_to_newick(z, labels, support_pct=support_pct) + '\n')


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


def draw_sample_group_dots(ax, dendr, labels, group_colors, show_labels=True, fontsize=8):
    if not group_colors:
        return
    y_top = ax.get_ylim()[1]
    y_pad = y_top * 0.04 if y_top > 0 else 0.04
    dot_y = -0.15 * y_pad
    text_y = -0.38 * y_pad
    for leaf_pos, leaf_idx in enumerate(dendr['leaves']):
        sample = labels[leaf_idx]
        group = sample.split('_')[-1]
        x = 5 + 10 * leaf_pos
        ax.scatter(
            [x],
            [dot_y],
            s=28,
            color=[group_colors[group]],
            edgecolors='none',
            linewidths=0,
            clip_on=False,
            zorder=5,
        )
        if show_labels:
            ax.text(x, text_y, sample, ha='center', va='top', rotation=90, fontsize=fontsize, clip_on=False)
    y_bottom = text_y - 0.15 * y_pad if show_labels else dot_y - 0.35 * y_pad
    ax.set_ylim(y_bottom, y_top)
    ax.set_xticks([])


def branch_group_color_func(z, labels, group_colors, mixed_color='black'):
    clusters = clusters_from_linkage(z, labels)
    color_by_node = {}
    for node_id, leaves in clusters.items():
        groups = {label.split('_')[-1] for label in leaves}
        if len(groups) == 1:
            color_by_node[node_id] = to_hex(group_colors[next(iter(groups))])
        else:
            color_by_node[node_id] = mixed_color
    return lambda node_id: color_by_node.get(node_id, mixed_color)


def color_terminal_leaf_branches(ax, dendr, labels, group_colors):
    if not group_colors:
        return

    x_to_color = {}
    for leaf_pos, leaf_idx in enumerate(dendr['leaves']):
        sample = labels[leaf_idx]
        group = sample.split('_')[-1]
        x_to_color[5 + 10 * leaf_pos] = to_hex(group_colors[group])

    segments = []
    colors = []
    for xs, ys in zip(dendr['icoord'], dendr['dcoord']):
        for x0, y0, x1, y1 in ((xs[0], ys[0], xs[1], ys[1]), (xs[3], ys[3], xs[2], ys[2])):
            if y0 == 0 and x0 == x1 and x0 in x_to_color:
                segments.append([(x0, y0), (x1, y1)])
                colors.append(x_to_color[x0])

    if not segments:
        return

    linewidth = None
    if ax.collections:
        linewidths = ax.collections[0].get_linewidths()
        if len(linewidths):
            linewidth = linewidths[0]
    ax.add_collection(LineCollection(segments, colors=colors, linewidths=linewidth, zorder=4))


def style_normal_axis(ax, y_axis_step=None):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_position(('outward', 3))
    if y_axis_step is None:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=12))
    else:
        ax.yaxis.set_major_locator(MultipleLocator(y_axis_step))
    ax.tick_params(axis='y', which='major', length=4)


def add_sample_group_legend(fig, group_colors, show_labels=True, branch_colors=False):
    # If sample labels are drawn under the dots, the legend is redundant.
    # Keep the legend only for label-free plots, e.g. with --no-axis.
    if not group_colors or show_labels:
        return
    handles = []
    for group, color in sorted(group_colors.items()):
        if branch_colors:
            handles.append(Line2D([0], [0], color=color, linewidth=2, label=group))
        else:
            handles.append(
                Line2D(
                    [0],
                    [0],
                    marker='o',
                    linestyle='None',
                    markerfacecolor=color,
                    markeredgecolor='none',
                    markersize=6,
                    label=group,
                )
            )
    fig.legend(
        handles=handles,
        labels=[h.get_label() for h in handles],
        loc='lower center',
        bbox_to_anchor=(0.5, 0.08),
        frameon=False,
        borderaxespad=0.0,
        handletextpad=0.15,
        labelspacing=0.2,
        ncol=len(handles),
        columnspacing=1.2,
    )


def draw_dendrogram(
    ax,
    z,
    labels,
    threshold,
    sample_group_colors=None,
    no_axis=False,
    no_labels=False,
    color_branches_by_group=False,
    y_axis_step=None,
):
    link_color_func = None
    if color_branches_by_group:
        link_color_func = branch_group_color_func(z, labels, sample_group_colors)
    dendr = dendrogram(
        z,
        labels=None if sample_group_colors and not color_branches_by_group else labels,
        color_threshold=None if color_branches_by_group else threshold,
        leaf_rotation=90,
        no_labels=no_labels or (sample_group_colors is not None and not color_branches_by_group),
        link_color_func=link_color_func,
        ax=ax,
    )
    if color_branches_by_group:
        color_terminal_leaf_branches(ax, dendr, labels, sample_group_colors)
    if not color_branches_by_group:
        draw_sample_group_dots(ax, dendr, labels, sample_group_colors, show_labels=not no_labels)
    if not no_axis:
        style_normal_axis(ax, y_axis_step=y_axis_step)
    return dendr


def save_clusters_file(path, dendr, m, Z, threshold):
    clusters = fcluster(Z, t=threshold, criterion='distance')
    ordered_samples = [m.index[i] for i in dendr['leaves']]
    ordered_clusters = [clusters[i] for i in dendr['leaves']]
    with open(path, 'w') as out:
        for sample, cluster in zip(ordered_samples, ordered_clusters):
            out.write(f"{sample}\t{cluster}\n")


def plot_single_dendrogram(
    m,
    out_path,
    threshold=None,
    no_axis=False,
    no_labels=False,
    width=None,
    height=None,
    height_per_y_unit=None,
    y_axis_step=None,
    sample_group_colors=None,
    color_branches_by_group=False,
):
    n = m.shape[0]
    dw = width or max(int(round(0.15 * n)), 5)
    z = compute_linkage(m)
    dh = height or (height_per_y_unit * z[:, 2].max() if height_per_y_unit is not None else int(round(dw / 3)))
    thr = threshold if threshold is not None else 0.7 * z[:, 2].max()
    fig, ax = plt.subplots(figsize=(dw, dh))
    dendr = draw_dendrogram(
        ax,
        z,
        m.index.tolist(),
        thr,
        sample_group_colors=sample_group_colors,
        no_axis=no_axis,
        no_labels=no_labels,
        color_branches_by_group=color_branches_by_group,
        y_axis_step=y_axis_step,
    )
    add_sample_group_legend(
        fig,
        sample_group_colors,
        show_labels=not no_axis and not color_branches_by_group,
        branch_colors=color_branches_by_group,
    )
    if no_axis:
        ax.axis('off')
    fig.savefig(out_path, bbox_inches='tight', dpi=200)
    plt.close(fig)
    return dendr, z


def heatmap_legend_path(heatmap_out):
    root, ext = os.path.splitext(heatmap_out)
    return root + '.legend' + (ext or '.png')


def plot_heatmap(m, out_path, legend_out_path, leaf_order, size=None, no_axis=False, cmap='coolwarm', vmin=None, vmax=None):
    n = m.shape[0]
    hs = size or max(int(round(0.3 * n)), 5)
    ordered = m.iloc[leaf_order, leaf_order]

    fig, ax = plt.subplots(figsize=(hs, hs))
    sns.heatmap(
        ordered,
        cmap=cmap,
        square=True,
        cbar=False,
        xticklabels=False,
        yticklabels=False,
        vmin=vmin,
        vmax=vmax,
        ax=ax,
    )
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(axis='both', which='both', length=0)
    ax.set_axis_off()
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    if no_axis:
        ax.axis('off')
    fig.savefig(out_path, bbox_inches='tight', pad_inches=0, dpi=300)
    plt.close(fig)

    legend_vmin = float(np.nanmin(ordered.values)) if vmin is None else vmin
    legend_vmax = float(np.nanmax(ordered.values)) if vmax is None else vmax
    legend_fig, legend_ax = plt.subplots(figsize=(4, 0.6))
    norm = Normalize(vmin=legend_vmin, vmax=legend_vmax)
    legend_fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap),
        cax=legend_ax,
        orientation='horizontal',
    )
    legend_fig.savefig(legend_out_path, bbox_inches='tight', dpi=300)
    plt.close(legend_fig)


def plot_support_dendrograms(
    support_dir,
    labels,
    threshold=None,
    no_axis=False,
    no_labels=False,
    width=None,
    height=None,
    height_per_y_unit=None,
    y_axis_step=None,
    sample_group_colors=None,
    color_branches_by_group=False,
):
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
        plot_single_dendrogram(
            m,
            out_path,
            threshold=threshold,
            no_axis=no_axis,
            no_labels=no_labels,
            width=width,
            height=height,
            height_per_y_unit=height_per_y_unit,
            y_axis_step=y_axis_step,
            sample_group_colors=sample_group_colors,
            color_branches_by_group=color_branches_by_group,
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw dendrogram and optional heatmap from one distance matrix, with optional branch support from support matrices.')
    parser.add_argument('matrix', help="TSV matrix with first column 'sample'.")
    parser.add_argument('--dendrogram-out', default='dendrogram.png', help='Dendrogram PNG [default: %(default)s]')
    parser.add_argument('--heatmap-out', default='heatmap.png', help='Heatmap PNG [default: %(default)s]')
    parser.add_argument('--clustermap-out', dest='heatmap_out', help=argparse.SUPPRESS)
    parser.add_argument('--heatmap-legend-out', default=None, help='Heatmap color legend PNG [default: <heatmap-out>.legend.png]')
    parser.add_argument('--support-dir', default=None, help='Directory with support matrices (*.tsv).')
    parser.add_argument('--support-min-height', type=float, default=0.0, help='Annotate support only for branches with y >= this value [default: %(default)s]')
    parser.add_argument('--support-fontsize', type=int, default=9, help='Font size for support labels.')
    parser.add_argument('--nwk', default=None, help='Write dendrogram tree to this Newick file. Internal nodes include support values when --support-dir is used.')
    parser.add_argument('--sample-group-colors', default=None, help="TSV with group names and RGB colors as group<TAB>R,G,B. Group names must match last parts of sample names split on '_' (sample.split('_')[-1]). Adds a color dot between each sample name and branch tip.")
    parser.add_argument('--color-branches-by-group', action='store_true', help='Use --sample-group-colors to color dendrogram branches instead of drawing sample group dots. Mixed-group branches are black and threshold colors are ignored.')
    parser.add_argument('--plot-support-dendrograms', action='store_true', help='Also draw a dendrogram for each support matrix and save it into that same directory.')
    parser.add_argument('--threshold', type=float, default=None, help='Dendrogram color threshold; default=0.7*max(height).')
    parser.add_argument('--skip-heatmap', action='store_true', help='Skip heatmap drawing.')
    parser.add_argument('--skip-clustermap', dest='skip_heatmap', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--no-axis', action='store_true', help='Hide axes.')
    parser.add_argument('--no-labels', action='store_true', help='Hide labels.')
    parser.add_argument('--dendrogram-width', type=int, default=None, help='Dendrogram width in inches.')
    parser.add_argument('--dendrogram-height', type=int, default=None, help='Dendrogram height in inches.')
    parser.add_argument('--dendrogram-height-per-y-unit', type=float, default=None, help='Set dendrogram height as this many inches per y-axis unit; ignored when --dendrogram-height is set. Example: value 10 with max height 0.5 gives a 5 inch figure.')
    parser.add_argument('--y-axis-step', type=float, default=None, help='Fixed spacing between y-axis tick labels; default chooses ticks automatically.')
    parser.add_argument('--heatmap-size', type=int, default=None, help='Heatmap size in inches.')
    parser.add_argument('--clustermap-size', dest='heatmap_size', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--heatmap-vmin', type=float, default=None, help='Heatmap color minimum; default is automatic.')
    parser.add_argument('--heatmap-vmax', type=float, default=None, help='Heatmap color maximum; default is automatic.')
    parser.add_argument('--heatmap-cmap', default='coolwarm', help='Matplotlib colormap for heatmap [default: %(default)s]')
    parser.add_argument('--legend', action='store_true', help='Show legend.')

    args = parser.parse_args()
    if args.heatmap_vmin is not None and args.heatmap_vmax is not None and args.heatmap_vmin >= args.heatmap_vmax:
        sys.exit('--heatmap-vmin must be smaller than --heatmap-vmax.')
    if args.color_branches_by_group and args.sample_group_colors is None:
        sys.exit('--color-branches-by-group requires --sample-group-colors.')
    if args.dendrogram_height_per_y_unit is not None and args.dendrogram_height_per_y_unit <= 0:
        sys.exit('--dendrogram-height-per-y-unit must be > 0.')
    if args.y_axis_step is not None and args.y_axis_step <= 0:
        sys.exit('--y-axis-step must be > 0.')

    m = read_matrix(args.matrix)
    sample_group_colors = None
    if args.sample_group_colors is not None:
        sample_group_colors = read_sample_group_colors(args.sample_group_colors)
        validate_sample_group_colors(m.index.tolist(), sample_group_colors, args.sample_group_colors)
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
        out_tsv = args.dendrogram_out[:-4] + '.support.tsv'
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
            plot_support_dendrograms(
                args.support_dir,
                m.index.tolist(),
                threshold=args.threshold,
                no_axis=args.no_axis,
                no_labels=args.no_labels,
                width=args.dendrogram_width,
                height=args.dendrogram_height,
                height_per_y_unit=args.dendrogram_height_per_y_unit,
                y_axis_step=args.y_axis_step,
                sample_group_colors=sample_group_colors,
                color_branches_by_group=args.color_branches_by_group,
            )

    n = m.shape[0]
    dw = args.dendrogram_width or max(int(round(0.15 * n)), 5)
    dh = args.dendrogram_height or (
        args.dendrogram_height_per_y_unit * z[:, 2].max()
        if args.dendrogram_height_per_y_unit is not None
        else int(round(dw / 3))
    )
    print('  Plotting dendrogram →', args.dendrogram_out)
    fig, ax = plt.subplots(figsize=(dw, dh))
    dendr = draw_dendrogram(
        ax,
        z,
        m.index.tolist(),
        thr,
        sample_group_colors=sample_group_colors,
        no_axis=args.no_axis,
        no_labels=args.no_labels,
        color_branches_by_group=args.color_branches_by_group,
        y_axis_step=args.y_axis_step,
    )
    if args.legend:
        print('  Adding sample group legend')
        add_sample_group_legend(
            fig,
            sample_group_colors,
            show_labels=not args.no_axis and not args.color_branches_by_group,
            branch_colors=args.color_branches_by_group,
        )
    if support_pct:
        annotate_support(ax, dendr, support_pct, selected_nodes, min_height=args.support_min_height, fontsize=args.support_fontsize)
    if args.no_axis:
        ax.axis('off')
    fig.savefig(args.dendrogram_out, bbox_inches='tight', dpi=200)
    plt.close(fig)

    save_clusters_file(args.dendrogram_out[:-4] + '.clusters', dendr, m, z, thr)
    print('  Saved cluster order →', args.dendrogram_out[:-4] + '.clusters')

    if args.nwk is not None:
        save_newick(args.nwk, z, m.index.tolist(), support_pct=support_pct)
        print('  Saved Newick tree →', args.nwk)

    if not args.skip_heatmap:
        print('  Plotting heatmap →', args.heatmap_out)
        heatmap_legend_out = args.heatmap_legend_out or heatmap_legend_path(args.heatmap_out)
        plot_heatmap(
            m,
            args.heatmap_out,
            heatmap_legend_out,
            dendr['leaves'],
            size=args.heatmap_size,
            no_axis=args.no_axis,
            cmap=args.heatmap_cmap,
            vmin=args.heatmap_vmin,
            vmax=args.heatmap_vmax,
        )
        print('  Saved heatmap legend →', heatmap_legend_out)
