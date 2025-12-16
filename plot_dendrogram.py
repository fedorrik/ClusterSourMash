#!/usr/bin/env python3
import argparse, sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

parser = argparse.ArgumentParser(
    description="Dendrogram (and optional clustermap) from a pairwise distance matrix."
)
# I/O
parser.add_argument("matrix", help="TSV with 'sample' column + square matrix.")
parser.add_argument("-do", "--dendrogram-out", default="dendrogram.png", help="Dendrogram output file.")
parser.add_argument("-co", "--clustermap-out", default="clustermap.png", help="Clustermap output file.")
# Behavior
parser.add_argument("-t", "--threshold", type=float, default=None, help="color_threshold; default is 0.7*max(Z[:,2]).")
parser.add_argument("-s", "--skip-clustermap", action="store_true", help="Skip clustermap.")
parser.add_argument("-na", "--no-axis", action="store_true", help="Plot heatmap w/o axis.")

# Sizes
parser.add_argument("-dw", "--dendrogram-width", type=int, default=None, help="Dendrogram width (in).")
parser.add_argument("-dh", "--dendrogram-height", type=int, default=None, help="Dendrogram height (in).")
parser.add_argument("-cs", "--clustermap-size", type=int, default=None, help="Clustermap size WxH (in).")
args = parser.parse_args()

# Load + sanity
try:
    M = pd.read_csv(args.matrix, sep="\t").set_index("sample")
except Exception as e:
    sys.exit(f"Failed to read matrix: {e}")
if M.shape[0] != M.shape[1]: sys.exit("Matrix must be square.")
if (M.index != M.columns).any(): sys.exit("Row/column labels must match and be in same order.")
if (M.values.diagonal() != 0).any(): print("  Warning: diagonal not all zeros.", file=sys.stderr)

print("  Matrix shape:", M.shape)
n = M.shape[0]
dw = args.dendrogram_width or max(int(round(0.15 * n)), 5)
dh = args.dendrogram_height or int(round(dw / 3))
cs = args.clustermap_size or max(int(round(0.3 * n)), 5)

# Linkage
print("  Computing linkage (average)")
Z = linkage(squareform(M.values), method="average")
thr = args.threshold if args.threshold is not None else 0.7 * Z[:, 2].max()

# Dendrogram
print("  Plotting dendrogram →", args.dendrogram_out, f"(size: {dw}x{dh} in)")
plt.figure(figsize=(dw, dh))
dendrogram(Z, labels=M.index.tolist(), color_threshold=thr, leaf_rotation=90)
if args.no_axis:
	plt.gca().axis('off')
plt.savefig(args.dendrogram_out, bbox_inches="tight", dpi=100)
plt.close()

# Clustermap (optional)
if not args.skip_clustermap:
    print("  Plotting clustermap →", args.clustermap_out, f"(size: {cs}x{cs} in)")
    sns.clustermap(M, row_linkage=Z, col_linkage=Z, figsize=(cs, cs),
                   cmap="vlag", annot=False, fmt=".2f", annot_kws={"size": 0})
    plt.savefig(args.clustermap_out, bbox_inches="tight", dpi=100)
    plt.close()
print("  Done.")
