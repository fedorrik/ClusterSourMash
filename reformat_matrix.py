#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Reformat sourmash pairwise matrix in-place and convert distance to similarity (1 - d).")
    parser.add_argument("matrix_file", help="Path to the pairwise matrix file (tab-delimited) to rewrite in-place.")
    args = parser.parse_args()
    path = args.matrix_file

    # Load matrix
    M = pd.read_csv(path, sep=",")
    # Normalize sample names and set both columns and index
    sample_names = [col.split('/')[-1].replace('.fa', '') for col in M.columns]
    M.columns = sample_names
    M.index = sample_names
    M.index.name = 'sample'
    # Convert distance to similarity
    M = 1 - M
    # Rewrite the same file
    M.to_csv(path, sep="\t")

if __name__ == "__main__":
    main()
