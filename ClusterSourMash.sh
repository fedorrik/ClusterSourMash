#!/bin/bash

#SBATCH --job-name=ClusterSourMash
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8


set -euo pipefail

# Path to sourmash executable (change if needed)
sourmash_exec=sourmash
script_dir=/private/home/fryabov/soft/ClusterSourMash

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input_dir_with_fasta> [N_JOBS]"
    exit 1
fi

# Directory with FASTA files (.fa)
input_dir=$1
# Number of parallel jobs (default: 8)
N_JOBS=${2:-8}

# Directory where this script lives (for Python scripts)
# Doesn't work with slurm
# script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

output_dir="output"
mkdir -p "$output_dir"

echo "Creating sourmash abundance-aware signatures in parallel (jobs: $N_JOBS)..."

# Find all .fa files and process them in parallel
find "$input_dir" -maxdepth 1 -type f -name '*.fa' -print0 \
  | xargs -0 -P "$N_JOBS" -I {} bash -c '
      file="$1"
      sourmash_exec="$2"
      output_dir="$3"

      base_name=$(basename "$file" .fa)
      echo "  -> $base_name"

      "$sourmash_exec" sketch dna -f -p k=21,scaled=1,abund \
        -o "$output_dir/$base_name.sig" "$file" 2>/dev/null
    ' _ {} "$sourmash_exec" "$output_dir"

echo "Calculating pairwise distances with sourmash (abundance-aware)..."

"$sourmash_exec" compare "$output_dir"/*.sig --csv pairwise_matrix.txt > /dev/null 2> /dev/null

echo "Pairwise comparison completed. Matrix saved to pairwise_matrix.txt"

# Reformat pairwise_matrix.txt (in-place, via your script)
echo "Reformatting matrix..."
python3 "$script_dir/reformat_matrix.py" pairwise_matrix.txt
echo "Matrix reformatted and saved to pairwise_matrix.txt"

# Plot dendrogram
echo "Plotting dendrogram..."
python3 "$script_dir/plot_dendrogram.py" pairwise_matrix.txt

# Cleanup temporary signatures
rm -r "$output_dir"

echo "All done."
