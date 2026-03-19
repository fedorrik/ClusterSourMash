#!/bin/bash

#SBATCH --job-name=ClusterSourMash
#SBATCH --partition=medium
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

set -euo pipefail

sourmash_exec=${SOURMASH_EXEC:-sourmash}
#script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script_dir=/private/home/fryabov/soft/ClusterSourMash

usage() {
    cat <<USAGE
Usage:
  $0 --input-dir DIR --suffix SUFFIX [options]

Required:
  -i, --input-dir DIR          Directory with FASTA files
  -x, --suffix SUFFIX          File suffix to match and strip, e.g. .fa or .fasta

Optional:
  -j, --jobs INT               Parallel jobs [default: 8]
  -n, --n-support INT          Number of support replicates [default: 0]
  -s, --test-scaled INT        scaled value for support replicate sketches [default: 2]
  -t, --support-min-hight INT      Annotate support only for branches with y >= this value [default: 0.0]
  -p, --plot-support-dendrograms
                               Also draw a dendrogram for every support matrix
  -h, --help                   Show this help

Outputs always:
  pairwise_matrix.txt
  dendrogram.png
  clustermap.png

If --n-support > 0, also writes:
  sourmash_support_matrices/rep_XXXX.tsv
  support_matrix_comparison/
    support_matrix_variation.mean.tsv
    support_matrix_variation.std.tsv
    support_matrix_variation.cv.tsv
    support_matrix_variation.std.clustermap.png
    support_matrix_variation.cv.clustermap.png
  dendrogram.png.support.tsv
USAGE
}

input_dir=""
filename_ending=""
N_JOBS=8
N_SUPPORT=0
TEST_SCALED=2
SUPPORT_MIN_HIGHT=0
PLOT_SUPPORT_DENDROGRAMS=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input-dir)
            input_dir="$2"; shift 2 ;;
        -x|--suffix)
            filename_ending="$2"; shift 2 ;;
        -j|--jobs)
            N_JOBS="$2"; shift 2 ;;
        -n|--n-support)
            N_SUPPORT="$2"; shift 2 ;;
        -s|--test-scaled)
            TEST_SCALED="$2"; shift 2 ;;
        -t|--support-min-hight)
            SUPPORT_MIN_HIGHT="$2"; shift 2 ;;
        -p|--plot-support-dendrograms)
            PLOT_SUPPORT_DENDROGRAMS=1; shift ;;
        -h|--help)
            usage; exit 0 ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 1 ;;
    esac
done

if [[ -z "$input_dir" || -z "$filename_ending" ]]; then
    usage
    exit 1
fi
if [[ ! -d "$input_dir" ]]; then
    echo "Error: input directory not found: $input_dir" >&2
    exit 1
fi
if ! [[ "$N_JOBS" =~ ^[0-9]+$ ]] || [[ "$N_JOBS" -lt 1 ]]; then
    echo "Error: --jobs must be a positive integer." >&2
    exit 1
fi
if ! [[ "$N_SUPPORT" =~ ^[0-9]+$ ]] || [[ "$N_SUPPORT" -lt 0 ]]; then
    echo "Error: --n-support must be a non-negative integer." >&2
    exit 1
fi
if ! [[ "$TEST_SCALED" =~ ^[0-9]+$ ]] || [[ "$TEST_SCALED" -lt 1 ]]; then
    echo "Error: --test-scaled must be a positive integer." >&2
    exit 1
fi
if [[ -n "$SUPPORT_MIN_HIGHT" ]] && ! awk -v x="$SUPPORT_MIN_HIGHT" '
       BEGIN {
           exit !(x ~ /^([0-9]+(\.[0-9]*)?|\.[0-9]+)$/ && x >= 0 && x <= 1)
       }
   '
then
    echo "Error: --support-min-hight must be a number from 0 to 1." >&2
    exit 1
fi

for req in "$script_dir/reformat_matrix.py" "$script_dir/plot_dendrogram.py"; do
    if [[ ! -f "$req" ]]; then
        echo "Error: required script not found: $req" >&2
        exit 1
    fi
done
if [[ "$N_SUPPORT" -gt 0 && ! -f "$script_dir/compare_support_matrices.py" ]]; then
    echo "Error: required script not found: $script_dir/compare_support_matrices.py" >&2
    exit 1
fi

mapfile -d '' fasta_files < <(find "$input_dir" -maxdepth 1 -type f -name "*${filename_ending}" -print0 | sort -z)
if [[ ${#fasta_files[@]} -eq 0 ]]; then
    echo "Error: no files matching *${filename_ending} found in $input_dir" >&2
    exit 1
fi

main_sig_dir="sourmash_sig_scaled1"
support_dir="sourmash_support_matrices"
comparison_dir="support_matrix_comparison"
main_matrix="pairwise_matrix.txt"

rm -rf "$main_sig_dir" "$support_dir" "$comparison_dir" "$main_matrix"
mkdir -p "$main_sig_dir"
#rm -f  dendrogram.png clustermap.png dendrogram.png.clusters dendrogram.png.support.tsv

echo "Found ${#fasta_files[@]} input files"
echo "Sketching main signatures with scaled=1 (jobs: $N_JOBS)"
export sourmash_exec main_sig_dir filename_ending
printf '%s\0' "${fasta_files[@]}" | \
  xargs -0 -P "$N_JOBS" -n 1 bash -c '
    file="$1"
    base_name=$(basename "$file")
    base_name=${base_name%"$filename_ending"}
    echo "  -> $base_name"
    "$sourmash_exec" sketch dna -f -p "k=21,scaled=1,abund" -o "$main_sig_dir/$base_name.sig" "$file" 2>> sourmash.log
  ' _

echo "Calculating pairwise matrix from scaled=1 signatures"
"$sourmash_exec" compare "$main_sig_dir"/*.sig --csv "$main_matrix" >> sourmash.log 2>> sourmash.log

echo "Reformatting main matrix"
python3 "$script_dir/reformat_matrix.py" pairwise_matrix.txt ${filename_ending}

plot_cmd=(python3 "$script_dir/plot_dendrogram.py" "$main_matrix" \
  --dendrogram-out dendrogram.png \
  --clustermap-out clustermap.png)

if [[ "$N_SUPPORT" -gt 0 ]]; then
    mkdir -p "$support_dir"
    echo "Running $N_SUPPORT support replicate analyses with scaled=$TEST_SCALED"
    export TEST_SCALED filename_ending sourmash_exec

    for rep in $(seq 1 "$N_SUPPORT"); do
        seed=$((4342 + rep))
        rep_tag=$(printf 'rep_%04d' "$rep")
        rep_matrix="$support_dir/${rep_tag}.tsv"
        rep_sig_dir=$(mktemp -d "${TMPDIR:-/tmp}/sourmash_${rep_tag}_XXXXXX")

        echo "[$rep_tag] sketching with seed=$seed"
        export rep_sig_dir seed
        printf '%s\0' "${fasta_files[@]}" | \
          xargs -0 -P "$N_JOBS" -n 1 bash -c '
            file="$1"
            base_name=$(basename "$file")
            base_name=${base_name%"$filename_ending"}
            "$sourmash_exec" sketch dna -f -p "k=21,scaled=$TEST_SCALED,seed=$seed,abund" -o "$rep_sig_dir/$base_name.sig" "$file" 2>> sourmash.log
          ' _

        echo "[$rep_tag] comparing signatures"
        "$sourmash_exec" compare "$rep_sig_dir"/*.sig --csv "$rep_matrix" >> sourmash.log 2>> sourmash.log

        echo "[$rep_tag] reformatting matrix"
        python3 "$script_dir/reformat_matrix.py" "$rep_matrix" ${filename_ending}

        rm -rf "$rep_sig_dir"
    done

    echo "Comparing support matrices"
    python3 "$script_dir/compare_support_matrices.py" \
      --bootstrap-dir "$support_dir" \
      --reference-matrix "$main_matrix" \
      --out-dir "$comparison_dir" \
      --prefix support_matrix_variation \
      --mask-diagonal

    plot_cmd+=(--support-dir "$support_dir")
    if [[ -n "$SUPPORT_MIN_HIGHT" ]]; then
        plot_cmd+=(--support-min-height "$SUPPORT_MIN_HIGHT")
    fi
    if [[ "$PLOT_SUPPORT_DENDROGRAMS" -eq 1 ]]; then
        plot_cmd+=(--plot-support-dendrograms)
    fi
else
    echo "--n-support not provided or zero: main matrix only"
fi

echo "Plotting main dendrogram and clustermap"
"${plot_cmd[@]}"

rm -r sourmash_sig_scaled1/

echo "All done"
