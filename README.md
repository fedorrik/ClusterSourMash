# ClusterSourMash

ClusterSourMash makes UPGMA clustering of input sequences using kmer abundance based pairwise angular distance matrix. It doens't subsample kmers but takes all of them. Optianlly analog of Jackknife resampling might be applied by randomly subsampling 50% of kmers a lot of times to test branches support.

Config: in ./ClusterSourMash.sh specify path sourmash executable (line 12) and path to this dir with scripts (line 14)

Usage: ./ClusterSourMash.sh --help

Dependencies: [sourmash](https://github.com/sourmash-bio/sourmash), python libs [numpy, pandas, seaborn, matplotlib, scipy]

