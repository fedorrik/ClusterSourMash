# ClusterSourMash

ClusterSourMash performs UPGMA clustering of input sequences using a k-mer abundance–based pairwise angular distance matrix. It does not subsample k-mers, but instead uses all of them. Optionally, a jackknife-like resampling procedure can be applied by repeatedly subsampling 50% of k-mers to estimate branch support.

## Configuration

In `./ClusterSourMash.sh`, specify:
- the path to the `sourmash` executable (line 12)
- the path to this directory containing the scripts (line 14)

## Usagege

```bash
./ClusterSourMash.sh --help
```
Dependencies: 
- [sourmash](https://github.com/sourmash-bio/sourmash)
- Python libraries [numpy, pandas, seaborn, matplotlib, scipy]
