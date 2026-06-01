# iota-implication-graphs

Reproducibility code for the paper:

> Noel, Y. (2026). On dimensional implication graphs. Psychometrika, 1–29. doi:10.1017/psy.2026.10123

This repository contains the R code and TikZ source needed to reproduce the
empirical illustration of the paper (the analysis of the Inductive Reasoning
Developmental Test, TDRI).

The $\iota$-index is a measure of asymmetric statistical dependence between
binary events. It decomposes into an association component (the log odds
ratio) and an asymmetry component, supports two-step gatekeeping inference for
directed dependence, and yields a Rasch-metric interpretation that licenses a
dimensional analysis of the resulting Dimensional Implication Graph (DIG).

## Repository layout

```
iota-implication-graphs/
├── R/
│   ├── iota17.R         # Core DiG R6 class and methods
│   ├── TDRI.R           # End-to-end analysis of the TDRI dataset
│   └── iota_plots.R     # Plotting helpers (Figures 5 and 6 of the paper)
├── tikz/
│   └── TDRI-full-meta-graph-iota16-unrotated.tex   # TikZ source for Figure 6
├── LICENSE
└── README.md
```

## Dependencies

R (>= 4.0) with the following packages from CRAN:

- `R6`, `igraph`, `diagram`, `MASS`, `Matrix`
- `mvtnorm`, `mclust`, `mirt`
- `blockmodels`, `greed` (Stochastic Block Model fitting)

Install all of them in one call:

```r
install.packages(c("R6","igraph","diagram","MASS","Matrix",
                   "mvtnorm","mclust","mirt",
                   "blockmodels","greed"))
```

## Data

The TDRI dataset (Golino & Gomes, 2015) is **not** redistributed in this
repository. It is shipped with the
[`EGAnet`](https://CRAN.R-project.org/package=EGAnet) R package; the easiest
way to obtain a local copy in the format the scripts expect is:

```r
install.packages("EGAnet")
data(TDRI, package = "EGAnet")
write.csv2(TDRI, "TDRI.csv", row.names = FALSE)
```

Place the resulting `TDRI.csv` in the `R/` directory before running the
scripts.

## Reproducing the analysis

From the `R/` directory:

```r
setwd("R")        # if not already there
source("TDRI.R")  # runs the full analysis pipeline
```

`TDRI.R` will:

1. load the data and define the theoretical Q-matrix (7 developmental stages);
2. instantiate the `DiG` class from `iota17.R` and compute the $\iota^*$ index
   for every item pair (Section 4 of the paper);
3. apply the two-step gatekeeping inference (one-sided OR test, then two-sided
   directionality test on $\Delta$, both BY-corrected);
4. fit a Stochastic Block Model on the directed graph to recover developmental
   clusters;
5. relabel SBM cluster IDs to align with the theoretical Q-matrix
   (`remap_clusters_to_stages()` in `iota_plots.R`);
6. produce Figures 5 and 6 of the paper.

## Notes

- The SBM step in `greed`/`blockmodels` returns clusters with arbitrary
  integer IDs; `iota_plots.R::remap_clusters_to_stages()` realigns them to
  the theoretical developmental order using each cluster's modal Q-block.
  This relabeling also permutes the cached `metagraph$adjacency` and
  `metagraph$centroids`, ensuring that meta-arrows on the cluster-level graph
  point in the correct direction.
- Figure 6 is rendered from `tikz/TDRI-full-meta-graph-iota16-unrotated.tex`.
  Compile with `pdflatex` (TikZ + a recent TeX Live distribution).

## License

MIT — see `LICENSE`.

## Citation

If you use this code, please cite:

```bibtex
@article{Noel2026iota,
  author  = {Yvonnick Noel},
  title   = {On dimensional implication graphs},
  journal = {Psychometrika},
  year    = {2026},
  note    = {Under review}
}
```

## Contact

Yvonnick Noel — *yvonnick.noel@univ-rennes2.fr*
Department of Psychology, University of Rennes 2, France
