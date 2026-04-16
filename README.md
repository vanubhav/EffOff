# EffOff: Efficacy of Offset Forestry

**Version**: 1.0.0

## Overview
`EffOff` automates the retrieval and robust statistical analysis of vegetation indices (EVI) from Google Earth Engine. It harmonizes Landsat imagery, performs advanced trend analysis utilizing conditional `BFAST` breakpoints, `ShapeSelectForest` geometric algorithms, and `modifiedmk` autocorrelation-corrected Mann-Kendall statistics, while intelligently segmenting evaluation across distinctly defined lifecycles (Pre-management, Management, and Post-management).

## Installation
You can install the development version of EffOff from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("yourusername/EffOff")
```

## Features
- **Data Extraction**: Automated full-history EVI retrieval from GEE, dynamically gap-filling overlapping Landsat sensor networks.
- **Trend Analysis**:
    - Project-aware temporal partitioning (Pre/During/Post Management transitions).
    - Autocorrelation-corrected Mann-Kendall calculations via `modifiedmk` to rigorously extract true `Tau` and `Sen's Slopes` over temporally correlated environmental variables.
    - Dynamic adaptive `BFAST` deployment that gracefully shifts mathematical configurations (`y ~ 1`, `y ~ t`, and harmonic) based on pixel time-series limits natively avoiding sparse-regression constraints.
    - Vectorized Concordance modeling across structural vs continuous evaluations natively (i.e. `Agreement_Triple`).
- **Visualization**: Generates custom `ggplot2` semantic palettes out-of-the-box and transparently renders interactive HTML `trend_report` profiles encompassing automated summary histograms, diagnostic regression cross-evaluations, and phase transition graphs flawlessly.

## License
MIT
