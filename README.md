```markdown
# mclink: Metabolic Pathway Completeness and Abundance Analysis

## Overview

`mclink` provides comprehensive tools for analyzing metabolic pathway completeness and abundance using KEGG Orthology (KO) data from (meta)genomic and (meta)transcriptomic studies. Key features include:

- **Dual analysis modes**: Completeness (presence/absence) and abundance-weighted scoring
- **Flexible input**: Works with built-in KEGG references or custom datasets
- **Smart KO handling**: Specialized methods for plus-separated (subunits) and comma-separated (isoforms) KOs
- **Publication-ready outputs**: Pathway coverage metrics and detailed KO detection reports

## Installation

```r
# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("LiuyangLee/mclink")
```

## Quick Start

```r
library(mclink)

# Using built-in datasets
data(KO_pathway_ref)
data(KO_Sample_wide)

# Analyze selected pathways
selected_modules <- c("M00176", "M00165", "M00173", "M00374")
results <- mclink(
  ref = KO_pathway_ref[KO_pathway_ref$Module_Entry %in% selected_modules, ],
  data = KO_Sample_wide,
  table_feature = "completeness"
)

# Access results
head(results$coverage)      # Pathway coverage metrics
head(results$detected_KOs)  # Detected KOs per pathway
```

## Key Features

### Analysis Modes
- **Completeness analysis**: Binary presence/absence scoring
- **Abundance analysis**: Weighted by KO abundance levels

### Specialized KO Handling
| KO Type          | Scaling Methods                     | Typical Use Case          |
|------------------|-------------------------------------|---------------------------|
| Plus-separated   | mean (default), min (conservative), max (liberal) | Protein complexes/subunits |
| Comma-separated  | max (completeness), sum (abundance) | Gene isoforms/alternatives |

### Output Options
- **File exports**: TSV format for both coverage and KO detection
- **Pathway-split outputs**: Optional separate files per pathway

## Documentation

Full function reference:
```r
?mclink::mclink
```

## Citation

If you use `mclink` in your research, please cite:

> Li, L., Huang, D., Hu, Y., Rudling, N. M., Canniffe, D. P., Wang, F., & Wang, Y.
> "Globally distributed Myxococcota with photosynthesis gene clusters illuminate the origin and evolution of a potentially chimeric lifestyle."
> *Nature Communications* (2023), 14, 6450.
> https://doi.org/10.1038/s41467-023-42193-7

## Dependencies

- R (≥ 3.5)
- data.table (≥ 1.17.0)
- dplyr (≥ 1.1.4)
- stringr (≥ 1.5.1)
- tibble (≥ 3.2.1)

## License

GPL-3 © [Liuyang Li](https://orcid.org/0000-0001-6004-9437)

## Contact

- Maintainer: Liuyang Li <cyanobacteria@yeah.net>
- Bug reports: https://github.com/LiuyangLee/mclink/issues
```
