![GitHub Installs](https://img.shields.io/endpoint?url=https://r-pkg.github.io/install-stats/LiuyangLee/mclink/badge.json&label=GitHub%20Installs&style=flat-square)
![CRAN Downloads](https://img.shields.io/badge/dynamic/json?url=https://cranlogs.r-pkg.org/badges/grand-total/mclink&query=$.count&label=CRAN%20Downloads&color=blue&style=flat-square)

# mclink: Metabolic Pathway Completeness and Abundance Analysis

## Overview
`mclink` provides comprehensive tools for analyzing metabolic pathway completeness and abundance using KEGG Orthology (KO) data from (meta)genomic and (meta)transcriptomic studies. Key features include:
- **Dual analysis modes**: Completeness (presence/absence) and abundance-weighted scoring
- **Flexible input**: Works with built-in KEGG references or custom datasets
- **Smart KO handling**: Specialized methods for plus-separated (subunits) and comma-separated (isoforms) KOs
- **Publication-ready outputs**: Pathway coverage metrics and detailed KO detection reports

## Description
The distill analysis of metabolic pathway coverage is calculated based on the abundance or presence of KOs in a given module, as per the KEGG Module Definition. In detail, the coverage of a KEGG Module is determined by first dividing a set of KOs into distinct steps. The coverage for each step is then calculated separately, and summarized as the coverage of this KEGG Module. When calculating coverage, spaces or plus signs connecting KO numbers are interpreted as AND operators, while commas are interpreted as OR operators. For instance, to calculate the completeness for the module M00020: “K00058 K00831 (K01079, K02203, K22305)”:
- **1**: Convert the abundance table of KOs into a 0-1 matrix.
- **2**: Consider K00058 as step 1, K00831 as step 2, and (K01079, K02203, K22305) as step 3.
- **3**: Calculate the `maximum` (or `minimum`, `mean`) value of step 3 (K01079, K02203, K22305). If presence of any KO indicates completeness, calculate the `maximum` value; If all KOs must be resent for completeness, calculate the `minimum` value. For a moderate approach, calculate the `mean` value.
- **4**: Use this value along with the values of step 1 and step 2 to calculate the `mean` value, representing the completeness of the module M00020.


## Key Features
### Analysis Modes
- **Completeness analysis**: Binary presence/absence scoring
- **Abundance analysis**: Weighted by KO abundance levels

### Specialized KO Handling
| KO Type          | Scaling Method       | Description                                                                 | Typical Use Case          |
|------------------|----------------------|-----------------------------------------------------------------------------|---------------------------|
| Plus-separated   | mean (default)       | Moderate approach - calculates average value of all components              | Protein complexes         |
| (K1+K2+...)      | min (conservative)   | Strict requirement - uses lowest value (all components must be present)     | Enzyme subunits           |
|                  | max (liberal)        | Lenient approach - uses highest value (any component indicates completeness)|                           |
|                  |                      |                                                                             |                           |
| Comma-separated  | max                  | Completeness assessment - (any component indicates completeness)            | Gene isoforms             |
| (K1,K2,...)      | sum                  | Abundance quantification - sums all functionally equivalent variants        | Alternative pathways      |


### Output Options
- **R list**: Coverage and KO detection results as dataframes, along with a log string recording the analysis process.
- **File exports**: Coverage and KO detection results in TSV format, along with a log file recording the analysis process.
- **Pathway-split exports**: Optional per-pathway TSV files (if enabled).


## Installation
```R
# Install from CRAN
install.packages("mclink")

# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("LiuyangLee/mclink")
```

## Quick Start
```R
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
head(results$log)           # log
```


## Example for Input and Output Data
### 1 Input Data Preview (R dataframes `KO_Sample_wide` and `KO_pathway_ref`)
#### 1.1 Genome/Sample KO abundance/transcript/prevalence dataframe (e.g., head(`KO_Sample_wide`))
| KO     | Marinobacter_salarius | Pseudooceanicola_nanhaiensis | Alteromonas_australica | Henriciella_pelagia |
|--------|----------------------|----------------------------|----------------------|-------------------|
| K00001 | 61.44954            | 16.92329                  | 6.854643            | 9.592472         |
| K00002 | 0.00000             | 0.00000                   | 5.983655            | 0.000000         |
| K00014 | 49.84410            | 20.23131                  | 17.343312           | 24.083303        |
| K00015 | 0.00000             | 27.09820                  | 0.000000            | 0.000000         |
| K00018 | 43.10113            | 0.00000                   | 19.115001           | 3.344125         |
| K00019 | 49.84410            | 41.88975                  | 13.879528           | 9.634129         |
#### 1.2 (Optional) Pathway information from built-in or user-defined dataframe (e.g., head(`KO_pathway_ref`))
| Orthology_Entry | Module_Type       | Level_2           | Level_3         | Module_Entry | Module_Name                                      | Definition                                                                                                                                                                                                                                                                                                                                 | Orthology_Symbol | Orthology_Name                                                                 | KO_Symbol          |
|-----------------|-------------------|-------------------|-----------------|--------------|---------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------|-------------------------------------------------------------------------------|-------------------|
| K00855          | Pathway modules   | Energy metabolism | Carbon fixation | M00165       | Reductive pentose phosphate cycle (Calvin cycle) | K00855 (K01601-K01602) K00927 (K05298,K00150,K00134) K01803 (K01623,K01624) (K03841,K02446,K11532,K01086) K00615 (K01100,K11532,K01086) (K01807,K01808) K01783                                                                                                             | PRK, prkB        | phosphoribulokinase [EC:2.7.1.19]                                             | K00855; PRK, prkB |
| K01601          | Pathway modules   | Energy metabolism | Carbon fixation | M00165       | Reductive pentose phosphate cycle (Calvin cycle) | K00855 (K01601-K01602) K00927 (K05298,K00150,K00134) K01803 (K01623,K01624) (K03841,K02446,K11532,K01086) K00615 (K01100,K11532,K01086) (K01807,K01808) K01783                                                                                                             | rbcL, cbbL       | ribulose-bisphosphate carboxylase large chain [EC:4.1.1.39]                   | K01601; rbcL, cbbL |
| K01602          | Pathway modules   | Energy metabolism | Carbon fixation | M00165       | Reductive pentose phosphate cycle (Calvin cycle) | K00855 (K01601-K01602) K00927 (K05298,K00150,K00134) K01803 (K01623,K01624) (K03841,K02446,K11532,K01086) K00615 (K01100,K11532,K01086) (K01807,K01808) K01783                                                                                                             | rbcS, cbbS       | ribulose-bisphosphate carboxylase small chain [EC:4.1.1.39]                   | K01602; rbcS, cbbS |

### 2 Output Data Preview (R list `mc_list`)
#### 2.1 Pathway Coverage Matrix for Genomes or Samples (`mc_list$coverage`)
| Module_Entry | Level_2           | Level_3         | Pathway Name                                      | Definition (Key KOs)                                                                                                                                                                                                                                                                                                                                 | Marinobacter salarius | Pseudooceanicola nanhaiensis | Alteromonas australica | Henriciella pelagia |
|--------------|-------------------|-----------------|---------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------|------------------------------|------------------------|---------------------|
| M00165       | Energy metabolism | Carbon fixation | Reductive pentose phosphate cycle (Calvin cycle) | K00855 (K01601-K01602) K00927 (K05298,K00150,K00134) K01803 (K01623,K01624) (K03841,K02446,K11532,K01086) K00615 (K01100,K11532,K01086) (K01807,K01808) K01783                                                                                                                                                                                     | 0.4545455             | 0.6363636                    | 0.6363636             | 0.6363636          |
| M00173       | Energy metabolism | Carbon fixation | Reductive citrate cycle (Arnon-Buchanan cycle)   | (K00169+K00170+K00171+K00172,K03737) ((K01007,K01006) K01595,K01959+K01960,K01958) K00024 (K01676,K01679,K01677+K01678) (K00239+K00240-K00241-K00242,K18556+K18557+K18558+K18559+K18560) (K01902+K01903) (K00174+K00175-K00177-K00176) K00031 (K01681,K27802,K01682) (K15230+K15231,K15232+K15233 K15234)                                        | 0.5000000             | 0.7000000                    | 0.5000000             | 0.5000000          |
| M00376       | Energy metabolism | Carbon fixation | 3-Hydroxypropionate bi-cycle                     | (K02160+K01961+K01962+K01963) K14468 K14469 K15052 K05606 (K01847,K01848+K01849) (K14471+K14472) (K00239+K00240+K00241) K01679                                                                                                                                                                                                                      | 0.2222222             | 0.5555556                    | 0.0000000             | 0.3333333          |
| M00375       | Energy metabolism | Carbon fixation | Hydroxypropionate-hydroxybutylate cycle          | (K01964+K15037+K15036) K15017 K15039 K15018 K15019 K15020 K05606 (K01848+K01849) (K15038,K15017) K14465 (K14466,K18861,K25774) K14534 K15016 K00626                                                                                                                                                                                                 | 0.0000000             | 0.1428571                    | 0.0000000             | 0.1428571          |
| M00374       | Energy metabolism | Carbon fixation | Dicarboxylate-hydroxybutyrate cycle              | (K00169+K00170+K00171+K00172) K01007 K01595 K00024 (K01677+K01678) (K00239+K00240-K00241-K18860) (K01902+K01903) (K15038,K15017) K14465 (K14467,K18861,K25774) K14534 K15016 K00626                                                                                                                                                                | 0.2307692             | 0.3076923                    | 0.2307692             | 0.3846154          |
| M00377       | Energy metabolism | Carbon fixation | Reductive acetyl-CoA pathway (Wood-Ljungdahl)    | K00198 (K05299-K15022,K22015+K25123+K25124) K01938 K01491-K01500 K00297-K25007-K25008 K15023 K14138+K00197+K00194                                                                                                                                                                                                                                   | 0.4285714             | 0.4285714                    | 0.2857143             | 0.2857143          |

##### Interpretation Guide
- *Autotrophic potential:* Genomes with >80% completeness in carbon fixation pathways may represent autotrophs. None of the tested genomes met the 80% completeness threshold for canonical carbon fixation pathways, suggesting they likely utilize heterotrophic metabolism.
- *Verification:* Always check presence of pathway-specific marker genes (e.g., RuBisCO for Calvin cycle, aclAB for rTCA cycle)
- *Threshold adjustment:* Some pathways have alternative enzyme implementations (e.g., Form II RuBisCO) requiring lower completeness thresholds

#### 2.2 Detected KOs per Genome/Sample (`mc_list$detected_KOs`)
| Module_Entry | Pathway Name                                      | Marinobacter salarius                        | Pseudooceanicola nanhaiensis                  | Alteromonas australica                     | Henriciella pelagia               |
|--------------|---------------------------------------------------|----------------------------------------------|-----------------------------------------------|--------------------------------------------|-----------------------------------|
| M00165       | Calvin cycle                                      | K00615 K00855 K01783 K01803 K01807           | K00134 K00615 K01623 K01783 K01803 K01808...  | K00134 K00615 K00855 K01623 K01783 K01803 | K00615 K01623 K01783 K01803...    |
| M00173       | Arnon-Buchanan cycle                              | K00031 K00239-K00242 K01007 K01595 K01676... | K00024 K00239-K00242 K01006 K01679 K01902...  | K00024 K00239 K00241 K01007 K01595 K01676 | K00024 K00239-K00242 K01006...    |
| M00376       | 3-Hydroxypropionate bi-cycle                      | K00239-K00241 K01847 K01962-K01963 K02160    | K00239-K00241 K01679 K01847 K01961-K01963...  | K00239 K00241 K01962-K01963 K02160        | K00239-K00241 K01847 K01962...    |
| M00375       | Hydroxypropionate-hydroxybutylate cycle           | K01848                                      | K00626 K05606                                 | -                                          | K00626 K05606                     |
| M00374       | Dicarboxylate-hydroxybutyrate cycle               | K00239-K00241 K01007 K01595 K01902          | K00024 K00239-K00241 K00626 K01902-K01903    | K00024 K00239 K00241 K01007 K01595 K01902 | K00024 K00239-K00241 K00626...    |
| M00377       | Wood-Ljungdahl pathway                            | K00297 K01491 K01938 K22015                 | K00297 K01491 K01938                         | K00297 K01491                              | K00297 K01491                     |

#### 2.3 Process Log (`mc_list$log`)
```
[1] "[2025-08-21 17:35:12] mclink started!"
[2] "[2025-08-21 17:35:12] Input Sample-KO table type: completeness"
[3] "[2025-08-21 17:35:12] Scale method for plus: min"
[4] "[2025-08-21 17:35:12] Scale method for comma: max"
[5] "[2025-08-21 17:35:12] Pathway information dataframe successfully imported."
[6] "[2025-08-21 17:35:12] There are 7 Modules in the Pathway information dataframe: M00165 M00173 M00376 M00375 M00374 M00377 M00176"
[7] "[2025-08-21 17:35:12] Genome-KO File successfully imported."
[8] "[2025-08-21 17:35:12] There are 50 intersect KOs between the KO list and input dataframe: K00024 K00031 K00134 K00239 K00240 K00241 K00242 K00297 K00380 K00381 K00390 K00615 K00626 K00855 K00860 K00955 K00956 K00957 K00958 K01006 K01007 K01491 K01595 K01623 K01676 K01679 K01681 K01682 K01783 K01803 K01807 K01808 K01847 K01848 K01902 K01903 K01938 K01958 K01960 K01961 K01962 K01963 K02160 K02446 K03841 K05606 K08691 K09709 K11532 K22015"
[9] "[2025-08-21 17:35:12] Converting abundance table to completeness table..."
[10] "[2025-08-21 17:35:12] Starting Module: M00165"
[11] "[2025-08-21 17:35:12] Processing module steps: K00855 (K01601-K01602) K00927 (K05298,K00150,K00134) K01803 (K01623,K01624) (K03841,K02446,K11532,K01086) K00615 (K01100,K11532,K01086) (K01807,K01808) K01783"
[12] "[2025-08-21 17:35:12] After omitting minus KOs: K00855 K01601 K00927 (K05298,K00150,K00134) K01803 (K01623,K01624) (K03841,K02446,K11532,K01086) K00615 (K01100,K11532,K01086) (K01807,K01808) K01783"
[13] "[2025-08-21 17:35:12] Nested steps include: "
[14] "[2025-08-21 17:35:12]     (K05298,K00150,K00134)"
[15] "[2025-08-21 17:35:12]     (K01623,K01624)"
[16] "[2025-08-21 17:35:12]     (K03841,K02446,K11532,K01086)"
[17] "[2025-08-21 17:35:12]     (K01100,K11532,K01086)"
[18] "[2025-08-21 17:35:12]     (K01807,K01808)"
[19] "[2025-08-21 17:35:12] Start processing nested steps..."
[20] "[2025-08-21 17:35:12] Analyzing M00165_1: (K05298,K00150,K00134)"
[21] "[2025-08-21 17:35:12] Bracket level: 1"
[22] "[2025-08-21 17:35:12]    Running KOs comma: M00165_1 = K05298"
[23] "[2025-08-21 17:35:12]    Running KOs comma: M00165_1 = K00150"
[24] "[2025-08-21 17:35:12]    Running KOs comma: M00165_1 = K00134"
[25] "[2025-08-21 17:35:12] Analyzing M00165_2: (K01623,K01624)"
[26] "[2025-08-21 17:35:12] Bracket level: 2"
[27] "[2025-08-21 17:35:12]    Running KOs comma: M00165_2 = K01623"
[28] "[2025-08-21 17:35:12]    Running KOs comma: M00165_2 = K01624"
[29] "[2025-08-21 17:35:12] Analyzing M00165_3: (K03841,K02446,K11532,K01086)"
[30] "[2025-08-21 17:35:12] Bracket level: 3"
[31] "[2025-08-21 17:35:12]    Running KOs comma: M00165_3 = K03841"
[32] "[2025-08-21 17:35:12]    Running KOs comma: M00165_3 = K02446"
[33] "[2025-08-21 17:35:12]    Running KOs comma: M00165_3 = K11532"
[34] "[2025-08-21 17:35:12]    Running KOs comma: M00165_3 = K01086"
[35] "[2025-08-21 17:35:12] Analyzing M00165_4: (K01100,K11532,K01086)"
[36] "[2025-08-21 17:35:12] Bracket level: 4"
[37] "[2025-08-21 17:35:12]    Running KOs comma: M00165_4 = K01100"
[38] "[2025-08-21 17:35:12]    Running KOs comma: M00165_4 = K11532"
[39] "[2025-08-21 17:35:12]    Running KOs comma: M00165_4 = K01086"
[40] "[2025-08-21 17:35:12] Analyzing M00165_5: (K01807,K01808)"
[41] "[2025-08-21 17:35:12] Bracket level: 5"
[42] "[2025-08-21 17:35:12]    Running KOs comma: M00165_5 = K01807"
[43] "[2025-08-21 17:35:12]    Running KOs comma: M00165_5 = K01808"
[44] "[2025-08-21 17:35:12] Start processing final step..."
[45] "[2025-08-21 17:35:12] Analyzing M00165: K00855 K01601 K00927 M00165_1 K01803 M00165_2 M00165_3 K00615 M00165_4 M00165_5 K01783"
[46] "[2025-08-21 17:35:12]    Running KOs space: M00165 = K00855"
[47] "[2025-08-21 17:35:12]    Running KOs space: M00165 = K01601"
[48] "[2025-08-21 17:35:12]    Running KOs space: M00165 = K00927"
[49] "[2025-08-21 17:35:12]    Running KOs space: M00165 = M00165_1"
[50] "[2025-08-21 17:35:12]    Running KOs space: M00165 = K01803"
[51] "[2025-08-21 17:35:12]    Running KOs space: M00165 = M00165_2"
[52] "[2025-08-21 17:35:12]    Running KOs space: M00165 = M00165_3"
[53] "[2025-08-21 17:35:12]    Running KOs space: M00165 = K00615"
[54] "[2025-08-21 17:35:12]    Running KOs space: M00165 = M00165_4"
[55] "[2025-08-21 17:35:12]    Running KOs space: M00165 = M00165_5"
[56] "[2025-08-21 17:35:12]    Running KOs space: M00165 = K01783"
[57] "[2025-08-21 17:35:12] Completed Module: M00165"
```


## Documentation
Full function reference:
```R
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
