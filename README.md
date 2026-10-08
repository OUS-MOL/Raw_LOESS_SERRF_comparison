# Raw vs LOESS vs SERRF normalization comparison (untargeted LC-MS)

R workflow that compares **uncorrected**, **LOESS-corrected** and **SERRF-corrected**
peak areas from Compound Discoverer 3.4 exports across multiple analytical batches,
using pooled QC (PQC) injections to drive the correction and independent reference
materials (QC/PQC, a long-term reference material "Oslo-RM"/LTQC, and NIST SRM 1950)
to assess it.

## What the script does

1. Imports peak areas (`Area: ... .raw (F#)` columns) and the acquisition sequence,
   drops blanks and negative-batch (identification-only) runs, and log2-transforms.
2. **LOESS**: per-batch LOESS fit to the PQC injections vs injection order
   (span chosen per batch by 5-fold cross-validation, range 0.25-1.00).
3. **SERRF**: per-batch random-forest model (`ranger`, 200 trees) per feature, trained on
   the PQC injections using the 10 most correlated features as predictors
   (SERRF; Fan et al., 2019 - see [References](#references)).
4. Both corrections are followed by PQC-anchored between-batch median scaling.
5. Assesses each condition by feature-wise CV within each QC type and by PCA
   (points and 95% ellipses **coloured by batch**).

## Input format

| File | Sheet | Required columns |
|---|---|---|
| Compound Discoverer export (`.xlsx`) | `Compounds` | `Name`, `Formula`, `RT [min]`, `m/z`, `Calc. MW`, and one `Area: <file>.raw (F<id>)` column per injection |
| Sequence (`.xlsx`) | `InputFiles` | `Study File ID` (`F<id>`), `File Name`, `Creation Date` (`dd.mm.yyyy hh:mm:ss`), `Batch`, `Sample Type`, `Group` |

Sample types are derived as: `Sample Type == "Quality Control"` -> PQC; `Group == "LTQC"` -> LTQC
(Oslo-RM); `Group == "NIST"` -> NIST; everything else -> study sample. Rows with
`Sample Type == "Blank"` or a negative batch number are excluded.

## Running

```bash
Rscript normalization_comparison.R data/compounds.xlsx data/sequence.xlsx output
```

Missing packages (`readxl`, `tidyverse`, `ggplot2`, `ggsci`, `patchwork`, `cowplot`,
`ranger`, `data.table`, `ggrepel`, and Bioconductor `pcaMethods`) are installed on first run.
Developed with R 4.3.3.

## Outputs

`output/figures/*_for_publication.pdf|tiff` is the main figure (A: PCA, B: median log10
intensity vs injection order, C: CV per QC type, one column per condition).
`output/tables/` contains `CV_summary.csv`, `CV_per_feature.csv` and the corrected data
(wide and long) per condition.

## Try it with the synthetic example

`example_data/` holds a small, fully simulated dataset (150 features, 6 batches, with
batch offsets and within-batch drift) in the input format above; `make_example_data.R`
regenerates it with a fixed seed. It contains no real data.

```bash
Rscript normalization_comparison.R example_data/compounds_example.xlsx example_data/sequence_example.xlsx output
```

## Data

Raw study data are not included in this repository and are not publicly available.
Use your own Compound Discoverer export and sequence file in the formats described above.

## References

Fan S, Kind T, Cajka T, Hazen SL, Tang WHW, Kaddurah-Daouk R, Irvin MR, Arnett DK,
Barupal DK, Fiehn O. Systematic Error Removal Using Random Forest for Normalizing
Large-Scale Untargeted Lipidomics Data. *Anal. Chem.* 2019, 91 (5), 3590-3596.
https://doi.org/10.1021/acs.analchem.8b05592

## License

MIT - see [LICENSE](LICENSE).
