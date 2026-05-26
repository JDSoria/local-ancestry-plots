# Local ancestry plotting scripts

R scripts for plotting local ancestry segments across autosomes. The inputs are tab-separated files with one row per genomic segment and separate ancestry/color assignments for haplotype 1 and haplotype 2.

## Scripts

- `LocalAncestryMosaic.R`: plots the local ancestry mosaic for one dataset, either `LC` or `HC`.
- `LocalAncestryComparison.R`: plots a side-by-side comparison of `LC` and `HC` for the same individual.

The plots include chromosomes `chr1` to `chr22`, both haplotypes, and automatically filled grey regions for missing intervals or chromosome-end gaps.

## Requirements

Install R packages:

```r
install.packages(c("ggplot2", "dplyr"))
```

## Input data

Input files must be tab-separated `.tsv` files with these columns:

```text
Chr	Inicio	Fin	Num_SNPs	Ancestriahap1	Ancestriahap2	colorhap1	colorhap2
```

See [`DATA_FORMAT.md`](DATA_FORMAT.md) for the full input specification.

## Expected input naming

Files are expected to follow this naming pattern:

```text
InputG<individual_id>.tsv
```

For example:

```text
InputGHG01976.tsv
InputG100.tsv
```

By default, the scripts look for files inside dataset-specific directories:

```text
LC0.9_GAP/InputG<individual_id>.tsv
HC0.9_GAP/InputG<individual_id>.tsv
```

The exact repository layout is flexible, but the scripts need to know where those input folders are located. Both scripts accept the base directory as an optional argument.

## Single-dataset mosaic

Usage:

```bash
Rscript LocalAncestryMosaic.R <individual_id> <LC|HC> [base_directory]
```

Example:

```bash
Rscript LocalAncestryMosaic.R HG01976 LC /media/daniel/Espacio/NuevosBED
Rscript LocalAncestryMosaic.R HG01976 HC /media/daniel/Espacio/NuevosBED
```

Output examples:

```text
LocalAncestry_LC_HG01976.tiff
LocalAncestry_HC_HG01976.tiff
```

## LC vs HC comparison

Usage:

```bash
Rscript LocalAncestryComparison.R <individual_id> [base_directory]
```

Example:

```bash
Rscript LocalAncestryComparison.R HG01976 /media/daniel/Espacio/NuevosBED
```

Output example:

```text
LocalAncestry_Comparison_HG01976.tiff
```

## Colors

The scripts use the color values already present in the input files. The current plotting legend assumes:

| Ancestry | Color |
|---|---|
| Europe | `#1B4F72` |
| Native American indigenous | `#C62828` |
| Sub-Saharian Africa | `#F4D03F` |
| Unknown | `black` |
| GAP | `gray` |

## Notes for GitHub

Large generated figures and local output files are ignored by default through `.gitignore`. If you want to include selected example figures, add them explicitly with `git add -f <file>`.
