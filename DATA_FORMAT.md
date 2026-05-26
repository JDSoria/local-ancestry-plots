# Input data format

The scripts expect tab-separated `.tsv` files with one row per local ancestry segment.

## File name

Each input file must be named as:

```text
InputG<individual_id>.tsv
```

Examples:

```text
InputG100.tsv
InputGHG01976.tsv
InputGBEP18.tsv
```

The individual ID passed to the script must match the part after `InputG` and before `.tsv`.

For example:

```bash
Rscript LocalAncestryMosaic.R HG01976 LC /media/daniel/Espacio/NuevosBED
```

looks for:

```text
/media/daniel/Espacio/NuevosBED/LC0.9_GAP/InputGHG01976.tsv
```

## Separator

Files must be tab-separated, not comma-separated.

## Required columns

Expected header:

```text
Chr	Inicio	Fin	Num_SNPs	Ancestriahap1	Ancestriahap2	colorhap1	colorhap2
```

| Column | Type | Description |
|---|---|---|
| `Chr` | string | Chromosome name. Expected values are `chr1`, `chr2`, ..., `chr22`. |
| `Inicio` | integer | Segment start position in base pairs. |
| `Fin` | integer | Segment end position in base pairs. |
| `Num_SNPs` | integer | Number of SNPs in the segment. This column is not used for plotting. |
| `Ancestriahap1` | string | Local ancestry assigned to haplotype 1. |
| `Ancestriahap2` | string | Local ancestry assigned to haplotype 2. |
| `colorhap1` | string | Plot color for haplotype 1. Hex colors and R-recognized color names are accepted. |
| `colorhap2` | string | Plot color for haplotype 2. Hex colors and R-recognized color names are accepted. |

## Minimal example

```text
Chr	Inicio	Fin	Num_SNPs	Ancestriahap1	Ancestriahap2	colorhap1	colorhap2
chr1	6786703	9351317	0	European	unknown	#1B4F72	black
chr1	9351317	9840876	0	unknown	unknown	black	black
chr1	9840876	12121615	0	unknown	European	black	#1B4F72
chr1	12121615	12759228	0	Native American	European	#C62828	#1B4F72
```

## Coordinate and interval rules

- Chromosomes must be named `chr1` to `chr22`.
- Positions are assumed to be hg38 coordinates because the scripts use hg38 chromosome sizes to fill missing regions.
- Each segment must satisfy `Inicio < Fin`.
- Segments should be ordered or orderable by `Inicio` within each chromosome.
- Missing intervals between consecutive segments are plotted automatically as grey `GAP` regions.
- Missing chromosome starts and ends are also filled as grey `GAP` regions.
- The plotting colors are taken directly from `colorhap1` and `colorhap2`.

## Dataset folders

For the comparison plot, both files are needed for the same individual:

```text
LC0.9_GAP/InputG<individual_id>.tsv
HC0.9_GAP/InputG<individual_id>.tsv
```

For the single-dataset mosaic, only the selected dataset is required:

```text
LC0.9_GAP/InputG<individual_id>.tsv
```

or:

```text
HC0.9_GAP/InputG<individual_id>.tsv
```
