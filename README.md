# CODfinder pipeline

<p align="center"> <strong>Cancer Bioinformatics · Transcriptomics · Gene Co-expression Domains</strong> </p> <p align="center"> <img src="https://img.shields.io/badge/Research-Cancer%20Bioinformatics-6C8EBF?style=flat-square" /> <img src="https://img.shields.io/badge/Data-RNA--seq-4B8B8B?style=flat-square" /> <img src="https://img.shields.io/badge/Data-CPTAC-64748B?style=flat-square" /> <img src="https://img.shields.io/badge/Focus-Co--expression%20Domains-9B7EBD?style=flat-square" /> <img src="https://img.shields.io/badge/Cancer-Pancreatic%20Adenocarcinoma-B07A7A?style=flat-square" /> <img src="https://img.shields.io/badge/Language-R-8A9BA8?style=flat-square" /> </p>

An R script that finds **COD (co-expression domains)**: stretches of neighbouring genes on a chromosome whose expression is correlated across samples. It works on any bulk gene-expression cohort (any tumor type, tissue, or condition, and any expression platform), as long as the input file matches the format below.

The COD-search algorithm is based on Jose Guerrero's CODfinder: https://github.com/joseguem/CODfinder

## Requirements

- R 4.x
- Bioconductor packages `GenomicFeatures` and a `TxDb` annotation package. They are installed automatically on first run. The default is `TxDb.Hsapiens.UCSC.hg38.knownGene` (human, hg38).

## Usage

```bash
Rscript codfinder.R <input_file> <output_dir> [gene_of_interest]

# example
Rscript codfinder.R data/expression.txt results SCAP
```

## Input file

A tab-delimited text file, **one row per gene, one column per sample**:

| Gene ID column | Symbol column | Sample 1 | Sample 2 | ... |
|---|---|---|---|---|
| 7157 | TP53 | 120.5 | 98.0 | ... |

- Two annotation columns plus any number of numeric sample columns. No other non-numeric columns are allowed.
- Default column names are `Entrez_Gene_Id` and `Hugo_Symbol` (the cBioPortal layout). Other names can be set with `id_column` and `symbol_column` at the top of the script.
- Gene IDs must be **Entrez IDs**, because coordinates are looked up in the TxDb package by Entrez ID.
- Values are expression levels such as RSEM, TPM, FPKM or normalized counts. By default they are transformed with `log2(x + 1)`; set `log_transform <- FALSE` if the data is already on a log scale.
- Tumor and normal samples should be analysed in separate runs.

## Logic

1. **Load** the expression table and drop rows without a gene ID.
2. **Annotate** each gene with chromosome, start, end and strand from the TxDb package. Genes without coordinates are dropped.
3. **Sort** genes by chromosome and genomic position. The rank of a gene along its chromosome (gene index) is the unit of distance used everywhere below.
4. **Filter**: `log2(x + 1)`, then remove genes whose mean expression is at or below the `expression_quantile` quantile (default 0.30) and genes with zero variance. Removed genes are listed in `removed_genes.tsv`.
5. **Correlate**: for each chromosome, a gene x gene Pearson correlation matrix across samples.
6. **Scan window sizes** (`window_sizes`, default 3-20). For each window size `w`:
   - *Binsignal* of gene `i`: mean correlation in the block between the `w` genes ending at `i` and the `w` genes following it. High values mean that the genes on both sides of `i` are co-expressed.
   - A t-test compares the `w` binsignal values on the left of `i` with the `w` on the right, giving a p-value for a boundary at `i`.
   - Threshold: the mean binsignal over chr1-chr22 and chrX.
   - A COD starts where binsignal exceeds the threshold (at least 3 of the next 4 genes) with a significant boundary, extends while the signal stays above the threshold, and ends near a significant boundary. CODs must contain at least 4 genes.
   - Each COD gets its intra-COD co-expression: the mean r over all gene pairs inside it.
7. **Pick the best window**: the one with the highest mean intra-COD co-expression on chr1. The script stops if chr1 has no CODs for any window.
8. **Final tables and figures** for the best window: per-chromosome heatmaps with binsignal, per-chromosome COD counts and densities, and distributions of COD size, length in kb and intra-COD co-expression.
9. **Gene of interest** (optional): a zoomed triangular heatmap of its COD, or of its neighbourhood if it is in no COD.
10. **COD x COD correlation**: mean r between all gene pairs of two CODs, within and between chromosomes. The diagonal is the intra-COD value.
11. **Gene-of-interest analysis** (optional): the 100 genes most correlated with it, plus the 10 CODs on the same chromosome and the 10 on other chromosomes most correlated with its COD. Significance comes from a permutation test (`n_permutations`, default 1000) that shuffles the COD's expression profile across samples and compares each observed r with the permuted maximum over all candidate CODs (max-T correction). It is exploratory: selection and testing use the same data.
12. **Mirrored control**: CODs are mirrored along each chromosome (same sizes, different positions) and the same statistics are recomputed. It is a single control, not a random ensemble, so no p-values are given.

## Output (`<output_dir>/`)

```
matrix/
  <chr>.correlation.tsv      gene x gene correlation matrix
  <chr>.annotation.tsv       gene order and coordinates
  <chr>.binsignal.tsv        binsignal and p-value per gene (best window)
  <chr>.COD.tsv              CODs on the chromosome
  all_COD.tsv                all CODs
  COD_stats_by_chromosome.tsv
  mean_binsignal_*.tsv
  removed_genes.tsv
  COD_x_COD/                 COD annotation, full and per-chromosome matrices,
                             pairwise tables (within / between chromosomes)
  <GENE>_analysis/           tables from step 11
  mirrored_control/          tables from step 12
window_test/
  all_COD_all_windows.tsv    every COD for every window
  window_summary.tsv         per-window statistics (chr1 and genome)
  thresholds.tsv, COD_count_by_chromosome.tsv
  window_selection_chr1.svg, window_selection_genome.svg
figures/
  expression_before_filtering.svg
  gene_pair_correlation_distribution.svg
  heatmap_<chr>.svg
  COD_by_chromosome.svg, COD_distributions.svg
  <GENE>_COD_zoom.svg, <GENE>_COD_genes.tsv
  COD_x_COD/*.svg
  mirrored_control/*.svg
```

Key COD columns: `START` / `END` (gene indices on the chromosome), `n_genes`, `first_gene_name`, `last_gene_name`, `start_bp`, `end_bp`, `size_kb`, `intra_coexpression`.

## Main parameters (top of the script)

| Parameter | Default | Meaning |
|---|---|---|
| `window_sizes` | `3:20` | Window sizes scanned |
| `expression_quantile` | `0.30` | Genes at or below this mean-expression quantile are removed |
| `log_transform` | `TRUE` | Apply `log2(x + 1)` |
| `id_column`, `symbol_column` | `Entrez_Gene_Id`, `Hugo_Symbol` | Input column names |
| `txdb_package` | hg38 human | Annotation package |
