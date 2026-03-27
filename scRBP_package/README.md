<table>
  <tr>
    <td>
      <img src="https://raw.githubusercontent.com/mayunlong89/scRBP/main/Examples/scRBP_logo.png" width="300">
    </td>
    <td>
      <h1>scRBP</h1>
      <p>A scalable framework for inferring RNA-binding protein regulons from single-cell data</p>
    </td>
  </tr>
</table>

![pypi](https://img.shields.io/badge/pypi-0.1.3.4-green)
![python](https://img.shields.io/badge/python-3.9--3.11-blue)
![license](https://img.shields.io/badge/license-MIT-yellow)

**scRBP (Single-cell RNA Binding Protein Regulon Inference)** is a command-line toolkit for comprehensive analysis of RNA-binding proteins (RBPs) in single-cell RNA-seq data. scRBP provides a systematic, scalable and integrative framework to infer RBP-mediated gene and isoform regulatory networks ("regulons") from single-cell transcriptomes and prioritize networks underlying complex genetic traits and disorders. scRBP is comprised of six main modules: (i) developing a comprehensive compendium of RBPs and their associated motif clusters from diverse public resources; (ii) systematic, motif-guided transcriptome-wide inference of RBP targets at both gene- and isoform-level resolution; (iii) construction of RBP-gene and/or RBP-isoform co-expression networks from short- or long-read single-cell transcriptomic data, respectively; (iv) defining high-fidelity regulons by integrating RBP-target interactions, and quantifying cell type-specific regulon activity scores (RAS); (v) integrating GWAS results to compute regulon-level genetic association scores (RGS); and (vi) constructing a unified trait-relevance score (TRS) by combining RAS and RGS for each regulon in a given cellular context, with statistical significance assessed using Monte Carlo (MC) sampling.

---

## What scRBP Does

RBPs are key post-transcriptional regulators that control mRNA splicing, stability, and translation. scRBP enables you to:

- **Construct** which RBPs regulate which genes or isoforms in your single-cell data
- **Prune** raw RBP–gene associations using motif-binding evidence to obtain high-confidence regulons
- **Score** each cell or cell type for regulon activity score (RAS) using the AUCell algorithm
- **Link** RBP regulons to human disease through GWAS genetic enrichment (RGS via MAGMA)
- **Integrate** RAS and RGS into a unified Trait Relevance Score (TRS) that ranks disease-relevant RBPs

---

## Pipeline at a Glance

```
Raw single-cell data (.h5ad / .feather)
          │
          ▼
[Step 1]  scRBP getSketch        ── Stratified GeoSketch cell downsampling
          │
          ▼
[Step 2]  scRBP getGRN           ── GRNBoost2/GENIE3 RBP→Gene/Isoform inference
          │                          (run N seeds for robustness, default 30 times)
          ▼
[Step 3]  scRBP getMerge_GRN     ── Merge N-seed GRNs → consensus network
          │
          ▼
[Step 4]  scRBP getModule        ── Extract regulon candidates (Top-N / percentile)
          │
          ▼
[Step 5]  scRBP getPrune         ── Motif-enrichment pruning via ctxcore
          │
          ▼
[Step 6]  scRBP getRegulon       ── Export pruned regulons to GMT format
          │
          ▼
[Step 7]  scRBP mergeRegulons    ── Merge region-specific GMT files
          │                          (3'UTR / 5'UTR / CDS / Introns)
          ▼
[Step 8]  scRBP ras              ── Regulon Activity Score (AUCell) per cell / cell type
          │                          (--mode sc | --mode ct)
          ▼
[Step 9]  scRBP rgs              ── Regulon Gene-Set analysis (MAGMA GWAS enrichment)
          │                          (--mode sc | --mode ct)
          ▼
[Step 10] scRBP trs              ── Trait Relevance Score (RAS × RGS integration)
                                     (--mode sc | --mode ct)
```

---

## Installation

### Requirements

- Python **3.9, 3.10, or 3.11** (Python 3.12+ not yet supported by `pyscenic`/`arboreto`)
- MAGMA binary (external, required only for Step 9 — `scRBP rgs`)

### Option 1 — Install from PyPI (recommended)

```bash
pip install scRBP
```

This installs scRBP together with all Python dependencies in one step.

### Option 2 — Install from source (development)

```bash
git clone https://github.com/mayunlong89/scRBP.git
cd scRBP/scRBP_package
pip install -e .
```

### Option 3 — Install via conda (recommended for HPC / cluster)

```bash
git clone https://github.com/mayunlong89/scRBP.git
cd scRBP/scRBP_package

conda env create -f environment.yml
conda activate scrbp

pip install -e .
```

### Install MAGMA (for Step 9 only)

MAGMA is a standalone binary not available on PyPI. Download from https://cncr.nl/research/magma and make it executable:

```bash
# Install MAGMA (v1.10, Linux static)

# 1. Create installation directory
mkdir -p ~/tools/magma
cd ~/tools/magma

# 2. Download MAGMA (note: must include /download)
wget -O magma_v1.10_static_linux.zip \
"https://vu.data.surf.nl/index.php/s/lxDgt2dNdNr6DYt/download"

# 3. Unzip the package
unzip magma_v1.10_static_linux.zip

# 4. Check extracted files
ls

# 5. Enter the extracted directory (name may vary)
cd magma*

# 6. Make the binary executable
chmod +x magma

# 7. Verify installation
./magma --version

# 8. Optionally add MAGMA to PATH
echo 'export PATH=~/tools/magma/magma_v1.10_static:$PATH' >> ~/.bashrc
source ~/.bashrc

```

### Verify installation

```bash
scRBP --help
scRBP getGRN --help
scRBP --version
```

---

## Step 1 — `scRBP getSketch`  (Optional)

### Purpose

Performs **stratified GeoSketch downsampling** on large single-cell datasets. For `.h5ad` input, cells are sampled per cell-type so that each type contributes at least a user-specified minimum number of cells, while keeping the total close to a global target. For `.feather` input, a global GeoSketch is applied.

GeoSketch uses a PCA-space geometric sketch to retain transcriptional diversity while making downstream analysis faster and statistically representative.

> For datasets exceeding 300,000 cells, direct analysis may incur substantial computational costs. We recommend downsampling to ~50,000 cells.

### Usage

```bash
scRBP getSketch \
    --input  <input_file>   \
    --output <output_file>  \
    [options]
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--input` | str | — | **Yes** | Input expression file. Accepted formats: `.h5ad` (AnnData) or `.feather` (gene × cell matrix). `.rds` (Seurat) must be converted to `.h5ad` first using `SeuratDisk` in R. |
| `--output` | str | — | **Yes** | Output file path. Format auto-detected from extension: `.h5ad`, `.csv`, `.feather`, `.npz`. |
| `--n_cells` | int | 50000 | No | Target total number of cells to retain after sketching. |
| `--n_pca` | int | 100 | No | Number of PCA (TruncatedSVD) components before running GeoSketch. |
| `--celltype_col` | str | `"celltype"` | No | Column name in `adata.obs` for cell-type labels. Used for stratified sampling (`.h5ad` only). |
| `--min_cells_per_type` | int | 50 | No | Minimum number of cells to retain per cell type. If a type has fewer cells, all are kept. |
| `--seed` | int | 42 | No | Random seed for reproducibility. |

### Output Files

| File | Description |
|------|-------------|
| `<output>` | Downsampled expression matrix in specified format. |
| `<output_prefix>_cell_to_celltype.csv` | *(`.h5ad` input only)* Two-column CSV mapping each retained cell barcode to its cell type. |

### Example

```bash
scRBP getSketch \
    --input  PBMC_healthy.h5ad \
    --output PBMC_sketch_50K.feather \
    --n_cells 50000 \
    --celltype_col annotation_broad3 \
    --min_cells_per_type 500 \
    --n_pca 100 \
    --seed 42
```

---

## Step 2 — `scRBP getGRN`

### Purpose

Infers a **gene or isoform regulatory network (GRN)** between RBPs and their target genes/isoforms using tree-based machine learning — either **GRNBoost2** (gradient boosting, faster) or **GENIE3** (random forest, more conservative). Optionally computes Spearman correlation and assigns a regulatory mode (`activating` / `repressing` / `unknown`) for each edge.

Two inference modes:
- **`gene` mode**: RBPs in the expression matrix are used directly as regulators; all genes are targets.
- **`isoform` mode**: Multiple isoforms of the same RBP gene are first aggregated into a single regulator signal, then individual isoforms serve as targets. Requires `--isoform_annotation`.

> **Recommended practice:** run `getGRN` with 30 different `--seed` values, then merge with `getMerge_GRN` to obtain a stable consensus network.

### Usage

```bash
scRBP getGRN \
    --matrix   <expression_matrix> \
    --rbp_list <rbp_list_file>     \
    --output   <output_prefix>     \
    --mode     gene|isoform        \
    [options]
```

### Output filename auto-suffix

| Mode | Output filename |
|------|----------------|
| `gene` | `<prefix>_scRBP_gene_GRNs.tsv` |
| `isoform` | `<prefix>_scRBP_isoform_GRNs.tsv` |

### Parameters — Core (all modes)

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--matrix` | str | — | **Yes** | Expression matrix file. Formats: `.csv`, `.csv.gz`, `.feather`, `.loom`. Rows = features, columns = cells. |
| `--rbp_list` | str | — | **Yes** | Plain-text file of RBP gene symbols (one per line). |
| `--output` | str | — | **Yes** | Output prefix. Mode suffix appended automatically. |
| `--mode` | str | `gene` | No | `gene` — RBP→Gene network; `isoform` — RBP→Isoform network. |
| `--method` | str | `grnboost2` | No | `grnboost2` (recommended) or `genie3`. |
| `--n_workers` | int | all CPUs | No | Number of parallel worker processes. |
| `--correlation` | bool | `True` | No | Compute Spearman correlation and assign regulatory mode per edge. |
| `--threshold` | float | 0.03 | No | Absolute Spearman correlation threshold; edges with `|r| ≤ threshold` are removed. |
| `--seed` | int | 1234 | No | Random seed for the tree ensemble. Change across runs for consensus merging. |

### Parameters — Isoform mode only

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--isoform_annotation` | str | `None` | **Yes** (isoform) | TSV/CSV mapping isoform/transcript IDs → gene symbols. |
| `--rbp_agg_method` | str | `sum` | No | How to aggregate multiple isoforms of the same RBP: `sum`, `mean`, `max`. |
| `--remove_self_targets` | bool | `True` | No | Remove edges where the target isoform belongs to the same gene as the RBP. |
| `--min_target_cells_expressed` | int | 10 | No | Minimum number of cells with expression > 0 to keep a target isoform. |
| `--min_target_mean_expr` | float | 0.01 | No | Minimum mean expression across all cells to keep a target isoform. |

### Output Files

| File | Description |
|------|-------------|
| `<prefix>_scRBP_gene_GRNs.tsv` | Gene mode GRN. Columns: `RBP`, `Gene`, `Importance`, `Correlation`, `Mode` |
| `<prefix>_scRBP_isoform_GRNs.tsv` | Isoform mode GRN. Columns: `RBP`, `Isoform`, `Importance`, `Correlation`, `Mode` |
| `<prefix>_rbp_aggregated_expr.tsv` | Aggregated RBP expression matrix (isoform mode, if `--save_rbp_agg_matrix True`) |
| `<prefix>_rbp_isoform_map.tsv` | RBP → isoform membership table (isoform mode) |
| `<prefix>_target_isoform_stats.tsv` | Per-isoform expression statistics used for filtering (isoform mode) |

### Examples

**Gene mode (30 seeds):**

```bash
for SEED in $(seq 1 30); do
  scRBP getGRN \
      --matrix    PBMC_sketch_15K.feather \
      --rbp_list  human_RBP_list.txt \
      --output    grn_seed${SEED} \
      --mode      gene \
      --method    grnboost2 \
      --n_workers 20 \
      --correlation True \
      --seed      ${SEED}
done
# Output: grn_seed1_scRBP_gene_GRNs.tsv, grn_seed2_scRBP_gene_GRNs.tsv, ...
```

**Isoform mode:**

```bash
for SEED in 1 2 3; do
  scRBP getGRN \
      --matrix                     PBMC_isoform_sketch.feather \
      --rbp_list                   human_RBP_list.txt \
      --output                     iso_grn_seed${SEED} \
      --mode                       isoform \
      --isoform_annotation         gencode_v44_isoform_gene_map.tsv \
      --rbp_agg_method             sum \
      --remove_self_targets        True \
      --min_target_cells_expressed 10 \
      --min_target_mean_expr       0.01 \
      --method                     grnboost2 \
      --n_workers                  20 \
      --seed                       ${SEED}
done
```

---

## Step 3 — `scRBP getMerge_GRN`

### Purpose

Merges multiple GRN output files (from different random seeds) into a single **consensus network**. For each RBP–Gene edge: importance and correlation are averaged across all runs, `n_present` counts how many runs contained the edge, and `presence_rate` = n_present / total_runs indicates edge stability. Optional filters retain only high-confidence edges.

### Usage

```bash
scRBP getMerge_GRN \
    --pattern "<glob_pattern>" \
    --output  <merged_output>  \
    [options]
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--pattern` | str | — | **Yes** | Glob pattern matching all GRN seed files, e.g. `"grn_seed*.tsv"`. Must be quoted. |
| `--output` | str | — | **Yes** | Output path for the merged TSV file. |
| `--corr-threshold` | float | 0.0 | No | Remove edges where averaged `|Correlation|` ≤ this value. |
| `--n_present` | int | 0 | No | Minimum number of seed runs an edge must appear in to be retained. |
| `--present_rate` | float | 0.0 | No | Minimum presence rate (n_present / N_runs). Typical cutoff: `0.5`. |

### Example

```bash
scRBP getMerge_GRN \
    --pattern "grn_seed*.tsv" \
    --output  grn_consensus.tsv \
    --n_present 15 \
    --present_rate 0.5
```

Edges appearing in fewer than 50% of seeds are discarded, yielding a stable consensus network.

---

## Step 4 — `scRBP getModule`

### Purpose

Extracts **RBP regulon candidate modules** from the consensus GRN using multiple complementary strategies:

- **TopN-per-gene** (`topN_per_gene_N`): For each gene, keep only its top-N highest-importance RBPs.
- **TopN-per-RBP** (`top_target_N`): For each RBP, keep its top-N highest-importance target genes.
- **Percentile** (`pctX_per_rbp`): For each RBP, keep targets above the Xth importance percentile.

### Usage

```bash
scRBP getModule \
    --input         <consensus_grn.tsv> \
    --output_merged <modules_output.tsv> \
    [options]
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--input` | str | — | **Yes** | Consensus GRN TSV file. Required columns: `RBP`, `Gene`, `Importance`. |
| `--output_merged` | str | — | **Yes** | Output TSV file containing all strategies merged. |
| `--importance_threshold` | float | 0.005 | No | Global minimum importance cutoff applied before any strategy. |
| `--top_n_list` | str | `"5,10,50"` | No | Comma-separated N values for TopN-per-gene strategy. |
| `--target_top_n` | str | `"50"` | No | Comma-separated N values for TopN-per-RBP strategy. |
| `--percentile` | str | `"0.75,0.9"` | No | Comma-separated percentile thresholds (0–1) for percentile-per-RBP strategy. |
| `--verbose` | flag | off | No | Print per-strategy edge counts to stdout. |

### Example

```bash
scRBP getModule \
    --input              grn_consensus.tsv \
    --output_merged      modules.tsv \
    --importance_threshold 0.005 \
    --top_n_list         "5,10,50" \
    --target_top_n       "50" \
    --percentile         "0.75,0.9" \
    --verbose
```

---

## Step 5 — `scRBP getPrune`

### Purpose

**Prunes candidate RBP regulons** by testing whether the potential target genes of each RBP are statistically enriched for known RBP-binding motifs using the `ctxcore` enrichment scoring engine. This step filters out spurious GRN edges that lack motif support and produces pruned regulons with enrichment statistics (AUC, NES, leading-edge genes).

For each RBP × Strategy combination:
1. Maps the RBP to its associated sequence motifs via `--motif_rbp_links`.
2. Scores the ranked list of target genes against the motif ranking matrix.
3. Retains RBPs passing all threshold filters (rank, AUC, NES).

### Usage

```bash
scRBP getPrune \
    --rbp_targets        <modules.tsv>         \
    --motif_rbp_links    <motif2rbp.csv>        \
    --motif_target_ranks <rankings.feather>    \
    --save_dir           <output_directory/>    \
    [options]
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--rbp_targets` | str | — | **Yes** | TSV/CSV of RBP regulon candidates from `getModule`. |
| `--motif_rbp_links` | str | — | **Yes** | CSV/TSV linking motif IDs to RBP names (from pySCENIC / cistarget databases). |
| `--motif_target_ranks` | str | — | **Yes** | Feather file containing the motif × gene ranking matrix (cistarget database). |
| `--save_dir` | str | — | **Yes** | Output directory. Created if it does not exist. |
| `--rank_threshold` | int | 1500 | No | Maximum rank for a gene to be considered in enrichment scoring. |
| `--auc_threshold` | float | 0.05 | No | AUC threshold; RBPs with AUC < this value are discarded. |
| `--nes_threshold` | float | 3.0 | No | NES threshold; RBPs with NES < this value are discarded. |
| `--min_genes` | int | 20 | No | Minimum number of target genes required for an RBP regulon to be scored. |
| `--n_jobs` | int | all CPUs | No | Number of parallel worker processes. |
| `--weight_mode` | str | `importance` | No | Gene weighting mode. `importance` uses GRN importance scores; `equal` weights all genes equally. |

### Example

```bash
scRBP getPrune \
    --rbp_targets        modules.tsv \
    --motif_rbp_links    motif_to_rbp_hg38_v10.csv \
    --motif_target_ranks hg38_3UTR_rankings.feather \
    --save_dir           pruned_results/ \
    --rank_threshold     1500 \
    --auc_threshold      0.05 \
    --nes_threshold      3.0 \
    --min_genes          20 \
    --n_jobs             16
```

---

## Step 6 — `scRBP getRegulon`

### Purpose

Converts pruned enrichment results (from `getPrune`) into **GMT regulon files** in two formats:
- **Symbol GMT**: Gene names as HGNC symbols (for use with `ras`).
- **Entrez GMT**: Gene names as NCBI Entrez IDs (for use with `rgs` / MAGMA).

### Usage

```bash
scRBP getRegulon \
    --input      <ctx_scores.csv>        \
    --out-symbol <regulons_symbol.gmt>   \
    --out-entrez <regulons_entrez.gmt>   \
    [options]
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--input` | str | — | **Yes** | Input CSV file with pruned enrichment scores from `getPrune`. |
| `--out-symbol` | str | — | **Yes** | Output path for the Symbol-format GMT file. |
| `--out-entrez` | str | — | **Yes** | Output path for the Entrez-format GMT file. |
| `--min_genes` | int | 1 | No | Minimum number of target genes required to include an RBP in the output GMT. |
| `--map-custom` | str | `None` | No | Custom gene-to-Entrez mapping file (e.g. NCBI `Homo_sapiens.gene_info` or a `gene.loc` file). |
| `--map-hgnc` | str | `None` | No | HGNC complete set file for Symbol → Entrez mapping. |
| `--map-ncbi` | str | `None` | No | NCBI gene2refseq or gene_info file for Symbol → Entrez mapping. |
| `--taxid` | int | 9606 | No | NCBI taxonomy ID for filtering when using `--map-ncbi`. Default `9606` = human. |
| `--drop-unmapped-genes` | flag | off | No | Remove genes that could not be mapped to an Entrez ID. |
| `--drop-empty-sets` | flag | off | No | Remove RBP regulons that become empty after gene mapping. |

### Example

```bash
scRBP getRegulon \
    --input       pruned_results/ctx_scores_merged.csv \
    --out-symbol  regulons_symbol.gmt \
    --out-entrez  regulons_entrez.gmt \
    --map-custom  NCBI38.gene.loc \
    --min_genes   5 \
    --drop-unmapped-genes \
    --drop-empty-sets
```

---

## Step 7 — `scRBP mergeRegulons`

### Purpose

Merges **region-specific GMT regulon files** (e.g., 3′UTR, 5′UTR, CDS, Introns) generated from separate `getPrune` + `getRegulon` runs into a single combined GMT file. The function traverses a directory tree, finds GMT files matching a specified name pattern within region subdirectories, and concatenates them in a defined region priority order.

### Usage

```bash
scRBP mergeRegulons \
    --base_dir <base_directory/>  \
    --input    <gmt_filename>     \
    --output   <merged_output.gmt>\
    [options]
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--base_dir` | str | — | **Yes** | Root directory to search for region subdirectories. |
| `--input` | str | — | **Yes** | GMT filename to look for inside each region subdirectory (e.g. `regulons_symbol.gmt`). |
| `--output` | str | — | **Yes** | Filename for the merged GMT output. |
| `--recursive` | flag | off | No | Process multiple tissue/dataset directories under `--base_dir`. |
| `--tissue_glob` | str | `"z_GRNBoost2_*_30times"` | No | Glob pattern for tissue/dataset directories when `--recursive` is used. |
| `--region_glob` | str | `"Results_final_*_RBP_top1500_*"` | No | Glob pattern for region subdirectories. |
| `--append_region_to_setname` | flag | off | No | Append region tag (e.g. `_3UTR`) to each regulon set name. |
| `--dedup_lines` | flag | off | No | Remove exact duplicate GMT lines after merging. |
| `--region_order` | list | `["3UTR","5UTR","CDS","Introns"]` | No | Priority order for regions. |
| `--overwrite` | flag | off | No | Overwrite the output file if it already exists. |
| `--summary_out` | str | `None` | No | Optional path to write a TSV summary table of per-region regulon counts. |

### Example

```bash
# Simple (single directory)
scRBP mergeRegulons \
    --base_dir ./analysis/ \
    --input    regulons_symbol.gmt \
    --output   regulons_combined.gmt

# Recursive (multiple tissues / datasets)
scRBP mergeRegulons \
    --base_dir /data/scRBP_results/ \
    --input    regulons_symbol.gmt \
    --output   regulons_all_regions.gmt \
    --recursive \
    --tissue_glob "z_GRNBoost2_*_30times" \
    --region_glob "Results_final_*_RBP_top1500_*" \
    --append_region_to_setname \
    --region_order 3UTR 5UTR CDS Introns \
    --summary_out region_summary.tsv
```

---

## Step 8 — `scRBP ras`

### Purpose

Computes **Regulon Activity Scores (RAS)** using the **AUCell** algorithm — for each cell or cell type, AUCell scores how active each RBP regulon is based on whether its target genes are at the top of that cell's expression ranking. Optionally computes **Regulon Specificity Scores (RSS)** based on Jensen-Shannon divergence.

**Two modes are supported:**
- **`--mode sc`**: RAS computed at single-cell resolution.
- **`--mode ct`**: RAS aggregated to cell-type level; also computes RSS and requires `--celltypes-csv`.

### Usage

```bash
scRBP ras \
    --mode     <sc|ct>             \
    --matrix   <expression_file>   \
    --regulons <regulons.gmt>      \
    --out      <output_prefix/>    \
    [options]
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--mode` | str | `ct` | No | Analysis mode. `sc` = single-cell AUCell scores; `ct` = cell-type-aggregated scores + RSS. |
| `--matrix` | str | — | Cond. | Expression matrix file. Formats: `.csv`, `.csv.gz`, `.feather`, `.loom`. Rows = genes, columns = cells. Required unless `--aucell-in` is provided. |
| `--regulons` | str | — | Cond. | Regulon file. Formats: `.gmt` / `.gmt.gz` (Symbol) or `.pkl` (pySCENIC). Required unless `--aucell-in` is provided. |
| `--aucell-in` | str | — | No | Precomputed AUCell CSV (cells × regulons). If provided, skips AUCell computation. |
| `--out` | str | — | **Yes** | Output path prefix or directory. |
| `--out_format` | str | `csv` | No | Output format for AUCell scores. Choices: `csv`, `loom`, `both`. |
| `--n_workers` | int | 4 | No | Number of parallel workers for AUCell computation. |
| `--celltypes-csv` | str | — | Cond. | CSV file with cell barcode and cell type columns. **Required when `--mode ct`.** |
| `--cell-col` | str | auto | No | Column name for cell barcodes in `--celltypes-csv`. Auto-detected if not specified. |
| `--ctype-col` | str | auto | No | Column name for cell-type labels in `--celltypes-csv`. Auto-detected if not specified. |
| `--min_genes` | int | 1 | No | Drop regulons with fewer than this many genes in the expression matrix. |
| `--emit-expr-stats` | flag | on | No | Compute and save per-gene expression statistics (`mean_expr`, `pct_detected`). Used as input to `rgs`. |
| `--expr-stats-out` | str | auto | No | Output path for expression statistics TSV. Auto-generated from `--out` if not specified. |

### Output Files

| File | Description |
|------|-------------|
| `aucell_sc.csv` | AUCell scores — cells × regulons (`--mode sc`) |
| `aucell_ct.csv` | AUCell scores — cell-types × regulons (`--mode ct`) |
| `rss.tsv` | Regulon Specificity Scores per cell type (`--mode ct` only) |
| `expr_stats.tsv` | Per-gene expression statistics: `mean_expr`, `pct_detected` |

### Examples

**Single-cell mode (`--mode sc`):**

```bash
scRBP ras \
    --mode      sc \
    --matrix    PBMC_sketch_15K.feather \
    --regulons  regulons_symbol.gmt \
    --out       ras_sc_output/ \
    --n_workers 8 \
    --emit-expr-stats
```

**Cell-type mode (`--mode ct`):**

```bash
scRBP ras \
    --mode          ct \
    --matrix        PBMC_sketch_15K.feather \
    --regulons      regulons_symbol.gmt \
    --out           ras_ct_output/ \
    --celltypes-csv cell_to_celltype.csv \
    --n_workers     8 \
    --emit-expr-stats
```

---

## Step 9 — `scRBP rgs`

### Purpose

Computes **Regulon Gene-Set (RGS) scores** by running **MAGMA gene-set analysis** on each RBP regulon using GWAS summary statistics. This step assesses whether the target genes of each RBP regulon are enriched for GWAS disease-associated genes.

A key innovation is the automatic generation of **matched null regulons** — for each real regulon, `n-null` null gene sets are sampled from the genome, matched on four confounding dimensions: number of SNPs, number of parameters, mean gene expression, and percent detected. This controls for size and expression-level biases, enabling fair statistical comparison.

**Two modes are supported:**
- **`--mode sc`**: Per-regulon MAGMA analysis.
- **`--mode ct`**: Cell-type-stratified analysis; requires precomputed expression stats from `ras --emit-expr-stats`.

### Usage

```bash
scRBP rgs \
    --mode      <sc|ct>           \
    --magma     <magma_binary>    \
    --genes-raw <gwas.genes.raw>  \
    --sets      <regulons.gmt>    \
    --out       <output_prefix>   \
    [options]
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--mode` | str | — | **Yes** | Analysis mode. `sc` = per-regulon analysis; `ct` = cell-type-stratified analysis. |
| `--magma` | str | — | **Yes** | Full path to the MAGMA binary executable. |
| `--genes-raw` | str | — | **Yes** | MAGMA `.genes.raw` file from MAGMA gene analysis step. |
| `--sets` | str | — | **Yes** | GMT regulon file. Use Symbol GMT (`--id-type symbol`) or Entrez GMT (`--id-type entrez`, default). |
| `--id-type` | str | `entrez` | No | Gene ID type in GMT file: `entrez` or `symbol`. |
| `--out` | str | — | **Yes** | Output file prefix. |
| `--gene-loc` | str | — | No | MAGMA NCBI gene location file. Required when `--id-type symbol`. |
| `--n-null` | int | 1000 | No | Number of matched null regulons per real regulon. Higher = more stable empirical p-values. |
| `--seed` | int | 2025 | No | Random seed for null regulon sampling. |
| `--q-bins` | int | 10 | No | Number of quantile bins for null matching dimensions. |
| `--min_genes` | int | 0 | No | Minimum number of genes a regulon must have (after overlap with MAGMA universe) to be tested. |
| `--threads` | int | `None` | No | Number of threads for MAGMA. |
| `--expr-stats` | str | — | No | Precomputed expression statistics TSV from `ras --emit-expr-stats`. Used for 4D null matching. |

### Output Files

| File | Description |
|------|-------------|
| `<out>_real.csv` | Parsed MAGMA results for real regulons: `VARIABLE`, `NGENES`, `BETA`, `BETA_STD`, `SE`, `P`, `z`, `mlog10p` |
| `<out>_null_summary.tsv` | Empirical null distribution summary per regulon |
| `<out>.gsa.out` | Raw MAGMA gene-set analysis output |

### Examples

**Single-cell mode (`--mode sc`):**

```bash
scRBP rgs \
    --mode      sc \
    --magma     /tools/magma \
    --genes-raw scz_gwas.genes.raw \
    --sets      regulons_entrez.gmt \
    --id-type   entrez \
    --out       rgs_output/scz_sc_rgs \
    --n-null    1000
```

**Cell-type mode (`--mode ct`):**

```bash
scRBP rgs \
    --mode       ct \
    --magma      /tools/magma \
    --genes-raw  scz_gwas.genes.raw \
    --sets       regulons_entrez.gmt \
    --id-type    entrez \
    --out        rgs_output/scz_ct_rgs \
    --gene-loc   NCBI38.gene.loc \
    --n-null     1000 \
    --expr-stats ras_ct_output/expr_stats.tsv \
    --threads    16
```

---

## Step 10 — `scRBP trs`

### Purpose

Integrates the **Regulon Activity Score (RAS)** and the **Regulon Gene-Set Score (RGS)** into a final **Trait Relevance Score (TRS)** for each RBP regulon across cells or cell types. The TRS formula balances the two signals while penalizing discordance:

```
TRS = norm(RAS) + norm(RGS) − λ × |norm(RAS) − norm(RGS)|
```

Where:
- `norm()` applies robust min-max normalization using upper quantile `q_hi` to reduce outlier influence.
- `λ` (lambda_penalty) controls how strongly discordance between RAS and RGS is penalized.

RBPs with high TRS are **both** highly active in the cell type **and** genetically linked to the trait.

Empirical p-values and z-scores are computed using matched null distributions from `rgs`. BH-FDR correction is applied in `ct` mode.

**Two modes are supported:**
- **`--mode sc`**: TRS computed per cell.
- **`--mode ct`**: TRS computed per cell type, with FDR correction.

### Usage

```bash
scRBP trs \
    --mode       <sc|ct>             \
    --ras        <aucell_scores.csv>  \
    --rgs-csv    <rgs_results.csv>    \
    --out-prefix <output_prefix>      \
    [options]
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--mode` | str | — | **Yes** | Analysis mode. `sc` = single-cell TRS; `ct` = cell-type-level TRS with FDR correction. |
| `--ras` | str | — | **Yes** | RAS file from `ras`. Formats: `.csv` / `.csv.gz` (cells × regulons) or `.loom`. |
| `--rgs-csv` | str | — | **Yes** | RGS results CSV from `rgs` (`<out>_real.csv`). Must contain `VARIABLE` column and a score column (`mlog10p` or `z`). |
| `--out-prefix` | str | — | **Yes** | Output file prefix. |
| `--rgs-score` | str | `mlog10p` | No | RGS score column to use. `mlog10p` = −log10(MAGMA p-value, recommended); `z` = MAGMA z-score. |
| `--lambda-penalty` | float | 1.0 | No | Penalty coefficient λ for the discordance term. Higher values more strongly penalize RBPs where RAS and RGS disagree. |
| `--q-hi-ras` | float | 0.99 | No | Upper quantile for robust min-max normalization of RAS scores. |
| `--q-hi-rgs` | float | 0.99 | No | Upper quantile for robust min-max normalization of RGS scores. |
| `--celltypes-csv` | str | — | Cond. | CSV with columns `cell_id`, `cell_type`. **Required when `--mode ct`.** |
| `--min_cells_pert_ct` | int | 25 | No | Minimum number of cells a cell type must have to be included in `ct` mode. |
| `--do-fdr` | int | 1 | No | Apply BH-FDR correction in `ct` mode. `1` = apply (default); `0` = skip. |

### Output Files

| File | Description |
|------|-------------|
| `<prefix>.sc.TRS_matrix.csv` | TRS scores matrix (`--mode sc`): rows = cells, columns = regulons |
| `<prefix>.ct.TRS_matrix.csv` | TRS scores matrix (`--mode ct`): rows = cell types, columns = regulons |
| `<prefix>_trs_long.csv` | Long-format TRS table: `regulon`, `cell_type` / `cell`, `TRS`, `p_empirical`, `z_score`, `FDR` |

### Examples

**Single-cell mode (`--mode sc`):**

```bash
scRBP trs \
    --mode         sc \
    --ras          ras_sc_output/aucell_sc.csv \
    --rgs-csv      rgs_output/scz_sc_rgs_real.csv \
    --out-prefix   trs_output/scz_sc_trs \
    --rgs-score    mlog10p \
    --lambda-penalty 1.0
```

**Cell-type mode (`--mode ct`):**

```bash
scRBP trs \
    --mode             ct \
    --ras              ras_ct_output/aucell_ct.csv \
    --rgs-csv          rgs_output/scz_ct_rgs_real.csv \
    --out-prefix       trs_output/scz_ct_trs \
    --rgs-score        mlog10p \
    --lambda-penalty   1.0 \
    --q-hi-ras         0.99 \
    --q-hi-rgs         0.99 \
    --celltypes-csv    cell_to_celltype.csv \
    --min_cells_pert_ct 25 \
    --do-fdr           1
```

---

## Complete Pipeline Example

Below is a full end-to-end example using PBMC data with schizophrenia GWAS:

```bash
# ── Inputs ──────────────────────────────────────────────────
MATRIX=PBMC_healthy.h5ad
RBP_LIST=human_RBP_list.txt
MOTIF_LINKS=motif_to_rbp_hg38_v10.csv
RANKINGS_3UTR=hg38_3UTR_v10_rankings.feather
MAGMA_BIN=/tools/magma
GWAS_RAW=scz_gwas.genes.raw
GENE_LOC=NCBI38.gene.loc

# ── Step 1: Downsample (optional) ───────────────────────────────────────
scRBP getSketch \
    --input $MATRIX --output pbmc_sketch.feather \
    --n_cells 50000 --celltype_col celltype --min_cells_per_type 500

# ── Step 2: Infer GRN (30 seeds) ────────────────────────────
for SEED in $(seq 1 30); do
  scRBP getGRN \
      --matrix pbmc_sketch.feather --rbp_list $RBP_LIST \
      --output grn/grn_seed${SEED} --n_workers 20 --seed $SEED
done

# ── Step 3: Merge GRN runs ───────────────────────────────────
scRBP getMerge_GRN \
    --pattern "grn/grn_seed*_scRBP_gene_GRNs.tsv" --output grn_consensus.tsv \
    --n_present 15 --present_rate 0.5

# ── Step 4: Extract modules ──────────────────────────────────
scRBP getModule \
    --input grn_consensus.tsv --output_merged modules.tsv \
    --top_n_list "5,10,50" --target_top_n "50"

# ── Step 5: Prune with motifs ────────────────────────────────
scRBP getPrune \
    --rbp_targets modules.tsv --motif_rbp_links $MOTIF_LINKS \
    --motif_target_ranks $RANKINGS_3UTR --save_dir pruned/ \
    --rank_threshold 1500 --nes_threshold 3.0 --n_jobs 16

# ── Step 6: Build GMT regulons ───────────────────────────────
scRBP getRegulon \
    --input pruned/ctx_scores.csv \
    --out-symbol regulons_symbol.gmt \
    --out-entrez regulons_entrez.gmt \
    --map-custom $GENE_LOC --min_genes 5

# ── Step 7: Merge regions (if run per-region) ────────────────
scRBP mergeRegulons \
    --base_dir results/ --input regulons_symbol.gmt \
    --output regulons_all.gmt --recursive

# ── Step 8: Compute RAS (cell-type mode) ─────────────────────
scRBP ras \
    --mode ct --matrix pbmc_sketch.feather \
    --regulons regulons_symbol.gmt --out ras_out/ \
    --celltypes-csv pbmc_sketch_cell_to_celltype.csv \
    --emit-expr-stats

# ── Step 9: Run MAGMA RGS (cell-type mode) ───────────────────
scRBP rgs \
    --mode ct --magma $MAGMA_BIN --genes-raw $GWAS_RAW \
    --sets regulons_entrez.gmt --id-type entrez \
    --out rgs_out/scz --gene-loc $GENE_LOC \
    --n-null 1000 --expr-stats ras_out/expr_stats.tsv

# ── Step 10: Compute TRS (cell-type mode) ────────────────────
scRBP trs \
    --mode ct --ras ras_out/aucell_ct.csv \
    --rgs-csv rgs_out/scz_real.csv \
    --out-prefix trs_out/scz_trs \
    --celltypes-csv pbmc_sketch_cell_to_celltype.csv
```

---

## Command Reference

| Step | Command | Key Inputs | Key Output |
|------|---------|-----------|------------|
| 1 | `scRBP getSketch` | `.h5ad` / `.feather` | Downsampled cells |
| 2 | `scRBP getGRN` | Expression matrix, RBP list | `*_scRBP_gene_GRNs.tsv` or `*_scRBP_isoform_GRNs.tsv` |
| 3 | `scRBP getMerge_GRN` | Multiple GRN TSV files (glob) | Consensus GRN TSV |
| 4 | `scRBP getModule` | Consensus GRN TSV | Modules TSV |
| 5 | `scRBP getPrune` | Modules TSV, motif files | Pruned scores (Parquet) |
| 6 | `scRBP getRegulon` | Pruned scores | Regulons GMT (symbol + Entrez) |
| 7 | `scRBP mergeRegulons` | Multiple GMT files | Merged GMT |
| 8 | `scRBP ras` (`--mode sc\|ct`) | Expression matrix, GMT | AUCell scores, RSS matrix |
| 9 | `scRBP rgs` (`--mode sc\|ct`) | MAGMA `.genes.raw`, GMT | RGS scores CSV |
| 10 | `scRBP trs` (`--mode sc\|ct`) | RAS CSV, RGS CSV | TRS scores CSV |

Use `scRBP <command> --help` to see all parameters for any step.

---

## Dependencies

| Category | Packages |
|----------|---------|
| Core numerics | `numpy`, `pandas`, `scipy`, `scikit-learn` |
| Single-cell I/O | `anndata`, `scanpy`, `loompy` |
| Fast I/O | `polars`, `pyarrow` |
| Cell downsampling | `geosketch` |
| GRN inference | `arboreto` (GRNBoost2 / GENIE3) |
| Motif enrichment | `ctxcore`, `pyscenic` |
| Progress display | `tqdm` |
| GWAS enrichment | **MAGMA** binary (external, user-provided) |

---



## Links

- **GitHub**: https://github.com/mayunlong89/scRBP
- **Issues**: https://github.com/mayunlong89/scRBP/issues
- **Website**: [https://mayunlong89.github.io/scRBP.github.io/](https://mayunlong89.github.io/scRBP.github.io/)



## Citation

If you use scRBP in your research, please cite:

> Ma Y. *et al.* *Decoding disease-associated RNA-binding protein-mediated regulatory networks through polygenic enrichment across diverse cellular contexts.* (2026)

---



## License

MIT License. See [LICENSE](LICENSE) for details.

---

