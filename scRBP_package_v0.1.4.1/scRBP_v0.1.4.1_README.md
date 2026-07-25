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




![pypi](https://img.shields.io/badge/pypi-0.1.4.1-green) ![python](https://img.shields.io/badge/python-3.9--3.11-blue) ![license](https://img.shields.io/badge/license-MIT-yellow)

**scRBP (single-cell RNA-binding protein regulon inference)** is a command-line toolkit for reconstructing RNA-binding protein (RBP)-mediated regulatory programs from single-cell transcriptomic data and prioritizing regulons associated with complex traits and disorders. It supports both gene-level networks from short-read data and isoform-level networks from long-read data.

The framework integrates six major analytical components: (i) a curated compendium of RBPs and clustered RBP-binding motifs assembled from public resources; (ii) motif-guided transcriptome-wide rankings of candidate RBP targets at gene and isoform resolution; (iii) inference of RBP–gene or RBP–isoform association networks from single-cell expression data; (iv) motif-based refinement of candidate edges to define high-confidence regulons and quantify regulon activity scores (RASs); (v) parallel common- and rare-variant models that map genetic association signals to regulons and derive regulon-level genetic association scores (RGSs); and (vi) integration of RAS and RGS into a trait-relevance score (TRS) for each regulon within each cellular context, with significance evaluated against matched null regulons using Monte Carlo sampling.

---

## What scRBP Does

RBPs regulate multiple layers of post-transcriptional gene control, including RNA splicing, localization, stability, and translation. scRBP provides an end-to-end workflow to:

- **Infer** candidate RBP–gene or RBP–isoform association networks from single-cell transcriptomes.
- **Refine** candidate RBP–target edges using sequence-motif evidence to define high-confidence regulons.
- **Quantify** regulon activity at single-cell or cell type resolution using the Regulon Activity Score (RAS).
- **Evaluate** common- or rare-variant enrichment for each regulon using the Regulon-level Genetic Association Score (RGS).
- **Prioritize** trait-relevant RBP regulons in specific cellular contexts by integrating RAS and RGS into the Trait-Relevance Score (TRS).

---

## New in v0.1.4.1

Version `0.1.4.1` expands scRBP with a denser GRN inference workflow and parallel common- and rare-variant genetic models:

- **`scRBP getMetacell`** constructs mini-metacells from transcriptionally similar cells to reduce sparsity and improve RBP–target network inference.
- **`scRBP ras`** supports single-cell and cell type-level activity analysis and can emit per-gene expression statistics for matched-null construction.
- **`scRBP rgs`** implements the common-variant model using MAGMA competitive gene-set analysis.
- **`scRBP rgs_rare`** implements a parallel rare-variant model using gene-level TADA, logBF, burden, SAIGE-GENE+, REGENIE, STAAR-O, or compatible association statistics.
- **`scRBP trs`** accepts either common-variant or rare-variant RGS results, enabling calculation of corresponding common-variant and rare-variant Trait-Relevance Scores.

---

## Pipeline at a Glance

```
Raw single-cell data (.h5ad / .feather)
          │
          ▼
[Step 1]  scRBP getSketch        ── Stratified GeoSketch downsampling (optional)
          │
          ▼
[Step 2]  scRBP getMetacell      ── Aggregate similar cells into mini-metacells for GRN inference
          │                          (optional; recommended for the GRN branch when
          │                           dropout is heavy. Complementary to getSketch)
          ▼
[Step 3]  scRBP getGRN           ── GRNBoost2/GENIE3 RBP→gene or RBP→isoform inference
          │                          (run multiple random seeds; 30 runs recommended)
          ▼
[Step 4]  scRBP getMerge_GRN     ── Merge multi-seed GRNs into a consensus network
          │
          ▼
[Step 5]  scRBP getModule        ── Extract regulon candidates (Top-N / percentile)
          │
          ▼
[Step 6]  scRBP getPrune         ── Motif-enrichment pruning
          │
          ▼
[Step 7]  scRBP getRegulon       ── Export pruned regulons to GMT format
          │
          ▼
[Step 8]  scRBP mergeRegulons    ── Merge region-specific GMT files
          │                          (3'UTR / 5'UTR / CDS / Introns)
          ▼
[Step 9]  scRBP ras              ── Regulon Activity Score (RAS) per cell / cell type
          │                          (--mode sc | --mode ct)
          ▼
[Step 10] scRBP rgs              ── Common-variant Regulon-level Genetic Association Score (RGS)
          │                          via MAGMA (--mode sc | --mode ct)
          ▼
[Step 11] scRBP rgs_rare         ── Rare-variant RGS via competitive gene-set regression
          │                          (TADA / logBF / burden; --mode sc | --mode ct)
          ▼
[Step 12] scRBP trs              ── Trait-Relevance Score (RAS–RGS integration)
                                     (--mode sc | --mode ct)
```

> Steps 2 (**getMetacell**) and 11 (**rgs_rare**) are major additions in v0.1.4.1. `getMetacell`
> is optional and is used only for the GRN inference branch when the input is
> particularly sparse; `getSketch` remains the recommended thinning strategy for
> the activity (RAS) branch. `rgs_rare` runs in parallel to `rgs` and can be fed
> back into `scRBP trs` in exactly the same way to obtain a **rare-variant** TRS.

---

## Installation

### Requirements

- Python **3.9, 3.10, or 3.11** (Python 3.12+ not yet supported by `arboreto`)
- MAGMA binary (external; required only for Step 10, `scRBP rgs`)

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

### Install MAGMA (for Step 10 only)

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

## Step 1 — `scRBP getSketch` (optional)

### Purpose

Performs diversity-preserving downsampling of large single-cell expression datasets using **GeoSketch**. For `.h5ad` input, sketching is stratified by the cell-type labels stored in `adata.obs`, allowing each annotated population to contribute at least a user-defined minimum number of cells while keeping the total close to the requested target. For `.feather` input, GeoSketch is applied globally because cell-level annotations are not embedded in the matrix.

GeoSketch operates on a reduced-dimensional representation generated with PCA or TruncatedSVD and preferentially retains cells that span the transcriptional state space. This reduces the computational cost of downstream analyses while preserving heterogeneous and relatively rare cellular states more effectively than uniform random sampling.

> For datasets with several hundred thousand cells or more, sketching to approximately 50,000 cells is a practical starting point. The appropriate target should be adjusted according to dataset complexity, the abundance of rare populations, and available computing resources.

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
| `--output` | str | — | **Yes** | Output file path. The output format is inferred from extension: `.h5ad`, `.csv`, `.feather`, `.npz`. |
| `--n_cells` | int | 50000 | No | Target total number of cells to retain after sketching. |
| `--n_pca` | int | 100 | No | Number of PCA (TruncatedSVD) components before running GeoSketch. |
| `--celltype_col` | str | `"celltype"` | No | Column name in `adata.obs` for cell-type labels. Used for stratified sampling (`.h5ad` only). |
| `--min_cells_per_type` | int | 50 | No | Minimum number of cells to retain per cell type. If a type has fewer cells, all are kept. |
| `--seed` | int | 42 | No | Random seed for reproducibility. |

### Outputs

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

## Step 2 — `scRBP getMetacell` (optional)

### Purpose

Aggregates transcriptionally similar single cells into **mini-metacells** to reduce expression sparsity and strengthen the co-expression signal used for network inference. By default, approximately 10 neighboring cells are combined into each metacell. This is particularly useful for RBP network reconstruction, where regulatory effects may have a smaller expression dynamic range than conventional transcription-factor programs.

Cell similarity is computed from a library-size-normalized, log1p-transformed low-dimensional embedding, whereas the metacell expression profile is calculated from the original linear-scale counts using either summation or averaging. When cell-type annotations are available, metacells are constructed within each cell type by default to prevent pooling across biologically distinct populations.

`getSketch` and `getMetacell` serve complementary purposes:

- **`getSketch`** selects a diversity-preserving subset of real cells and is recommended for the RAS branch, where single-cell resolution should be retained.
- **`getMetacell`** aggregates related cells into denser pseudo-profiles and is intended for the GRN branch (`getGRN → getModule → getPrune`).

For very large annotated datasets, the two operations can be combined by first creating an `.h5ad` sketch that preserves `adata.obs`, followed by within-cell-type metacell construction.

### Usage

```bash
scRBP getMetacell \
    --input  <input_file>  \
    --output <output_file> \
    [options]
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--input` | str | — | **Yes** | Input expression file. `.h5ad` (with a cell-type column for stratified pooling) or `.feather` (gene × cell). |
| `--output` | str | — | **Yes** | Output file (gene × metacell). The output format is inferred: `.h5ad`, `.csv`, `.feather`, `.npz`. Directly consumable by `getGRN`. |
| `--metacell_size` | int | 10 | No | Target number of cells pooled per metacell (atlas typical 10–15). |
| `--method` | str | `knn` | No | Cell-grouping strategy: `knn` (greedy nearest-neighbour, most faithful to the metacell concept), `kmeans` (MiniBatchKMeans, scales well), or `random` (fastest baseline; ignores similarity). |
| `--agg` | str | `sum` | No | Aggregate pooled cells by `sum` (default) or `mean` of the original counts. |
| `--within_celltype` | flag | on | No | Pool only within each cell type. On by default; disable via `--global_pooling`. |
| `--global_pooling` | flag | off | No | Disable within-cell-type pooling (pool across all cells). Not recommended when a cell-type annotation is available. |
| `--celltype_col` | str | `"celltype"` | No | Column in `adata.obs` holding cell-type labels (`.h5ad` input). |
| `--min_metacell_size` | int | 1 | No | Merge metacells smaller than this into the nearest one. Recommended `~metacell_size // 2` for cleaner pooling. |
| `--n_pca` | int | 50 | No | PCA components for the similarity embedding (`knn`/`kmeans`). |
| `--seed` | int | 42 | No | Random seed for reproducibility. |
| `--save_members` | flag | off | No | Also write `<out>_metacell_members.csv` mapping each metacell to its source cells. |

### Outputs

| File | Description |
|------|-------------|
| `<output>` | Metacell expression matrix (gene × metacell). Directly consumable by `getGRN`. |
| `<output_prefix>_metacell_to_celltype.csv` | Mapping table (`cell`, `cell_type`, `n_cells`). Column names match `getSketch`, so it is drop-in for `ras --celltypes-csv`. |
| `<output_prefix>_metacell_summary.csv` | Per-cell type size summary (`n_metacells`, `median_cells`, `min_cells`, `max_cells`, `total_cells`). |
| `<output_prefix>_metacell_members.csv` | *(optional, with `--save_members`)* Full source-cell membership table. |

### Examples

**A. Within-cell type metacells (recommended) → feed `getGRN`:**

```bash
scRBP getMetacell \
    --input             lung_lineage.h5ad \
    --output            lung_metacell.feather \
    --metacell_size     10 \
    --method            knn \
    --within_celltype \
    --celltype_col      cell_type \
    --min_metacell_size 5
```

**B. Very large data — combine with `getSketch` in the GRN branch:**

```bash
# First, geometric-sketch pre-thin
scRBP getSketch  --input huge.h5ad --output sketch.h5ad --n_cells 300000 --celltype_col cell_type
# Then, densify with within-cell-type metacells
scRBP getMetacell --input sketch.h5ad --output mc.feather --metacell_size 12 --celltype_col cell_type
```

---

## Step 3 — `scRBP getGRN`

### Purpose

Infers an RBP–target association network from gene- or isoform-level expression data using tree-based regression. The command supports **GRNBoost2**, a gradient-boosting implementation optimized for speed, and **GENIE3**, a random-forest-based alternative. For each RBP–target pair, the resulting importance score quantifies the contribution of the RBP expression profile to prediction of the target expression profile.

When correlation analysis is enabled, scRBP also calculates the Spearman correlation for each edge and assigns a correlation-based mode label (`activating`, `repressing`, or `neutral`). These labels describe the direction of expression covariation and should not by themselves be interpreted as proof of direct biochemical activation or repression.

Two inference modes are available:

- **`gene` mode** uses RBP gene expression as the regulator signal and genes as candidate targets.
- **`isoform` mode** aggregates all detected isoforms of each RBP into a gene-level regulator signal while retaining individual transcript isoforms as candidate targets. This mode requires `--isoform_annotation`.

> Because tree-based network inference is stochastic, running 30 independent seeds and combining them with `getMerge_GRN` is recommended for a stable consensus network.

### Usage

```bash
scRBP getGRN \
    --matrix   <expression_matrix> \
    --rbp_list <rbp_list_file>     \
    --output   <output_prefix>     \
    --mode     gene|isoform        \
    [options]
```

### Automatic output filenames

| Mode | Output filename |
|------|----------------|
| `gene` | `<prefix>_scRBP_gene_GRNs.tsv` |
| `isoform` | `<prefix>_scRBP_isoform_GRNs.tsv` |

### Core parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|:--------|----------|-------------|
| `--matrix` | str | — | **Yes** | Expression matrix in `.csv`, `.csv.gz`, `.feather`, or `.loom` format. Rows represent genes or isoforms and columns represent cells or metacells. |
| `--rbp_list` | str | — | **Yes** | Plain-text file of RBP gene symbols (one per line). |
| `--output` | str | — | **Yes** | Output prefix. Mode suffix appended automatically. |
| `--mode` | str | `gene` | No | `gene` — RBP→Gene network; `isoform` — RBP→Isoform network. |
| `--method` | str | `grnboost2` | No | `grnboost2` (recommended) or `genie3`. |
| `--n_workers` | int | all CPUs | No | Number of parallel worker processes. |
| `--correlation` | bool | `True` | No | Calculate Spearman correlation for each edge and assign a correlation-based mode label. |
| `--threshold` | float | 0.03 | No | Minimum absolute Spearman correlation. Edges with `|r| ≤ threshold` are removed when correlation filtering is enabled. |
| `--seed` | int | 1234 | No | Random seed for the tree ensemble. Change across runs for consensus merging. |
| `--batch_size` | int | 10 | No | Number of processing batches. Lower values can reduce peak resource usage but may increase overhead; `5` is a practical setting for memory-constrained runs. |



### Isoform-mode parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--isoform_annotation` | str | `None` | **Yes** (isoform) | TSV/CSV mapping isoform/transcript IDs → gene symbols. |
| `--rbp_agg_method` | str | `sum` | No | How to aggregate multiple isoforms of the same RBP: `sum`, `mean`, `max`. |
| `--remove_self_targets` | bool | `True` | No | Remove edges where the target isoform belongs to the same gene as the RBP. |
| `--min_target_cells_expressed` | int | 10 | No | Minimum number of cells with expression > 0 to keep a target isoform. |
| `--min_target_mean_expr` | float | 0.01 | No | Minimum mean expression across all cells to keep a target isoform. |



### Outputs

| File | Description |
|------|-------------|
| `<prefix>_scRBP_gene_GRNs.tsv` | Gene-level network with columns `RBP`, `Gene`, `Importance`, `Correlation`, and `Mode`. |
| `<prefix>_scRBP_isoform_GRNs.tsv` | Isoform-level network with columns `RBP`, `Isoform`, `Importance`, `Correlation`, and `Mode`. |
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
      --n_workers 10 \
      --correlation True \
      --batch_size 5 \
      --threshold 0.03 \
      --seed      ${SEED}
done
# Output: grn_seed1_scRBP_gene_GRNs.tsv, grn_seed2_scRBP_gene_GRNs.tsv, ...
```

**Isoform mode (30 seeds):**

```bash
for SEED in $(seq 1 30); do
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
      --batch_size 5 \
      --threshold 0.03 \
      --n_workers                  10 \
      --seed                       ${SEED}
done
```

---

## Step 4 — `scRBP getMerge_GRN`

### Purpose

Combines GRN files generated from multiple random seeds into a single **consensus network**. For each RBP–gene or RBP–isoform edge, the command averages the available importance and correlation values across runs, records the number of runs in which the edge was detected as `n_present`, and calculates its stability as:

```text
presence_rate = n_present / total_runs
```

Optional thresholds can then be applied to retain edges with sufficient correlation strength, recurrence across runs, or presence rate. When more than one threshold is specified, an edge must satisfy all requested criteria.

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
| `--pattern` | str | — | **Yes** | Glob pattern used to match all GRN files generated from different random seeds, e.g. `"grn_seed*.tsv"`. The pattern should be enclosed in quotation marks to prevent shell expansion. |
| `--output` | str | — | **Yes** | Output path for the merged consensus GRN file in TSV format. |
| `--corr-threshold` | float | 0.0 | No | Minimum absolute mean correlation required to retain an edge. Edges with averaged `|Correlation|` less than or equal to this threshold are removed. |
| `--n_present` | int | 0 | No | Minimum number of seed runs in which an edge must be detected. |
| `--present_rate` | float | 0.0 | No | Minimum proportion of seed runs in which an edge must be detected, calculated as `n_present / total_runs`. A value of `0.3` retains edges observed in at least 30% of runs. |

### Example

```bash
scRBP getMerge_GRN \
    --pattern "grn_seed*.tsv" \
    --output  grn_consensus.tsv \
    --n_present 10 \
    --present_rate 0.3
```

> This command retains edges detected in at least 10 runs and in at least 30% of all runs. Because both filters are active, an edge must satisfy both criteria.

---

## Step 5 — `scRBP getModule`

### Purpose

Extracts candidate RBP regulon modules from the consensus GRN using complementary edge-selection strategies. These strategies capture different aspects of network structure and reduce dependence on a single arbitrary cutoff:

- **Top-N RBPs per target** (`topN_per_gene_N`): for each target gene or isoform, retain the N RBPs with the highest importance scores.
- **Top-N targets per RBP** (`top_target_N`): for each RBP, retain its N highest-importance targets.
- **RBP-specific percentile** (`pctX_per_rbp`): for each RBP, retain targets above a specified within-RBP importance percentile.

A global importance threshold is applied before strategy-specific selection. The merged output can contain modules derived from multiple strategies, which are subsequently evaluated independently during motif-based pruning.

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
| `--input` | str | — | **Yes** | Consensus GRN in TSV format. Required columns are `RBP`, the target column (`Gene` or `Isoform`), and `Importance`. |
| `--output_merged` | str | — | **Yes** | Output TSV containing candidate modules from all requested selection strategies. |
| `--importance_threshold` | float | 0.005 | No | Global minimum importance cutoff applied before any strategy. Strategy-specific thresholds can be adjusted according to the desired network stringency. |
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

## Step 6 — `scRBP getPrune`

### Purpose

Refines candidate RBP modules by testing whether their predicted targets are enriched for sequence motifs associated with the corresponding RBP. The command uses the `ctxcore` rank-based enrichment framework to integrate the expression-derived GRN with independent motif-binding evidence.

For each RBP–strategy combination, `getPrune`:

1. maps the RBP to one or more associated motifs using `--motif_rbp_links`;
2. evaluates the candidate target set against the motif-specific gene or isoform ranking matrix;
3. calculates enrichment statistics, including AUC, normalized enrichment score (NES), and leading-edge targets; and
4. retains regulons that satisfy the specified rank, AUC, NES, and minimum-size thresholds.

This step removes candidate edges that lack motif support and yields higher-confidence regulons supported by both co-expression and sequence evidence.

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
| `--motif_rbp_links` | str | — | **Yes** | CSV/TSV linking motif IDs to RBP names (from [scRBP resource database](https://resources.gandallab.org/scRBP/)). |
| `--motif_target_ranks` | str | — | **Yes** | Feather file containing the motif × gene or motif × isoform ranking matrix ([scRBP resource database](https://resources.gandallab.org/scRBP/)). |
| `--save_dir` | str | — | **Yes** | Output directory. Created if it does not exist. |
| `--rank_threshold` | int | 1500 | No | Maximum target rank included in the recovery curve used for enrichment scoring. |
| `--auc_threshold` | float | 0.05 | No | Minimum enrichment AUC required to retain an RBP–strategy module. |
| `--nes_threshold` | float | 3.0 | No | Minimum normalized enrichment score required to retain a regulon. |
| `--min_genes` | int | 20 | No | Minimum number of target genes required for an RBP regulon to be scored. |
| `--n_jobs` | int | all CPUs | No | Number of parallel worker processes. |
| `--weight_mode` | str | `importance` | No | Gene weighting mode. `importance` uses GRN importance scores; `equal` weights all genes equally. |

### Example

```bash
scRBP getPrune \
    --rbp_targets        modules.tsv \
    --motif_rbp_links    rbp_motif_annotation.csv \
    --motif_target_ranks hg38_3UTR_rankings.feather \
    --save_dir           pruned_results/ \
    --rank_threshold     1500 \
    --auc_threshold      0.05 \
    --nes_threshold      3.0 \
    --min_genes          20 \
    --n_jobs             16
```

---

## Step 7 — `scRBP getRegulon`

### Purpose

Converts the pruned enrichment output from `getPrune` into standard **GMT gene-set files** for downstream activity and genetic analyses. Each GMT line represents one RBP regulon and contains the regulon identifier followed by its retained target genes.

Two identifier formats are generated:

- **Gene-symbol GMT**, typically used by `scRBP ras` because expression matrices are usually indexed by gene symbols.
- **Entrez-ID GMT**, typically used by `scRBP rgs` for MAGMA gene-set analysis.

The command supports custom, HGNC, or NCBI mapping resources and can optionally remove unmapped genes or regulons that become empty after identifier conversion.

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
| `--input` | str | — | **Yes** | Input CSV containing pruned enrichment results from `getPrune`. |
| `--out-symbol` | str | — | **Yes** | Output path for the gene-symbol GMT file. |
| `--out-entrez` | str | — | **Yes** | Output path for the Entrez-ID GMT file. |
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

## Step 8 — `scRBP mergeRegulons`

### Purpose

Combines transcript-region-specific GMT files—such as 3′UTR, 5′UTR, CDS, and intronic regulons—into a single GMT file. The command searches the specified directory structure for matching GMT files, orders them according to a user-defined region priority, and concatenates the regulons into one downstream-ready resource.

When the same RBP is represented in multiple transcript regions, `--append_region_to_setname` is recommended so that each region-specific regulon retains a unique identifier (for example, `RBFOX1_3UTR` and `RBFOX1_Introns`). Recursive mode can process multiple tissues or datasets in a single command, and an optional summary table reports the number of regulons contributed by each region.

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
| `--output` | str | — | **Yes** | Output path for the merged GMT file. |
| `--recursive` | flag | off | No | Process multiple tissue/dataset directories under `--base_dir`. |
| `--tissue_glob` | str | `"z_GRNBoost2_*_30times"` | No | Glob pattern for tissue/dataset directories when `--recursive` is used. |
| `--region_glob` | str | `"Results_final_*_RBP_top1500_*"` | No | Glob pattern for region subdirectories. |
| `--append_region_to_setname` | flag | off | No | Append the transcript-region label (for example, `_3UTR`) to each regulon name to preserve region identity. |
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

## Step 9 — `scRBP ras`

### Purpose

Computes **Regulon Activity Scores (RASs)** from a gene-expression matrix and a collection of RBP regulons. At single-cell resolution, scRBP uses **AUCell** to quantify whether the targets of each regulon are preferentially represented near the top of each cell’s ranked expression profile. Because AUCell is rank based, the score is relatively robust to differences in library size and expression scale.

Two analysis modes are supported:

- **`--mode sc`** returns the cell-by-regulon AUCell activity matrix.
- **`--mode ct`** summarizes single-cell activity across annotated cell types and calculates regulon specificity using an entropy-based Jensen–Shannon divergence framework. This mode requires a cell-to-cell-type mapping file.

The command can also emit per-gene expression statistics (`mean_expr` and `pct_detected`), which are used by the common- and rare-variant RGS modules to construct expression-matched null regulons.

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
| `--mode` | str | `ct` | No | Analysis mode: `sc` returns single-cell AUCell scores; `ct` additionally summarizes activity and specificity across annotated cell types. |
| `--matrix` | str | — | Cond. | Expression matrix file. Formats: `.csv`, `.csv.gz`, `.feather`, `.loom`. Rows = genes, columns = cells. Required unless `--aucell-in` is provided. |
| `--regulons` | str | — | Cond. | Regulon file. Formats: `.gmt` / `.gmt.gz` (Symbol) or `.pkl`. Required unless `--aucell-in` is provided. |
| `--aucell-in` | str | — | No | Precomputed cell-by-regulon AUCell matrix. When provided, the AUCell calculation step is skipped. |
| `--out` | str | — | **Yes** | Output path prefix or directory. |
| `--out_format` | str | `loom` | No | Output format for RAS scores. Choices: `csv`, `loom`, `both`. |
| `--no-csv` | flag | off | No | Disable CSV output even if `--out_format` includes csv. Useful when the LOOM output is enough for downstream `trs`. |
| `--no-loom` | flag | off | No | Disable LOOM output even if `--out_format` includes loom. |
| `--n_workers` | int | 4 | No | Number of parallel workers for RAS computation. |
| `--celltypes-csv` | str | — | Cond. | CSV file with cell barcode and cell-type columns. **Required when `--mode ct`.** |
| `--cell-col` | str | auto | No | Column name for cell barcodes in `--celltypes-csv`. Auto-detected if not specified. |
| `--ctype-col` | str | auto | No | Column name for cell-type labels in `--celltypes-csv`. Auto-detected if not specified. |
| `--min_genes` | int | 1 | No | Minimum number of regulon targets that must be present in the expression matrix. |
| `--emit-expr-stats` | flag | on | No | Compute and save per-gene expression statistics (`mean_expr`, `pct_detected`). Enabled by default. Used as input to `rgs`. |
| `--no-expr-stats` | flag | — | No | Disable expr-stats output. Mutually exclusive with `--emit-expr-stats`. |
| `--expr-stats-out` | str | auto | No | Output path for expression statistics TSV. Auto-generated from `--out` if not specified. |

### Outputs

| File | Description |
|------|-------------|
| `aucell_sc.csv` | Single-cell AUCell matrix: rows are cells and columns are regulons (`--mode sc`). |
| `aucell_ct.csv` | Cell type–level activity matrix: rows are cell types and columns are regulons (`--mode ct`). |
| `rss.tsv` | Regulon specificity scores for each cell type (`--mode ct` only). |
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

## Step 10 — `scRBP rgs`

### Purpose

Computes a **common-variant Regulon-level Genetic Association Score (RGS)** using MAGMA competitive gene-set analysis. For each RBP regulon, MAGMA tests whether its target genes show stronger GWAS association than genes outside the regulon while accounting for the gene-level covariates included in the MAGMA model.

To obtain an empirical reference distribution, scRBP can generate matched null regulons for each observed regulon. Null genes are matched on MAGMA-derived genomic features—including the number of SNPs and model parameters—and, when expression statistics are supplied, on mean expression and detection rate. This design reduces confounding by regulon size, genomic annotation density, and expression-related gene-selection biases.

The RGS is defined at the regulon level; cellular specificity is introduced by pairing the RGS with cell- or cell type-resolved RAS values in the downstream `trs` step. In `ct` mode, expression statistics from `ras --emit-expr-stats` are used for expression-aware matched-null construction.

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
| `--mode` | str | — | **Yes** | Analysis mode. `sc` performs the observed-regulon analysis; `ct` enables the cell type–aware matched-null workflow. |
| `--magma` | str | — | **Yes** | Full path to the MAGMA binary executable. |
| `--genes-raw` | str | — | **Yes** | MAGMA `.genes.raw` file from MAGMA gene analysis step. |
| `--sets` | str | — | **Yes** | GMT regulon file. Use gene-symbol GMT (`--id-type symbol`) or Entrez-ID GMT (`--id-type entrez`, default). |
| `--id-type` | str | `entrez` | No | Gene ID type in GMT file: `entrez` or `symbol`. |
| `--out` | str | — | **Yes** | Output file prefix. |
| `--gene-loc` | str | — | No | MAGMA NCBI gene location file. Required when `--id-type symbol`. |
| `--n-null` | int | 1000 | No | Number of matched null regulons per real regulon. Larger values provide a more stable empirical null distribution at increased computational cost. |
| `--seed` | int | 2025 | No | Random seed for null regulon sampling. |
| `--q-bins` | int | 10 | No | Number of quantile bins for null matching dimensions. |
| `--min_genes` | int | 0 | No | Minimum number of genes a regulon must have (after overlap with MAGMA universe) to be tested. |
| `--threads` | int | `None` | No | Number of threads for MAGMA. |
| `--expr-stats` | str | — | No | Expression statistics generated by `ras --emit-expr-stats`; used with genomic covariates for four-dimensional matched-null sampling. |

### Outputs

| File | Description |
|------|-------------|
| `<out>_real.csv` | Parsed MAGMA results for real regulons: `VARIABLE`, `NGENES`, `BETA`, `BETA_STD`, `SE`, `P`, `z`, `mlog10p` |
| `<out>_null_summary.tsv` | Empirical null distribution summary per regulon |
| `<out>.gsa.out` | Raw MAGMA gene-set analysis output |

### Examples

**Observed-regulon mode (`--mode sc`):**

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

**Matched-null cell type mode (`--mode ct`):**

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

## Step 11 — `scRBP rgs_rare`

### Purpose

Computes a **rare-variant Regulon-level Genetic Association Score (RGS)** from externally generated gene-level rare-variant association statistics. Supported inputs include Bayesian scores from TADA or extTADA and frequentist gene-level results from burden, SCHEMA-style, SAIGE-GENE+, REGENIE, STAAR-O, or related frameworks.

The command first converts each gene-level statistic into a common rare-variant score: `−log10(P)` for P-value inputs, `logBF` for Bayesian inputs, or a user-supplied direct score. The upper tail can be winsorized before z-standardization to reduce the influence of extreme values. For each regulon, scRBP then fits a competitive gene-set regression:

```text
rare_score_g = beta0 + beta_s × I(g in regulon_s) + gamma × C_g + error_g
```

Here, `C_g` represents coding opportunity, implemented as log-transformed and standardized union CDS length. Expression covariates are deliberately excluded from the primary regression and are used only for constructing matched null regulons in `ct` mode. The standardized regulon coefficient, `RGS_z = beta_s / SE(beta_s)`, provides the principal rare-variant RGS.

- **`--mode sc`** evaluates the observed regulons and writes a regulon-level RGS table.
- **`--mode ct`** additionally generates null regulons matched on CDS length, mean expression, and detection rate, producing REAL and NULL rows compatible with downstream empirical TRS analysis.

> `rgs_rare` does not perform primary rare-variant association testing. It requires gene-level summary statistics generated by an upstream rare-variant method.

### Usage

```bash
scRBP rgs_rare \
    --mode          <sc|ct>              \
    --rare-summary  <rare_gene_scores>   \
    --rare-gene-col <gene_column>        \
    --score-mode    <pvalue|logbf|direct>\
    --sets          <regulons.gmt>       \
    --out           <output_prefix>      \
    [options]
```

### Core parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--mode` | str | — | **Yes** | `sc` = real regulons only; `ct` = REAL + matched NULLs (feeds `trs`). |
| `--rare-summary` | str | — | **Yes** | Gene-level rare-variant summary CSV/TSV (TADA, burden, SAIGE-GENE+, STAAR-O, …). |
| `--rare-gene-col` | str | — | **Yes** | Gene identifier column in `--rare-summary`. |
| `--rare-id-type` | str | `symbol` | No | Gene ID type in the rare summary: `symbol` or `entrez`. |
| `--score-mode` | str | — | **Yes** | How to build the gene-level rare score: `pvalue` (uses `--p-col` → `-log10(P)`), `logbf` (uses `--logbf-col`), or `direct` (uses `--score-col` verbatim). |
| `--p-col` / `--logbf-col` / `--score-col` | str | — | Cond. | Column name in `--rare-summary` matching the chosen `--score-mode`. |
| `--top-winsor` | float | 0.01 | No | Upper-tail winsorization fraction for gene-level rare scores. |
| `--sets` | str | — | **Yes** | Regulon GMT (Symbol or Entrez). |
| `--id-type` | str | `symbol` | No | Gene ID type in `--sets`. |
| `--out` | str | — | **Yes** | Output file prefix. |
| `--gene-loc` | str | — | Cond. | MAGMA `NCBI*.gene.loc` for Symbol ↔ Entrez mapping. Required when input ID types differ. |
| `--min_genes` | int | 0 | No | Minimum overlap size for a regulon to be tested. |
| `--max-regulon-frac` | float | 0.5 | No | Skip regulons that cover more than this fraction of the gene universe. |

### CDS-covariate parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--cds-length` | str | — | Optional CDS length table (CSV/TSV). If not provided, the CDS column is looked up in `--rare-summary`. |
| `--cds-gene-col` | str | `symbol` | Gene column in `--cds-length`. |
| `--cds-id-type` | str | `symbol` | Gene ID type in `--cds-length`. |
| `--cds-col` | str | `union_cds_length` | CDS covariate column name. |
| `--cds-scale` | str | `raw` | `raw` — raw union CDS length, internally `log1p` + z-scored within the analyzed universe. `z_log` — already a z-standardized log-CDS covariate, used as-is (no double-standardisation). |

### Expression-statistics parameters (`--mode ct`)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--expr-stats` | str | — | Precomputed expression statistics TSV with columns `symbol`, `mean_expr`, `pct_detected`. Reuses the file emitted by `scRBP ras --emit-expr-stats`. |
| `--emit-expr-stats` | bool | `False` | If `True` and `--expr-stats` is missing, compute from `--matrix-stats` and save. |
| `--matrix-stats` | str | — | Expression matrix used only when computing expr-stats on the fly. |
| `--expr-stats-out` | str | `<out>_expr_stats.tsv` | Path for the computed expr-stats file. |
| `--dtype` | str | `float32` | Numeric data type used when loading `--matrix-stats`. |
| `--allow-missing-expr-stats` | flag | off | Fill missing `mean_expr`/`pct_detected` with 0 (not recommended). |

### Matched-null parameters (`--mode ct`)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--n-null` / `--null` | int | 1000 | Number of matched null regulons per real regulon. |
| `--seed` | int | 2025 | Seed for null sampling. |
| `--q-bins` | int | 5 | Quantile bins for matched-null construction. Use `10` for stricter matching on large universes. |
| `--exclude-self` | bool | `True` | Exclude real-regulon genes when sampling nulls. |
| `--min-bucket-size` | int | 5 | Minimum matched-bucket size before falling back to a coarser stratum. |
| `--fdr-method` | str | `BH` | FDR method for the empirical audit table (`BH` or `BY`). |
| `--save-null-gmt` | str | — | Path to save REAL+NULL GMT in the working ID type. |
| `--save-null-gmt-symbol` | str | — | Also save REAL+NULL gene-symbol GMT. |
| `--save-null-gmt-entrez` | str | — | Also save REAL+NULL Entrez-ID GMT. |

### Outputs

| File | Description |
|------|-------------|
| `<out>.gsa_RGS.csv` | Regulon-level RGS table. `sc`: `Regulon, GeneSet, NGENES, BETA, BETA_STD, SE, P, RGS_z, RGS_mlog10P`. `ct`: additionally `RBP, SET_KIND, NULL_ID`. |
| `<out>_gene_scores.tsv` | Per-gene rare score with winsorized / z-scored columns and matched CDS covariate. |
| `<out>_REAL_PLUS_NULLS.<idtype>.gmt` | *(ct only)* REAL + NULL regulons, compatible with `scRBP ras` if you want to recompute RAS on the null sets. |
| `<out>.null_index.tsv` | *(ct only)* Index of every REAL/NULL set with sizes. |
| `<out>_empirical.csv` | *(ct only)* Audit table summarizing REAL vs NULL RGS distributions (`P_empirical`, `z_empirical`, `FDR_empirical`, `P_param`, `FDR_param`). |

### Examples

**Observed-regulon mode with TADA-like logBF input:**

```bash
scRBP rgs_rare \
    --mode         sc \
    --rare-summary tada_asd_gene_scores.tsv \
    --rare-gene-col gene_symbol \
    --rare-id-type symbol \
    --score-mode   logbf \
    --logbf-col    logBF \
    --sets         regulons_symbol.gmt \
    --id-type      symbol \
    --out          rgs_rare_out/asd_sc_rare
```

**Matched-null cell type mode with a P-value input (SAIGE-GENE+ / burden), reusing
`ras --emit-expr-stats`:**

```bash
scRBP rgs_rare \
    --mode          ct \
    --rare-summary  saige_gene_burden.tsv \
    --rare-gene-col gene \
    --rare-id-type  symbol \
    --score-mode    pvalue \
    --p-col         P \
    --cds-length    gene_union_cds.tsv \
    --cds-col       union_cds_length \
    --cds-scale     raw \
    --sets          regulons_symbol.gmt \
    --id-type       symbol \
    --gene-loc      NCBI38.gene.loc \
    --expr-stats    ras_ct_output/scz_ct_expr_stats.tsv \
    --n-null        1000 \
    --q-bins        10 \
    --out           rgs_rare_out/scz_ct_rare
```

The resulting `<out>.gsa_RGS.csv` can be handed directly to `scRBP trs`
via `--rgs-csv` to obtain a **rare-variant TRS**.

---

## Step 12 — `scRBP trs`

### Purpose

Integrates the **Regulon Activity Score (RAS)** with the corresponding common- or rare-variant **Regulon-level Genetic Association Score (RGS)** to calculate a **Trait-Relevance Score (TRS)** for each regulon in each cell or cell type.

RAS and RGS are first robustly normalized, after which scRBP rewards concordant support from both modalities and penalizes imbalance between them:

```text
TRS = norm(RAS) + norm(RGS) − λ × |norm(RAS) − norm(RGS)|
```

The penalty coefficient `λ` controls the strength of the discordance penalty. A high TRS therefore identifies regulons that are both transcriptionally active in a given cellular context and genetically associated with the trait, rather than regulons driven by only one component.

Matched null regulons generated by `rgs` or `rgs_rare` are used to estimate empirical P values and z scores. In cell type mode, scRBP can additionally apply Benjamini–Hochberg false-discovery-rate correction across tested regulon–cell type combinations.

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
| `--mode` | str | — | **Yes** | Analysis mode. `sc` = single-cell TRS; `ct` = cell type–level TRS with FDR correction. |
| `--ras` | str | — | **Yes** | RAS file from `ras`. Formats: `.csv` / `.csv.gz` (cells × regulons) or `.loom`. |
| `--rgs-csv` | str | — | **Yes** | Regulon-level RGS table generated by `rgs` or `rgs_rare`. The file must contain a regulon identifier and the score column selected with `--rgs-score`. |
| `--out-prefix` | str | — | **Yes** | Output file prefix. |
| `--rgs-score` | str | `mlog10p` | No | RGS column used for integration. Common-variant outputs typically use `mlog10p` or `z`; rare-variant outputs use the corresponding RGS score column supported by the command. |
| `--lambda-penalty` | float | 1.0 | No | Penalty coefficient λ for the discordance term. Higher values more strongly penalize RBPs where RAS and RGS disagree. |
| `--q-hi-ras` | float | 0.99 | No | Upper quantile for robust min-max normalization of RAS scores. |
| `--q-hi-rgs` | float | 0.99 | No | Upper quantile for robust min-max normalization of RGS scores. |
| `--celltypes-csv` | str | — | Cond. | CSV with columns `cell_id`, `cell_type`. **Required when `--mode ct`.** |
| `--min_cells_pert_ct` | int | 25 | No | Minimum number of cells required for a cell type to be retained in `ct` mode. |
| `--do-fdr` | int | 1 | No | Apply BH-FDR correction in `ct` mode. `1` = apply (default); `0` = skip. |

### Outputs

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
MOTIF_LINKS=rbp_motif_annotation.csv
RANKINGS_3UTR=hg38_3UTR_v10_rankings.feather
MAGMA_BIN=/tools/magma
GWAS_RAW=scz_gwas.genes.raw
GENE_LOC=NCBI38.gene.loc

# ── Step 1: Downsample real cells for the RAS branch (optional) ─────────
scRBP getSketch \
    --input $MATRIX --output pbmc_sketch.feather \
    --n_cells 50000 --celltype_col celltype --min_cells_per_type 500

# ── Step 2: Build within-cell-type metacells for the GRN branch (optional)
# Use the annotated .h5ad input so that `adata.obs[celltype]` is available.
scRBP getMetacell \
    --input             $MATRIX \
    --output            pbmc_metacell.feather \
    --metacell_size     10 \
    --method            knn \
    --within_celltype \
    --celltype_col      celltype \
    --min_metacell_size 5
# `getGRN` uses the metacell matrix, whereas `ras` uses the sketch of real
# cells together with the cell-to-cell-type mapping emitted by `getSketch`.

# ── Step 3: Infer GRN (30 seeds; use metacells if produced) ─────────────
for SEED in $(seq 1 30); do
  scRBP getGRN \
      --matrix pbmc_metacell.feather --rbp_list $RBP_LIST \
      --output grn/grn_seed${SEED} --n_workers 20 --seed $SEED
done

# ── Step 4: Merge GRN runs ─────────────────────────────────────────────
scRBP getMerge_GRN \
    --pattern "grn/grn_seed*_scRBP_gene_GRNs.tsv" --output grn_consensus.tsv \
    --n_present 10 --present_rate 0.3

# ── Step 5: Extract modules ────────────────────────────────────────────
scRBP getModule \
    --input grn_consensus.tsv --output_merged modules.tsv \
    --top_n_list "5,10,50" --target_top_n "50" --percentile "0.75,0.90"

# ── Step 6: Prune with motifs ──────────────────────────────────────────
scRBP getPrune \
    --rbp_targets modules.tsv --motif_rbp_links $MOTIF_LINKS \
    --motif_target_ranks $RANKINGS_3UTR --save_dir pruned/ \
    --rank_threshold 1500 --nes_threshold 3.0 --n_jobs 16

# ── Step 7: Build GMT regulons ─────────────────────────────────────────
scRBP getRegulon \
    --input pruned/ctx_scores_merged.csv \
    --out-symbol regulons_symbol.gmt \
    --out-entrez regulons_entrez.gmt \
    --map-custom $GENE_LOC --min_genes 5

# ── Step 8: Merge transcript regions (if analyses were run per region) ─
scRBP mergeRegulons \
    --base_dir results/ --input regulons_symbol.gmt \
    --output regulons_all_symbol.gmt --recursive \
    --append_region_to_setname

scRBP mergeRegulons \
    --base_dir results/ --input regulons_entrez.gmt \
    --output regulons_all_entrez.gmt --recursive \
    --append_region_to_setname

# ── Step 9: Compute RAS (cell type mode; RAS uses real single cells) ──
scRBP ras \
    --mode ct --matrix pbmc_sketch.feather \
    --regulons regulons_all_symbol.gmt --out ras_out/scz \
    --celltypes-csv pbmc_sketch_cell_to_celltype.csv \
    --emit-expr-stats

# ── Step 10: Common-variant RGS via MAGMA (cell type mode) ────────────
scRBP rgs \
    --mode ct --magma $MAGMA_BIN --genes-raw $GWAS_RAW \
    --sets regulons_all_entrez.gmt --id-type entrez \
    --out rgs_out/scz --gene-loc $GENE_LOC \
    --n-null 1000 --expr-stats ras_out/scz_expr_stats.tsv

# ── Step 11: (Optional) Rare-variant RGS ──────────────────────────────
# Substitute your rare-variant summary (TADA / SAIGE-GENE+ / burden ...)
scRBP rgs_rare \
    --mode          ct \
    --rare-summary  scz_rare_gene_scores.tsv \
    --rare-gene-col gene_symbol \
    --score-mode    logbf --logbf-col logBF \
    --sets          regulons_all_symbol.gmt --id-type symbol \
    --gene-loc      $GENE_LOC \
    --cds-length    gene_union_cds.tsv --cds-col union_cds_length --cds-scale raw \
    --expr-stats    ras_out/scz_expr_stats.tsv \
    --n-null 1000 --q-bins 10 \
    --out rgs_rare_out/scz_ct_rare

# ── Step 12: Compute TRS (common- or rare-variant) ────────────────────
# Common-variant TRS
scRBP trs \
    --mode ct --ras ras_out/scz.loom \
    --rgs-csv rgs_out/scz_real.csv \
    --out-prefix trs_out/scz_trs_common \
    --celltypes-csv pbmc_sketch_cell_to_celltype.csv

# Rare-variant TRS — same TRS command, just point at the rare RGS CSV
scRBP trs \
    --mode ct --ras ras_out/scz.loom \
    --rgs-csv rgs_rare_out/scz_ct_rare.gsa_RGS.csv \
    --out-prefix trs_out/scz_trs_rare \
    --celltypes-csv pbmc_sketch_cell_to_celltype.csv
```

---

## Command Reference

| Step | Command | Key Inputs | Key Output |
|------|---------|-----------|------------|
| 1 | `scRBP getSketch` | `.h5ad` / `.feather` | Downsampled cells |
| 2 | `scRBP getMetacell` *(optional, GRN branch)* | `.h5ad` (with cell-type column) / `.feather` | Metacell matrix (gene × metacell) + metacell→cell type map |
| 3 | `scRBP getGRN` | Expression / metacell matrix, RBP list | `*_scRBP_gene_GRNs.tsv` or `*_scRBP_isoform_GRNs.tsv` |
| 4 | `scRBP getMerge_GRN` | Multiple GRN TSV files (glob) | Consensus GRN TSV |
| 5 | `scRBP getModule` | Consensus GRN TSV | Modules TSV |
| 6 | `scRBP getPrune` | Modules TSV and motif-ranking resources | Pruned motif-enrichment tables |
| 7 | `scRBP getRegulon` | Pruned enrichment table | Gene-symbol and Entrez-ID regulon GMT files |
| 8 | `scRBP mergeRegulons` | Multiple GMT files | Merged GMT |
| 9 | `scRBP ras` (`--mode sc\|ct`) | Expression matrix and gene-symbol GMT | Single-cell or cell type–level RAS matrix; optional expression-statistics TSV |
| 10 | `scRBP rgs` (`--mode sc\|ct`) | MAGMA `.genes.raw`, regulon GMT, and optional expression statistics | Common-variant RGS and matched-null results |
| 11 | `scRBP rgs_rare` (`--mode sc\|ct`) | Gene-level rare-variant statistics, regulon GMT, CDS covariate, and optional expression statistics | Rare-variant RGS; REAL/NULL regulons in `ct` mode |
| 12 | `scRBP trs` (`--mode sc\|ct`) | RAS matrix, RGS CSV (common **or** rare) | TRS scores CSV |

Run `scRBP <command> --help` to view the complete command-specific argument list.

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

> Ma Y. *et al.* **Single-cell maps of RNA-binding protein networks reveal post-transcriptional architecture across development and disease.** (2026)

---



## License

MIT License. See [LICENSE](LICENSE) for details.

---

