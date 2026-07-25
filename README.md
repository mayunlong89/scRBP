<table>
  <tr>
    <td>
      <img src="https://raw.githubusercontent.com/mayunlong89/scRBP/main/Examples/scRBP_logo.png" width="350">
    </td>
    <td>
      <h1>scRBP</h1>
      <p>A scalable framework for inferring RNA-binding protein regulons from single-cell data</p>
    </td>
  </tr>
</table>


![pypi](https://img.shields.io/badge/pypi-0.1.4.1-green)
![python](https://img.shields.io/badge/python-3.9--3.11-blue)
![license](https://img.shields.io/badge/license-MIT-yellow)

[scRBP](https://github.com/mayunlong89/scRBP) (single-cell RNA Binding Protein regulon inference) is a computational framework that reconstructs RNA-binding protein (RBP)-mediated regulatory networks (regulons) from single-cell transcriptomes. It integrates motif-binding evidence with either common-variant GWAS signals or gene-level rare-variant association statistics to prioritize trait-associated RBP programs across diverse cellular contexts.

## scRBP workflow
![Workflow](https://github.com/mayunlong89/scRBP/blob/main/Examples/Figure_1_copy.png)

> **Version 0.1.4.1** adds `getMetacell` for mini-metacell construction, extends `ras` with expression-statistics output for matched-null analyses, and introduces `rgs_rare`, a rare-variant model parallel to the MAGMA-based common-variant `rgs` workflow.

---

## Installation

We recommend installing [scRBP](https://pypi.org/project/scRBP/) via PyPI:

```bash
pip install scRBP
```

Or install from source (development):

```bash
git clone https://github.com/mayunlong89/scRBP.git
cd scRBP/scRBP_package
pip install -e .
```

> **Requirements**: Python 3.9, 3.10, or 3.11. Python 3.12+ is not yet supported by `arboreto`.

### Install via conda (recommended for HPC / cluster)

```bash
git clone https://github.com/mayunlong89/scRBP.git
cd scRBP/scRBP_package
conda env create -f environment.yml
conda activate scrbp
pip install -e .
```
> If `CXXABI_1.3.15 not found`, try using conda's bundled libstdc++.
> 
> ```bash
> export LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6
> 
> scRBP --help
> ```


---

## How to run scRBP

### Pipeline overview

```text
Raw single-cell data (.h5ad / .feather)
          │
[Step 1]  scRBP getSketch        ── Stratified GeoSketch cell downsampling (Optional)
          │
[Step 2]  scRBP getMetacell      ── Aggregate similar cells into mini-metacells
          │                          (Optional; recommended for sparse GRN inference)
[Step 3]  scRBP getGRN           ── GRNBoost2/GENIE3 RBP→Gene/Isoform inference (30 seeds)
[Step 4]  scRBP getMerge_GRN     ── Merge N-seed GRNs → consensus network
[Step 5]  scRBP getModule        ── Extract regulon candidates (Top-N / percentile)
[Step 6]  scRBP getPrune         ── Motif-enrichment pruning via ctxcore
[Step 7]  scRBP getRegulon       ── Export pruned regulons to GMT format
[Step 8]  scRBP mergeRegulons    ── Merge region-specific GMT files (3'UTR/5'UTR/CDS/Introns)
[Step 9]  scRBP ras              ── Regulon Activity Score (AUCell) --mode sc | --mode ct
          │
          ├── Common-variant model
[Step 10] scRBP rgs              ── MAGMA-based common-variant RGS
          │
          └── Rare-variant model
[Step 11] scRBP rgs_rare         ── Rare-variant RGS from gene-level association statistics
          │
[Step 12] scRBP trs              ── Trait Relevance Score using common- or rare-variant RGS
```

`getSketch` and `getMetacell` serve complementary purposes. `getSketch` retains a diversity-preserving subset of real cells and is suitable for downstream RAS analysis, whereas `getMetacell` aggregates transcriptionally similar cells into denser pseudo-profiles for the GRN inference branch.

### Step 1 — Downsample cells with GeoSketch (Optional)

```bash
scRBP getSketch \
    --input  PBMC_healthy.h5ad \
    --output PBMC_sketch_50K.feather \
    --n_cells 50000 \
    --celltype_col celltype \
    --min_cells_per_type 500 \
    --n_pca 100 \
    --seed 42
```

### Step 2 — Construct mini-metacells (Optional)

`getMetacell` aggregates transcriptionally similar cells into denser mini-metacells to reduce single-cell sparsity and strengthen co-expression signals for downstream GRN inference. Metacells are constructed within each annotated cell type by default to avoid pooling across heterogeneous cell populations.

Cell similarity is calculated from a library-size-normalized, log1p-transformed PCA embedding, while metacell expression profiles are generated from the original input values using sum aggregation by default.

```bash
scRBP getMetacell \
    --input             PBMC_healthy.h5ad \
    --output            PBMC_metacell.feather \
    --metacell_size     10 \
    --method            knn \
    --within_celltype \
    --celltype_col      celltype \
    --min_metacell_size 5
```

> `getMetacell` is recommended for the GRN branch (`getGRN → getModule → getPrune`) when dropout is substantial. For RAS analysis, use the original or sketched single-cell matrix to preserve single-cell and cell-type resolution.

### Step 3 — Infer gene regulatory networks (GRN)

We recommend running `getGRN` with 30 random seeds for robustness:

```bash
#Gene mode
for SEED in $(seq 1 30); do
  scRBP getGRN \
      --matrix    PBMC_sketch_50K.feather \
      --rbp_list  human_RBP_list.txt \
      --output    grn_seed${SEED} \
      --mode      gene \
      --method    grnboost2 \
      --n_workers 20 \
      --correlation True \
      --seed      ${SEED}
done
# Output: grn_seed1_scRBP_gene_GRNs.tsv, grn_seed2_scRBP_gene_GRNs.tsv, ...

#Isoform mode
for SEED in $(seq 1 30);  do
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

### Step 4 — Merge GRN seeds into a consensus network

```bash
scRBP getMerge_GRN \
    --pattern "grn/grn_seed*_scRBP_gene_GRNs.tsv" \
    --output  grn_consensus.tsv \
    --n_present 10 \
    --present_rate 0.3
```

### Step 5 — Extract regulon candidate modules

```bash
scRBP getModule \
    --input              grn_consensus.tsv \
    --output_merged      modules.tsv \
    --importance_threshold 0.005 \
    --top_n_list         "5,10,50" \
    --target_top_n       "50" \
    --percentile         "0.75,0.9"
```

### Step 6 — Prune with motif-binding evidence

```bash
scRBP getPrune \
    --rbp_targets        modules.tsv \
    --motif_rbp_links    motif_to_rbp_hg38_v10.csv \
    --motif_target_ranks hg38_3UTR_rankings.feather \
    --save_dir           pruned_results/ \
    --rank_threshold     1500 \
    --nes_threshold      3.0 \
    --n_jobs             16
```
> This step requires [motif-based target gene/isoform ranking.feather](https://resources.gandallab.org/scRBP/).

### Step 7 — Export regulons to GMT format

```bash
scRBP getRegulon \
    --input      pruned_results/ctx_scores_merged.csv \
    --out-symbol regulons_symbol.gmt \
    --out-entrez regulons_entrez.gmt \
    --map-custom NCBI38.gene.loc \
    --min_genes  5
```

### Step 8 — Merge region-specific GMT files

```bash
scRBP mergeRegulons \
    --base_dir /data/scRBP_results/ \
    --input    regulons_symbol.gmt \
    --output   regulons_all_regions.gmt \
    --recursive \
    --append_region_to_setname
```

### Step 9 — Compute Regulon Activity Scores (RAS)

`ras` quantifies regulon activity using AUCell and supports two modes: single-cell (`--mode sc`) and cell-type level (`--mode ct`). In cell-type mode, scRBP summarizes single-cell regulon activity and calculates cell-type-specific regulon scores using a Jensen–Shannon-divergence-based framework.

The `--emit-expr-stats` option writes per-gene mean expression and detection-rate statistics. These statistics can be reused by `rgs` and `rgs_rare` to construct expression-matched null regulons in cell-type mode.

```bash
# Single-cell mode
scRBP ras \
    --mode      sc \
    --matrix    PBMC_sketch_50K.feather \
    --regulons  regulons_symbol.gmt \
    --out       ras_sc_output/ \
    --n_workers 8 \
    --emit-expr-stats

# Cell-type mode
scRBP ras \
    --mode          ct \
    --matrix        PBMC_sketch_50K.feather \
    --regulons      regulons_symbol.gmt \
    --out           ras_ct_output/ \
    --celltypes-csv cell_to_celltype.csv \
    --n_workers     8 \
    --emit-expr-stats
```

### Step 10 — Common-variant Regulon Genetic Association Score (RGS)

`rgs` implements the **common-variant model** by applying MAGMA competitive gene-set analysis to each RBP regulon. It tests whether genes within a regulon are more strongly associated with a GWAS trait than other genes in the genome.

The command supports single-cell (`--mode sc`) and cell-type (`--mode ct`) workflows. In cell-type mode, matched null regulons can be generated while controlling for gene-set size, MAGMA gene properties, mean expression, and detection rate.

```bash
# Single-cell mode
scRBP rgs \
    --mode      sc \
    --magma     ~/tools/magma/magma \
    --genes-raw scz_gwas.genes.raw \
    --sets      regulons_entrez.gmt \
    --id-type   entrez \
    --out       rgs_output/scz_sc_rgs \
    --n-null    1000

# Cell-type mode
scRBP rgs \
    --mode       ct \
    --magma      ~/tools/magma/magma \
    --genes-raw  scz_gwas.genes.raw \
    --sets       regulons_entrez.gmt \
    --id-type    entrez \
    --out        rgs_output/scz_ct_rgs \
    --gene-loc   NCBI38.gene.loc \
    --n-null     1000 \
    --expr-stats ras_ct_output/expr_stats.tsv \
    --threads    16
```

### Step 11 — Rare-variant Regulon Genetic Association Score (`rgs_rare`)

`rgs_rare` implements a **rare-variant model** parallel to the common-variant `rgs` workflow. It accepts gene-level rare-variant association statistics generated by external methods such as TADA/extTADA, SCHEMA-style burden tests, SAIGE-GENE+, REGENIE, or STAAR-O.

Depending on the input, gene-level evidence can be represented as `-log10(P)`, log Bayes factors (`logBF`), or a user-provided score. scRBP then fits a competitive regulon-level regression while adjusting for coding opportunity using union CDS length. In `--mode ct`, expression statistics from `ras --emit-expr-stats` are used only for matched-null construction.

> `rgs_rare` does not perform primary rare-variant association testing. It uses gene-level rare-variant summary statistics generated upstream.

```bash
# TADA/extTADA-style logBF input
scRBP rgs_rare \
    --mode          ct \
    --rare-summary  asd_tada_gene_scores.tsv \
    --rare-gene-col gene_symbol \
    --rare-id-type  symbol \
    --score-mode    logbf \
    --logbf-col     logBF \
    --sets          regulons_symbol.gmt \
    --id-type       symbol \
    --gene-loc      NCBI38.gene.loc \
    --cds-length    gene_union_cds.tsv \
    --cds-col       union_cds_length \
    --expr-stats    ras_ct_output/expr_stats.tsv \
    --n-null        1000 \
    --q-bins        10 \
    --out           rgs_rare_output/asd_ct_rare
```

For frequentist burden-test results, use `--score-mode pvalue` together with the relevant P-value column:

```bash
scRBP rgs_rare \
    --mode          sc \
    --rare-summary  saige_gene_burden.tsv \
    --rare-gene-col gene \
    --score-mode    pvalue \
    --p-col         P \
    --sets          regulons_symbol.gmt \
    --id-type       symbol \
    --out           rgs_rare_output/scz_sc_rare
```

The resulting `<out>.gsa_RGS.csv` file can be supplied directly to `scRBP trs` to calculate a rare-variant TRS.

### Step 12 — Compute Trait Relevance Score (TRS)

`trs` integrates RAS with either common-variant RGS results from `rgs` or rare-variant RGS results from `rgs_rare`. It supports `--mode sc` and `--mode ct`:

```
TRS = norm(RAS) + norm(RGS) − λ × |norm(RAS) − norm(RGS)|
```

RBPs with high TRS are **both** highly active in the corresponding cellular context **and** genetically linked to the trait. Using `rgs` produces a common-variant TRS, whereas using `rgs_rare` produces a rare-variant TRS.

```bash
# Single-cell mode
scRBP trs \
    --mode         sc \
    --ras          ras_sc_output/aucell_sc.csv \
    --rgs-csv      rgs_output/scz_sc_rgs_real.csv \
    --out-prefix   trs_output/scz_sc_trs \
    --rgs-score    mlog10p \
    --lambda-penalty 1.0

# Cell-type mode
scRBP trs \
    --mode             ct \
    --ras              ras_ct_output/aucell_ct.csv \
    --rgs-csv          rgs_output/scz_ct_rgs_real.csv \
    --out-prefix       trs_output/scz_ct_trs \
    --rgs-score        mlog10p \
    --lambda-penalty   1.0 \
    --celltypes-csv    cell_to_celltype.csv \
    --do-fdr           1

# Rare-variant TRS
scRBP trs \
    --mode             ct \
    --ras              ras_ct_output/aucell_ct.csv \
    --rgs-csv          rgs_rare_output/asd_ct_rare.gsa_RGS.csv \
    --out-prefix       trs_output/asd_ct_trs_rare \
    --lambda-penalty   1.0 \
    --celltypes-csv    cell_to_celltype.csv \
    --do-fdr           1
```

---

## Command Reference

| Step | Command | Description |
|------|---------|-------------|
| 1 | `scRBP getSketch` | Stratified GeoSketch downsampling of large single-cell datasets |
| 2 | `scRBP getMetacell` | Construct mini-metacells for robust GRN inference |
| 3 | `scRBP getGRN` | RBP→Gene/Isoform network inference (GRNBoost2 or GENIE3) |
| 4 | `scRBP getMerge_GRN` | Merge multi-seed GRN runs into a consensus network |
| 5 | `scRBP getModule` | Extract regulon candidate modules (Top-N / percentile) |
| 6 | `scRBP getPrune` | Rank-based motif-enrichment pruning via ctxcore |
| 7 | `scRBP getRegulon` | Export pruned regulons to GMT format (symbol + Entrez) |
| 8 | `scRBP mergeRegulons` | Merge region-specific GMT files across transcript regions |
| 9 | `scRBP ras` | Regulon Activity Score (RAS) via AUCell — `--mode sc\|ct` |
| 10 | `scRBP rgs` | Common-variant RGS via MAGMA — `--mode sc\|ct` |
| 11 | `scRBP rgs_rare` | Rare-variant RGS from gene-level association statistics — `--mode sc\|ct` |
| 12 | `scRBP trs` | Trait Relevance Score integrating RAS with common- or rare-variant RGS |

Use `scRBP <command> --help` for full parameter details.

---

## Example input format

### 1) Single-cell expression data

Accepted formats: `.h5ad` (AnnData/Seurat-converted), `.feather` (gene × cell matrix), `.loom`, `.csv`.

```bash
# Convert Seurat RDS to h5ad in R
library(SeuratDisk)
SaveH5Seurat(seurat_obj, filename = "PBMC.h5seurat")
Convert("PBMC.h5seurat", dest = "h5ad")
```

### 2) RBP list

Plain-text file with one RBP gene symbol per line:

```
ELAVL1
FUS
HNRNPA1
HNRNPC
IGF2BP1
...
```

### 3) GWAS gene-level association results (MAGMA output)

```
GENE       CHR      START       STOP  NSNPS  NPARAM       N        ZSTAT            P
148398       1     854993     884961     76      20  482730       0.7726      0.21988
26155        1     874583     899679     58      13  482730       0.4058      0.34244
339451       1     890967     906099     34       8  482730      0.70319      0.24097
...
```

### 4) Rare-variant gene-level association results

`rgs_rare` accepts CSV/TSV files containing one row per gene. The required score column depends on the selected `--score-mode`.

```text
gene_symbol    P           logBF    union_cds_length
DDX3X          2.1e-08     6.82     5019
SCN2A          4.7e-07     5.44     6009
CHD8           1.3e-06     4.91     7749
...
```

Use `--score-mode pvalue` with a P-value column, `--score-mode logbf` with a logBF column, or `--score-mode direct` with a precomputed gene-level score.

---

## Generate MAGMA results

```bash
# MAGMA annotation (10 kb window around TSS)
magma \
    --snp-loc  GWAS_summary.hg38.location \
    --annotate window=25,25 \
    --gene-loc NCBI38.gene.loc \
    --out      GWAS_SNP_Gene_annotation

# Gene-based association analysis
magma \
    --bfile  1000G_data/g1000_eur \
    --pval   GWAS_summary.results_Pval N=482730 \
    --gene-annot GWAS_SNP_Gene_annotation.genes.annot \
    --out    GWAS_gene_analysis
# Output: GWAS_gene_analysis.genes.raw (used as input to scRBP rgs)
```

### Install MAGMA (required for the common-variant `rgs` model only)

MAGMA is a standalone binary not available on PyPI. Download from [CNCR](https://cncr.nl/research/magma) and make it executable:

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
For more detailed MAGMA usage, see [GWASTutorial](https://cloufield.github.io/GWASTutorial/09_Gene_based_analysis/).

---

## More information
For more details, please refer to [scRBP protocols](https://mayunlong89.github.io/scRBP.github.io/), or [scRBP description](https://pypi.org/project/scRBP), [scRBP.github.io](https://github.com/mayunlong89/scRBP.github.io), or [scRBP reproduce](https://github.com/mayunlong89/scRBP_reproduce).

## scRBP motif collection
Please refer to the [Gandal Lab resource](https://resources.gandallab.org/scRBP/). We provide 616 RBPs and 20,746 corresponding motifs with quality metrics [here](https://resources.gandallab.org/scRBP/motifs/). 


## Citation

If you use scRBP in your research, please cite:

> Ma Y. *et al.* ***Single-cell maps of RNA-binding protein networks reveal post-transcriptional architecture across development and disease***. (2026)

---

