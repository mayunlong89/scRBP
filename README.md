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


![pypi](https://img.shields.io/badge/pypi-0.1.2-green)
![python](https://img.shields.io/badge/python-3.9--3.11-blue)
![license](https://img.shields.io/badge/license-MIT-yellow)

[scRBP](https://github.com/mayunlong89/scRBP) (single-cell RNA Binding Protein regulon inference) is a computational framework that reconstructs RNA-binding protein (RBP)-mediated regulatory networks (regulons) from single-cell transcriptomes, integrating motif-binding evidence and polygenic signals from genome-wide association studies (GWAS) to prioritize trait-associated RBP programs across diverse cellular contexts.

## scRBP workflow
![Workflow](https://github.com/mayunlong89/scRBP/blob/main/Examples/Figure_1.png)


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

```
Raw single-cell data (.h5ad / .feather)
          │
[Step 1]  scRBP getSketch        ── Stratified GeoSketch cell downsampling (Optional)
[Step 2]  scRBP getGRN           ── GRNBoost2/GENIE3 RBP→Gene/Isoform inference (30 seeds)
[Step 3]  scRBP getMerge_GRN     ── Merge N-seed GRNs → consensus network
[Step 4]  scRBP getModule        ── Extract regulon candidates (Top-N / percentile)
[Step 5]  scRBP getPrune         ── Motif-enrichment pruning via ctxcore
[Step 6]  scRBP getRegulon       ── Export pruned regulons to GMT format
[Step 7]  scRBP mergeRegulons    ── Merge region-specific GMT files (3'UTR/5'UTR/CDS/Introns)
[Step 8]  scRBP ras              ── Regulon Activity Score (AUCell)  --mode sc | --mode ct
[Step 9]  scRBP rgs              ── Regulon Gene-Set analysis (MAGMA) --mode sc | --mode ct
[Step 10] scRBP trs              ── Trait Relevance Score (RAS + RGS) --mode sc | --mode ct
```

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

### Step 2 — Infer gene regulatory networks (GRN)

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

### Step 3 — Merge GRN seeds into a consensus network

```bash
scRBP getMerge_GRN \
    --pattern "grn/grn_seed*_scRBP_gene_GRNs.tsv" \
    --output  grn_consensus.tsv \
    --n_present 15 \
    --present_rate 0.5
```

### Step 4 — Extract regulon candidate modules

```bash
scRBP getModule \
    --input              grn_consensus.tsv \
    --output_merged      modules.tsv \
    --importance_threshold 0.005 \
    --top_n_list         "5,10,50" \
    --target_top_n       "50" \
    --percentile         "0.75,0.9"
```

### Step 5 — Prune with motif-binding evidence

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
> This step need to download [motif-based target gene/isoform ranking.feather](https://zenodo.org/records/19256061).

### Step 6 — Export regulons to GMT format

```bash
scRBP getRegulon \
    --input      pruned_results/ctx_scores_merged.csv \
    --out-symbol regulons_symbol.gmt \
    --out-entrez regulons_entrez.gmt \
    --map-custom NCBI38.gene.loc \
    --min_genes  5
```

### Step 7 — Merge region-specific GMT files

```bash
scRBP mergeRegulons \
    --base_dir /data/scRBP_results/ \
    --input    regulons_symbol.gmt \
    --output   regulons_all_regions.gmt \
    --recursive \
    --append_region_to_setname
```

### Step 8 — Compute Regulon Activity Scores (RAS)

`ras` supports two modes: single-cell (`--mode sc`) and cell-type level (`--mode ct`).

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

### Step 9 — Regulon genetic association score (RGS)

`rgs` supports two modes: single-cell (`--mode sc`) and cell-type stratified (`--mode ct`).

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

### Step 10 — Compute Trait Relevance Score (TRS)

`trs` integrates RAS and RGS into a unified score and supports `--mode sc` and `--mode ct`:

```
TRS = norm(RAS) + norm(RGS) − λ × |norm(RAS) − norm(RGS)|
```

RBPs with high TRS are **both** highly active in the cell type **and** genetically linked to the trait.

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
```

---

## Command Reference

| Step | Command | Description |
|------|---------|-------------|
| 1 | `scRBP getSketch` | Stratified GeoSketch downsampling of large single-cell datasets |
| 2 | `scRBP getGRN` | RBP→Gene/Isoform network inference (GRNBoost2 or GENIE3) |
| 3 | `scRBP getMerge_GRN` | Merge multi-seed GRN runs into a consensus network |
| 4 | `scRBP getModule` | Extract regulon candidate modules (Top-N / percentile) |
| 5 | `scRBP getPrune` | Rank-based motif-enrichment pruning via ctxcore |
| 6 | `scRBP getRegulon` | Export pruned regulons to GMT format (symbol + Entrez) |
| 7 | `scRBP mergeRegulons` | Merge region-specific GMT files across transcript regions |
| 8 | `scRBP ras` | Regulon Activity Score (RAS) via AUCell — `--mode sc\|ct` |
| 9 | `scRBP rgs` | Regulon Genetic Association Score (RGS) via MAGMA — `--mode sc\|ct` |
| 10 | `scRBP trs` | Trait Relevance Score (TRS) integrating RAS and RGS — `--mode sc\|ct` |

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

---

## Generate MAGMA results

```bash
# MAGMA annotation (10 kb window around TSS)
magma \
    --snp-loc  GWAS_summary.hg38.location \
    --annotate window=10,10 \
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

### Install MAGMA (required for Step 9 only)

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
Please refer to the [Gandallab resource](https://resources.gandallab.org/scRBP/), we provided 616 RBPs and their corresponding 20746 motifs with quality metrics in this [resource](https://resources.gandallab.org/scRBP/motifs/). 


## Citation

If you use scRBP in your research, please cite:

> Ma Y. *et al.* ***Decoding cell-specific RNA-binding protein regulatory networks across development and disease***. (2026)

---

