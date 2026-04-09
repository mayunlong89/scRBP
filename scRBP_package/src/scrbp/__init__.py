"""
scRBP — Single-cell RNA Binding Protein regulon inference

scRBP is a computational framework that reconstructs RNA-binding protein (RBP)-mediated
regulatory networks (regulons) from single-cell transcriptomes, integrating motif-binding
evidence and polygenic signals from genome-wide association studies (GWAS) to prioritize
trait-associated RBP programs across diverse cellular contexts.

The pipeline comprises ten steps:
  1. getSketch      – Stratified GeoSketch cell downsampling
  2. getGRN         – GRNBoost2/GENIE3 RBP → Gene/Isoform network inference
  3. getMerge_GRN   – Consensus network construction from multi-seed GRN runs
  4. getModule      – Regulon candidate extraction (Top-N / percentile strategies)
  5. getPrune       – Rank-based motif-enrichment pruning of RBP-target modules
  6. getRegulon     – Export high-confidence regulons to GMT format
  7. mergeRegulons  – Merge region-specific GMT files across transcript regions
  8. ras            – Regulon Activity Score (RAS) computation via AUCell
  9. rgs            – Regulon-level Genetic Association Score (RGS) via MAGMA
  10. trs           – Trait Relevance Score (TRS) integrating RAS and RGS

Reference:
    Ma Y. et al. Decoding cell-specific RNA-binding protein regulatory networks across development and disease. (2026)
"""

# ── Suppress anndata FutureWarnings BEFORE any import touches anndata ────
# pyscenic internally does ``from anndata import read_csv, ...`` which
# triggers FutureWarning in anndata ≥ 0.10.  The filter must be installed
# here (the first module Python loads for the scrbp package) so it takes
# effect before scanpy / pyscenic pull in anndata.
import warnings as _warnings
_warnings.filterwarnings(
    "ignore",
    message=r"Importing read_\w+ from `anndata` is deprecated",
    category=FutureWarning,
)
del _warnings
# ─────────────────────────────────────────────────────────────────────────

__version__ = "0.1.3.11"
__author__ = "Yunlong Ma et al., University of Pennsylvania"
__email__ = ""
__description__ = "Single-cell RNA Binding Protein regulon inference"

from scrbp.commands.get_sketch import stratified_geosketch
from scrbp.commands.merge_regulons import scRBP_mergeRegulons

__all__ = [
    "__version__",
    "__author__",
    "stratified_geosketch",
    "scRBP_mergeRegulons",
]
