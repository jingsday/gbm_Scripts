import scanpy as sc
import anndata as an
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
gbm_pre_pairs=sc.read_h5ad('/home/jyang/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scanorama_cor_full.h5ad')
sc.pp.neighbors(gbm_pre_pairs, use_rep="X_scanorama")
sc.tl.umap(gbm_pre_pairs)
sc.tl.leiden(
    gbm_pre_pairs, key_added="leiden01", n_iterations=2, flavor="igraph", directed=False,resolution=0.1
)
sc.pl.umap(
    gbm_pre_pairs, color=["leiden01", "source"], palette=sc.pl.palettes.default_20,save='gbm_pre_pairs_scanorama_umap_leiden01.svg')
sc.pp.neighbors(gbm_pre_pairs, use_rep="X_scanorama")
sc.tl.umap(gbm_pre_pairs)
sc.tl.leiden(
    gbm_pre_pairs, key_added="leiden001", n_iterations=2, flavor="igraph", directed=False,resolution=0.01
)
sc.pl.umap(
    gbm_pre_pairs, color=["leiden001", "source"], palette=sc.pl.palettes.default_20,save='gbm_pre_pairs_scanorama_umap_leiden001.svg')
sc.tl.rank_genes_groups(gbm_pre_pairs, "leiden001", method="t-test")
sc.pl.rank_genes_groups(gbm_pre_pairs, n_genes=25, sharey=False,save='gbm_pre_pairs_scanorama_markers_leiden001.svg')
