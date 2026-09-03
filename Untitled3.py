#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pandas as pd
from scipy.io import mmread
import os
from scipy.io import mmread

os.chdir('/gpfs/scratch/hmz255/python_export/')

# Load ALL 10x files
matrix = mmread('matrix.mtx').T.tocsr()  # Cells x Genes, CSR format
barcodes = pd.read_csv('barcodes.tsv', sep='\t', header=None, index_col=0, names=['barcode'])
features = pd.read_csv('features.tsv', sep='\t', header=None, names=['gene_ids', 'gene_names', 'type'])

# Create adata
adata = sc.AnnData(X=matrix)
adata.obs_names = barcodes.index
adata.var_names = features['gene_names']

# Load metadata (preserves mixed types)
metadata = pd.read_csv('metadata.tsv', sep='\t', index_col=0)
adata.obs = metadata.reindex(adata.obs_names)

print(f"✅ LOADED: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print("Groups:", adata.obs['Group'].unique())



# In[ ]:





# In[ ]:





# In[ ]:





# In[2]:


print("\n=== FINAL STATUS ===")
print(f"Cells: {adata.n_obs:,}")
print(f"Genes: {adata.n_vars:,}")
print(f"obs columns: {list(adata.obs.columns)}")
print(f"Clusters available: {[col for col in adata.obs.columns if 'cluster' in col.lower() or 'leiden' in col.lower()]}")


# In[2]:


import scanpy as sc
sc.set_figure_params(dpi=80, facecolor='white')

# 1. Basic QC & filtering (standard scRNA-seq workflow)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 20, :]

# 2. Normalize, HVG, PCA
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(adata, svd_solver='arpack')

# 3. Neighbors + UMAP (THIS creates X_umap!)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# 4. Leiden clustering (matches your 12 clusters)
sc.tl.leiden(adata, resolution=0.4)

# 5. Your beautiful plots work NOW!
sc.pl.umap(adata, color=['leiden', 'Group', 'Patient'], ncols=1)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
             jitter=0.4, multi_panel=True, save='_qc.png')


# In[3]:


import pandas as pd

# Reload ORIGINAL metadata with clusters
meta = pd.read_csv("/gpfs/scratch/hmz255/python_export/metadata.tsv", sep="\t")
meta = meta.set_index("barcode").loc[adata.obs.index]  # Match cell order

# Add clusters back (keep QC metrics too)
for col in meta.columns:
    if col not in adata.obs.columns:
        adata.obs[col] = meta[col]

print("✅ Clusters restored!")
print("Available columns now:", adata.obs.columns.tolist())


# In[4]:


sc.set_figure_params(dpi=80, facecolor='white')

# UMAP with your original clusters
sc.pl.umap(adata, color=['seurat_clusters'], save='_clusters.png')

# Key metadata
sc.pl.umap(adata, color=['Group', 'Patient'], save='_metadata.png', ncols=1)

# QC metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
             jitter=0.4, multi_panel=True, save='_qc.png')


# In[5]:


# Check how many NAs
print("NA count:", adata.obs['keloid_normal'].isna().sum())

# Fill NAs (remove or assign based on other columns)
adata.obs['keloid_normal'] = adata.obs['keloid_normal'].fillna('normal')  # or remove
# OR remove NA cells entirely
adata = adata[~adata.obs['keloid_normal'].isna()].copy()


# In[5]:


adata.obs['keloid_normal'] = adata.obs['Group'].replace({
    'Keloid': 'keloid',
    'keloid': 'keloid',
    'Normal': 'normal',
    'normal': 'normal',
    'CASE': 'keloid',
    'CTRL': 'normal'
}).astype('category')

print(adata.obs['keloid_normal'].value_counts())


# In[8]:


# All sample identifiers
print("Unique sample names:")
print(adata.obs['orig.ident'].unique())

print("\nCells per sample:")
print(adata.obs['orig.ident'].value_counts())

print("\nSamples by Group:")
print(adata.obs.groupby(['Group', 'orig.ident']).size())


# In[51]:


print("Current columns:", adata.obs.columns.tolist())


# In[13]:


# Map ALL your samples to Keloid/Normal
sample_to_group = {
    # Keloid samples (K*, KF*, Kd*)
    'Kd1': 'Keloid', 'Kd2': 'Keloid', 'Kd3': 'Keloid', 'Kd4': 'Keloid',
    'KF1': 'Keloid', 'KF2': 'Keloid', 'KF3': 'Keloid',
    'K007CASE': 'Keloid', 'K009CASE': 'Keloid', 'K012CASE': 'Keloid', 'K013CASE': 'Keloid',

    # Normal samples (N*, NF*)
    'Nsc1': 'Normal', 'Nsc2': 'Normal', 'Nsc3': 'Normal', 
    'Nsk1': 'Normal', 'NF1': 'Normal', 'NF2': 'Normal', 'NF3': 'Normal', 'K012CTRL': 'Normal','K007CTRL':'Normal','K009CTRL':'Normal' ,'K013CTRL':'Normal'          
}

# Create the new group column
adata.obs['sample_group'] = adata.obs['orig.ident'].map(sample_to_group).fillna('Other')

print("✅ Keloid vs Normal groups created:")
print(adata.obs['sample_group'].value_counts())


# In[6]:


# Show samples labeled as 'Other'
other_cells = adata.obs[adata.obs['sample_group'] == 'Other']
print("OTHER samples (15,639 cells):")
print(other_cells['orig.ident'].value_counts())


# In[56]:


import scanpy as sc
import pandas as pd

# Calculate QC metrics first (if not done)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], inplace=True)

# Plot per SAMPLE (not sample_group)
# Plots MT/ribo % for your EXACT 22 samples
sc.pl.violin(adata, ['pct_counts_mt', 'pct_counts_ribo'], 
             groupby='orig.ident', rotation=90, size=0.5)


# In[57]:


# Remove MT- (mito) and RPS/RPL- (ribo) genes
mito_genes = adata.var_names.str.startswith('MT-')
ribo_genes = adata.var_names.str.startswith(('RPS', 'RPL'))

adata = adata[:, ~mito_genes & ~ribo_genes].copy()

print(f"✅ After filtering: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"Removed: {mito_genes.sum()} mito + {ribo_genes.sum()} ribo genes")


# In[6]:


# First, check your current clusters
print("Current clusters:")
print(adata.obs['seurat_clusters'].value_counts().sort_index())

# Remove last 3 clusters (assuming you have 0-14, keep 0-11)
# UPDATE these numbers based on your print output above
keep_clusters = list(range(12))  # [0,1,2,3,4,5,6,7,8,9,10,11]

adata_filtered = adata[adata.obs['seurat_clusters'].isin(keep_clusters)].copy()

print(f"✅ Kept {adata_filtered.n_obs} cells from {len(keep_clusters)} clusters")
print(f"Removed {adata.n_obs - adata_filtered.n_obs} cells")


# In[7]:


# New cluster distribution
print("\nFiltered clusters:")
print(adata_filtered.obs['seurat_clusters'].value_counts().sort_index())

# Plot to confirm
sc.pl.umap(adata_filtered, color='seurat_clusters', save='_filtered_12clusters.png')
sc.pl.umap(adata_filtered, color='seurat_clusters', save='_samples_filtered.png')


# In[20]:


# Convert seurat_clusters to categorical first
adata_filtered.obs['seurat_clusters'] = adata_filtered.obs['seurat_clusters'].astype('category')

# Now rank genes will work
sc.tl.rank_genes_groups(adata_filtered, 'seurat_clusters', method='wilcoxon')


# In[21]:


# Skip dendrogram - use simple dotplot
sc.pl.rank_genes_groups_dotplot(adata_filtered, n_genes=5, 
                               dendrogram=False, save='_top5_markers.png')

# OR use heatmap (more reliable)
sc.pl.rank_genes_groups_heatmap(adata_filtered, groups=list(range(12)), 
                               n_genes=3, show_gene_labels='all',
                               save='_cluster_markers.png')


# In[ ]:


print("=== DIAGNOSTICS ===")
print("Shape:", adata.shape)
print("\n.obs columns:", adata.obs.columns.tolist())
print("\n.obsm keys:", list(adata.obsm.keys()))
print("\nFirst few .obs rows:")
print(adata.obs.head())
print("\nLeiden clusters exist?", 'leiden' in adata.obs.columns)
print("Sample leiden values:", adata.obs['leiden'].value_counts() if 'leiden' in adata.obs else "No leiden")


# In[8]:


# Merge clusters 2 and 5 into "Keratinocytes"
keratinocyte_clusters = ['2', '5']
adata_filtered.obs['merged_clusters'] = adata_filtered.obs['seurat_clusters'].astype(str)

# Set keratinocytes as single group
mask = adata_filtered.obs['merged_clusters'].isin(keratinocyte_clusters)
adata_filtered.obs.loc[mask, 'merged_clusters'] = 'Keratinocytes'

# Make categorical
adata_filtered.obs['merged_clusters'] = adata_filtered.obs['merged_clusters'].astype('category')

print("New cluster distribution:")
print(adata_filtered.obs['merged_clusters'].value_counts())


# In[9]:


# UMAP with merged keratinocytes
sc.pl.umap(adata_filtered, color=['seurat_clusters', 'cell_type', 'sample_group'], 
           ncols=3, frameon=False)

# Confirm with keratinocyte markers
sc.pl.dotplot(adata_filtered, ['KRT5', 'KRT14', 'KRT1', 'KRT10'], 
              groupby='merged_clusters', save='_kc_merged.png')


# In[11]:


print("Original Seurat clusters:", adata_filtered.obs['seurat_clusters'].nunique())
print("New Leiden clusters:", adata_filtered.obs['leiden'].nunique() if 'leiden' in adata_filtered.obs else "No leiden")
print("Your cell_type names:", adata_filtered.obs['cell_type'].nunique())


# In[10]:


# FINAL FIX - create dict with STRING keys to match your category levels
cluster_names_str = {
    '0': 'Mesenchymal_fibroblasts',
    '1': 'MHCII+ endothelial cells', 
    '2': 'Keratinocytes',  # Merged 2+5
    '3': 'Proinflammatory_Fibroblasts',
    '4': 'Myofibroblasts',
    '5': 'Keratinocytes',
    '6': 'Immune cells',
    '7': 'Dendritic/Macrophage',
    '8': 'Endothelial cells',
    '9': 'Neural cells',

    '10': 'Myo-like_fibroblasts',
    '11': 'Mast cells'
}

# Map directly - NO conversion needed
adata_filtered.obs['cell_type'] = adata_filtered.obs['seurat_clusters'].astype(str).map(cluster_names_str)
adata_filtered.obs['cell_type'] = adata_filtered.obs['cell_type'].fillna(adata_filtered.obs['seurat_clusters'].astype(str)).astype('category')

print("AFTER mapping (success!):")
print(adata_filtered.obs['cell_type'].value_counts())


# In[12]:


sc.pl.umap(
    adata_filtered,
    color='cell_type',
    legend_loc='right margin',
    frameon=True,
    size=3,
    title='cell_type'
)


# In[12]:


sc.pl.umap(adata_filtered, color=['sample_group'], 
           ncols=3, frameon=False)


# In[15]:


# Check ALL your adata objects
print("Original adata obs:", 'seurat_clusters' in adata.obs)
print("Filtered adata obs:", 'seurat_clusters' in adata_filtered.obs)


# In[156]:


fibro_markers = [
    'COL1A1', 'COL1A2', 'DCN', 'LUM',
    'COL3A1', 'PDGFRA', 'THY1',
    'FAP', 'ACTA2', 'TAGLN', 'POSTN'
]

[g for g in fibro_markers if g in adata_filtered.var_names]


# In[158]:


sc.pl.dotplot(
    adata_filtered,
    fibro_markers,
    groupby='cell_type',
    standard_scale='var',
    cmap='Reds',
)


# In[ ]:





# In[ ]:





# In[ ]:





# In[13]:


print("=== DIAGNOSTICS ===")
print("Shape:", adata_filtered.shape)
print("\n.obs columns:", adata_filtered.obs.columns.tolist())
print("\n.obsm keys:", list(adata_filtered.obsm.keys()))
print("\nFirst few .obs rows:")
print(adata_filtered.obs.head())
print("\nLeiden clusters exist?", 'leiden' in adata_filtered.obs.columns)
print("Sample leiden values:", adata_filtered.obs['seurat_clusters'].value_counts() if 'leiden' in adata.obs else "No leiden")


# In[23]:


print("Available clusters in adata_filtered:")
print(adata_filtered.obs['seurat_clusters'].value_counts().sort_index())


# In[15]:


# Create paper source mapping based on your samples
paper_mapping = {
    # Paper 1 samples (e.g. Kd*, KF* samples)
    'KF1': 'Deng et al', 'KF2': 'Deng et al', 'KF3': 'Deng et al', 'NF1': 'Deng et al',
    'NF2': 'Deng et al', 'NF3': 'Deng et al',

    # Paper 2 samples (e.g. K00* CASE/CTRL)
    'K007CASE': 'Li et al', 'K007CTRL': 'Li et al',
    'K009CASE': 'Li et al', 'K009CTRL': 'Li et al', 
    'K012CASE': 'Li et al', 'K012CTRL': 'Li et al',
    'K013CASE': 'Li et al', 'K013CTRL': 'Li et al',

    # Paper 3 samples (e.g. N*, NF* normal)
    'Kd1': 'Direder et al', 'Kd2': 'Direder et al', 'Kd3': 'Direder et al',
    'Kd4': 'Direder et al', 'Nsc1': 'Direder et al', 
    'Nsc2': 'Direder et al', 'Nsc3': 'Direder et al','Nsk1':'Direder et al'
}

# Add to adata
adata.obs['paper_source'] = adata.obs['orig.ident'].map(paper_mapping).fillna('Unknown')

print("Paper distribution:")
print(adata.obs['paper_source'].value_counts())


# In[148]:


sc.pl.umap(adata_filtered, color='paper_source', 
           palette='Set1', frameon=False,
           legend_loc='right margin',
           save='_by_paper.png')


# In[151]:


import scanpy as sc

# Find top markers per labelled cell type
sc.tl.rank_genes_groups(
    adata_filtered,
    groupby='cell_type',
    method='wilcoxon'
)

# Extract top 5 markers per cell type
top5_celltype = {}

groups = adata_filtered.uns['rank_genes_groups']['names'].dtype.names

for group in groups:
    top5_celltype[group] = adata_filtered.uns['rank_genes_groups']['names'][group][:5].tolist()

print(top5_celltype)

# Dotplot with your labels
sc.pl.dotplot(
    adata_filtered,
    top5_celltype,
    groupby='cell_type',
    standard_scale='var',
    cmap='Reds',
    swap_axes=False,
    dendrogram=False
)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[11]:


fibro_clusters = ["0", "3", "4", "10"]
adata_fibro = adata_filtered[adata_filtered.obs['seurat_clusters'].astype(str).isin(fibro_clusters)].copy()
print(f"Fibroblast cells: {adata_fibro.n_obs:,}")


# In[12]:


sc.pp.scale(adata_fibro, max_value=10)
sc.tl.pca(adata_fibro, n_comps=30)


sc.pp.neighbors(adata_fibro, n_pcs=20)
sc.tl.umap(adata_fibro)
sc.tl.leiden(adata_fibro, resolution=0.2)


# In[54]:


# UMAP colored by Louvain clusters (your fibroblast subclusters)
sc.pl.umap(adata_fibro, 
           color='leiden', 
           frameon=False,
           legend_loc='right margin',
           save='_fibro_louvain_umap.png')

# Compare original vs new clusters
sc.pl.umap(adata_fibro, 
           color=['seurat_clusters', 'leiden'], 
           ncols=1,
           save='leiden.png')


# In[19]:


sc.tl.rank_genes_groups(adata_fibro, 'leiden',method='wilcoxon')
sc.pl.rank_genes_groups_dotplot(adata_fibro, n_genes=5, groupby='leiden',dendrogram=False)


# In[64]:


# Run DE
sc.tl.rank_genes_groups(
    adata_custom,
    groupby='leiden',
    method='wilcoxon'
)

# Extract results
result = adata_custom.uns['rank_genes_groups']
groups = result['names'].dtype.names

# Filter unwanted genes
top5_clean = {}

exclude_prefixes = (
    'RPL', 'RPS', 'MT-',
    'MALAT1', 'HSP'
)

for group in groups:

    genes = result['names'][group]

    filtered_genes = [
        gene for gene in genes
        if not gene.startswith(exclude_prefixes)
    ]

    top5_clean[group] = filtered_genes[:5]

print(top5_clean)

# Dotplot
sc.pl.dotplot(
    adata_custom,
    top5_clean,
    groupby='leiden',
    standard_scale='var',
    cmap='Reds',
    dendrogram=False
)


# In[55]:


import scanpy as sc
import matplotlib.pyplot as plt


# =========================================================
# Define genes
# =========================================================

adipocyte_genes = [
    'EBF2',
    'CEBPA',
    'PPARG',
    'ZNF423'
]

chondrocyte_genes = [
    'MSX2',
    'NKX2-6',
    'DLX5',
    'PITX1',
    'COMP',
    'SOX9',
    'RUNX1',
    'BMP2',
    'OGN'
]


# =========================================================
# Map Leiden clusters to fibroblast populations
# =========================================================

celltype_map = {
    '0': 'Mesenchymal',
    '1': 'Proinflammatory',
    '2': 'Myofibroblast',
    '3': 'Secretory Papillary'
}


adata_fibro_no45.obs['fibroblast_type'] = (
    adata_fibro_no45.obs['leiden']
    .astype(str)
    .map(celltype_map)
)


# =========================================================
# Define order
# =========================================================

celltype_order = [
    'Mesenchymal',
    'Proinflammatory',
    'Myofibroblast',
    'Secretory Papillary'
]


# Convert to ordered categorical
adata_fibro_no45.obs['fibroblast_type'] = (
    adata_fibro_no45.obs['fibroblast_type']
    .astype('category')
)

adata_fibro_no45.obs['fibroblast_type'] = (
    adata_fibro_no45.obs['fibroblast_type']
    .cat.reorder_categories(
        celltype_order,
        ordered=True
    )
)


# =========================================================
# Check genes exist
# =========================================================

adipocyte_present = [
    gene for gene in adipocyte_genes
    if gene in adata_fibro_no45.var_names
]

chondrocyte_present = [
    gene for gene in chondrocyte_genes
    if gene in adata_fibro_no45.var_names
]


print("Adipocyte genes found:")
print(adipocyte_present)

print("\nChondrocyte genes found:")
print(chondrocyte_present)


# =========================================================
# DOT PLOT 1 — ADIPOCYTE GENES
# =========================================================

sc.pl.dotplot(
    adata_fibro_no45,
    var_names=adipocyte_present,
    groupby='fibroblast_type',

    categories_order=celltype_order,

    standard_scale='var',

    cmap='Reds',

    dot_max=0.8,
    dot_min=0,

    title='Adipocyte genes'
)


# =========================================================
# DOT PLOT 2 — CHONDROCYTE GENES
# =========================================================

sc.pl.dotplot(
    adata_fibro_no45,
    var_names=chondrocyte_present,
    groupby='fibroblast_type',

    categories_order=celltype_order,

    standard_scale='var',

    cmap='Reds',

    dot_max=0.8,
    dot_min=0,

    title='Chondrocyte-associated genes'
)


# In[59]:


import scanpy as sc

# =========================
# Gene sets
# =========================

adipocyte_genes = [
    'EBF2',
    'CEBPA',
    'PPARG',
    'ZNF423'
]

chondrocyte_genes = [
    'MSX2',
    'NKX2-6',
    'DLX5',
    'PITX1',
    'COMP',
    'SOX9',
    'RUNX1',
    'BMP2',
    'OGN'
]


# =========================
# Check genes are present
# =========================

adipocyte_present = [
    g for g in adipocyte_genes
    if g in adata_keloid_only.var_names
]

chondrocyte_present = [
    g for g in chondrocyte_genes
    if g in adata_keloid_only.var_names
]

print("Adipocyte genes found:", adipocyte_present)
print("Chondrocyte genes found:", chondrocyte_present)


# =========================
# Adipocyte genes – Keloid
# =========================

sc.pl.dotplot(
    adata_keloid_only,
    var_names=adipocyte_present,
    groupby='leiden',
    cmap='Reds',
    dot_max=0.8,
    dot_min=0,
    title='Adipocyte genes – Keloid fibroblasts'
)


# =========================
# Chondrocyte genes – Keloid
# =========================

sc.pl.dotplot(
    adata_keloid_only,
    var_names=chondrocyte_present,
    groupby='leiden',
    cmap='Reds',
    dot_max=0.8,
    dot_min=0,
    title='Chondrocyte-associated genes – Keloid fibroblasts'
)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[56]:


import scanpy as sc

proliferation_genes = [
    "MKI67", "TOP2A", "UBE2C", "BIRC5", "CENPF",
    "CDK1", "CCNB1", "CCNB2", "PCNA", "MCM2",
    "MCM5", "TYMS"
]

# Keep only genes present in the dataset
available_genes = [
    gene for gene in proliferation_genes
    if gene in adata.var_names
]

print("Genes found:", available_genes)


# In[54]:


sc.pl.dotplot(
    adata_custom,
    var_names=available_genes,
    groupby="leiden",   # change this to your cluster column
    standard_scale="var"
)


# In[63]:


import scanpy as sc

sc.tl.rank_genes_groups(
    adata_custom,
    groupby="leiden",
    method="wilcoxon"
)


# In[26]:


# Collagen/ECM-focused fibroblast markers
collagen_markers = [
    'COL1A1', 'COL1A2', 'COL3A1',    # Main collagens
    'SPARC', 'CTHRC1', 'POSTN',      # ECM/pathogenic  
    'DCN', 'LUM', 'FBLN2'            # Matrix organization
]

sc.pl.dotplot(adata_fibro, 
              collagen_markers, 
              groupby='leiden', 
              save='fibro_collagen_markers.png')


# In[20]:


# Remove ONLY cluster 2, KEEP original clustering/UMAP
adata_fibro_no2 = adata_fibro[adata_fibro.obs['leiden'] != '5',].copy()

# Original UMAP stays exactly the same (just missing cluster 2)
sc.pl.umap(adata_fibro_no2, color='leiden', legend_loc='right margin')


# In[18]:


# Remove clusters 4 and 5
adata_fibro_no45 = adata_fibro[
    ~adata_fibro.obs['leiden'].isin(['4', '5'])
].copy()

# Plot using the original UMAP coordinates
sc.pl.umap(
    adata_fibro_no45,
    color='leiden',
    legend_loc='right margin'
)


# In[19]:


fibro_names = {
    '0': 'Mesenchymal',
    '1': 'Proinflammatory',
    '2': 'Myofibroblasts',
    '3': 'Secretory papillary',
}

adata_fibro_no45.obs['fibro_type'] = (
    adata_fibro_no45.obs['leiden']
    .astype(str)
    .map(fibro_names)
)

print(adata_fibro_no45.obs[['leiden', 'fibro_type']].head())
print(adata_fibro_no45.obs['fibro_type'].value_counts(dropna=False))


# In[105]:


sc.pl.umap(
    adata_keloid_only,
    color='fibro_type',
    legend_loc='right margin',
    size=6,          # smaller points
    alpha=0.8,       # transparency
    edges=False
)



# In[82]:


comp_fibro = adata_custom[adata_custom[:, 'COMP'].X > 0.5].copy()


# In[83]:


mean_expr = comp_fibro.X.mean(axis=0).A1 if hasattr(comp_fibro.X.mean(axis=0), 'A1') else comp_fibro.X.mean(axis=0)
top50_genes = pd.Series(mean_expr, index=comp_fibro.var_names).nlargest(50).index.tolist()


# In[84]:


sc.tl.score_genes(adata_custom, top50_genes, score_name='CompFibro')


# In[85]:


sc.pl.violin(adata_custom, ['CompFibro'], groupby='leiden', 
             rotation=45, save='_COMP_signature_violin.png')


# In[86]:


sc.pl.umap(adata_custom, color='CompFibro', cmap='RdYlBu_r', 
           save='_COMP_signature_umap.png')


# In[87]:


top50_df = pd.DataFrame({'COMP_fibroblast_signature': top50_genes})
top50_df.to_csv('comp_fibroblast__genes.csv', index=False)
print("✅ Saved comp_fibroblast_genes.csv")
print("Top 10:", top50_genes[:10])


# In[36]:


sc.tl.score_genes(adata_fibro_no2, top50_genes, score_name='COMP_signature')


# In[101]:


import scipy  # Add this line
import numpy as np
import pandas as pd

# 1. Find EBF2-high fibroblasts (top 25%)
ebf2_expr = adata_fibro_no2[:, 'EBF2'].X.mean(axis=1)
if scipy.sparse.issparse(ebf2_expr):
    ebf2_expr = ebf2_expr.A1
ebf2_expr = ebf2_expr.squeeze()

ebf2_hi_threshold = np.quantile(ebf2_expr, 0.75)
ebf2_hi = adata_fibro_no2[ebf2_expr > ebf2_hi_threshold].copy()

# 2. Get top50 signature from EBF2-high cells
mean_expr = ebf2_hi.X.mean(axis=0)
if hasattr(mean_expr, 'A1'):
    mean_expr = mean_expr.A1
top50_ebf2 = pd.Series(mean_expr, index=ebf2_hi.var_names).nlargest(50).index.tolist()

# 3. Save EBF2 signature
pd.DataFrame({'EBF2_signature': top50_ebf2}).to_csv('ebf2_fibroblastFIB_signature.csv', index=False)



# In[102]:


# Print the EBF2 fibroblast signature
print("🔥 EBF2+ adipogenic fibroblast signature (top 50):")
print(top50_ebf2)

# Top 10 only
print("\nTop 10 EBF2 signature genes:")
for i, gene in enumerate(top50_ebf2[:50], 1):
    print(f"{i:2d}. {gene}")


# In[106]:


# Score EBF2 signature on keloid fibroblasts
sc.tl.score_genes(data_normal_only, top50_ebf2, score_name='EBF2_signature')

# UMAP colored by EBF2 signature strength
sc.pl.umap(data_normal_only, 
           color='EBF2_signature', 
           cmap='RdYlBu_r', 
           vmax=0.2,
           frameon=False,
           save='_ebf2_signature_umap.png')


# In[108]:


# Violin plot of EBF2 signature across clusters

sc.pl.violin(
    adata_keloid_only,
    keys='EBF2_signature',
    groupby='leiden',
    stripplot=True,
    rotation=45
)


# In[1]:


chondrocyte_sig = [
    "SOX9", "ACAN", "COL2A1", "COL9A1", "COL9A2",
    "COL11A1", "COMP", "MATN3", "PRG4", "HAPLN1"
]

sc.tl.score_genes(adata, gene_list=chondrocyte_sig, score_name="ChondrocyteScore")
sc.pl.umap(adata, color=["ChondrocyteScore"])


# In[110]:


# keep only keloid cells
adata_keloid = adata_fibro_no2[adata_fibro_no2.obs["sample_group"] == "Keloid"].copy()

# compute top 50 expressed genes only in keloid cells
keloid_mean = adata_keloid.X.mean(axis=0)
if hasattr(keloid_mean, "A1"):
    keloid_mean = keloid_mean.A1

compfibro_keloid = pd.Series(keloid_mean, index=adata_keloid.var_names).nlargest(50).index.tolist()

print("✅ CompFibro signature from keloid cells only:")
print(compfibro_keloid[:50])
pd.Series(compfibro_keloid).to_csv("CompFibro_keloid_only.csv", index=False)


# In[57]:


import seaborn as sns
import matplotlib.pyplot as plt

# make sure leiden is categorical
adata_keloid.obs["leiden"] = adata_keloid.obs["leiden"].astype(str).astype("category")

# signature genes
sig = [g for g in adata_keloid.uns["CompFibro_keloid"] if g in adata_keloid.var_names]

# per-cell signature score
Xsig = adata_keloid[:, sig].X
if hasattr(Xsig, "toarray"):
    Xsig = Xsig.toarray()

adata_keloid.obs["CompFibro_score"] = Xsig.sum(axis=1)

# violin plot by leiden cluster
plt.figure(figsize=(7, 4))
sns.violinplot(
    data=adata_keloid.obs,
    x="leiden",
    y="CompFibro_score",
    inner="box",
    cut=0,
    scale="width"
)
plt.xlabel("Leiden cluster")
plt.ylabel("CompFibro signature score")
plt.tight_layout()
plt.show()


# In[120]:


# Score the CompFibro signature
sc.tl.score_genes(
    data_normal_only,
    compfibro_keloid,
    score_name='CompFibro_signature'
)

# UMAP of overall signature
sc.pl.umap(
    data_normal_only,
    color='CompFibro_signature',
    cmap='RdYlBu_r',
    vmin='p1',
    vmax='p99',
    frameon=False
)


# In[16]:


# keloid only
adata_keloid = adata_fibro_no2[adata_fibro_no2.obs["sample_group"] == "Keloid"].copy()

# make sure clusters are categorical
adata_keloid.obs["leiden"] = adata_keloid.obs["leiden"].astype(str).astype("category")

# rank genes per cluster
sc.tl.rank_genes_groups(
    adata_keloid,
    groupby="leiden",
    method="wilcoxon"
)

# top genes for each cluster
for cl in adata_keloid.obs["leiden"].cat.categories:
    print(f"\nCluster {cl} top genes:")
    print(sc.get.rank_genes_groups_df(adata_keloid, group=cl).head(20)[["names", "logfoldchanges", "pvals_adj"]])


# In[20]:


import numpy as np
import pandas as pd

# 1. Keloid cells only
adata_keloid = adata_keloid_only[adata_keloid_only.obs["sample_group"] == "Keloid"].copy()

# 2. Adipocyte seed genes
adipo_genes = ['PPARG', 'CEBPA', 'FABP4', 'ADIPOQ', 'EBF2', 'LPL', 'LEP']
adipo_genes = [g for g in adipo_genes if g in adata_keloid_only.var_names]

# 3. FIXED: Calculate scores on CORRECT adata_keloid (not adata_custom!)
Xseed = adata_keloid_only[:, adipo_genes].X
if hasattr(Xseed, "toarray"):
    Xseed = Xseed.toarray()

# CORRECT: Assign to adata_keloid (matching dimensions)
adata_keloid.obs["Adipo_score"] = Xseed.sum(axis=1)

# 4. Define threshold from SAME adata
threshold = np.quantile(adata_keloid.obs["Adipo_score"], 0.8)

# 5. Adipo-high population
adipo_hi = adata_keloid_only[adata_keloid.obs["Adipo_score"] >= threshold].copy()

# 6. Mean expression
Xhi = adipo_hi.X
if hasattr(Xhi, "toarray"):
    Xhi = Xhi.toarray()

gene_means = pd.Series(Xhi.mean(axis=0), index=adipo_hi.var_names).sort_values(ascending=False)
top50 = gene_means.head(50)

print(top50)
top50.to_csv("adipo_high_population_top50_genes.csv", header=["mean_expression"])


# In[32]:


# keep only keloid cells
adata_keloid = adata_keloid_only[adata_keloid_only.obs["sample_group"] == "Keloid"].copy()

# mean expression across keloid cells
keloid_mean = adata_keloid_only.X.mean(axis=0)
if hasattr(keloid_mean, "A1"):
    keloid_mean = keloid_mean.A1

# top 50 genes
top50_genes = (
    pd.Series(keloid_mean, index=adata_keloid.var_names)
    .sort_values(ascending=False)
    .head(50)
)

# print them
print(top50_genes)

# save them
top50_genes.to_csv("keloid_top50_genes.csv", header=["mean_expression"])


# In[68]:


print("Top adipocyte genes in keloid cells:")
print(gene_means)

gene_means.to_csv("keloid_adipocyte_signature_gene_means.csv", header=["mean_expression"])


# In[72]:


import numpy as np
import pandas as pd

# keloid cells only
adata_keloid = adata_fibro_no2[adata_fibro_no2.obs["sample_group"] == "Keloid"].copy()

# make sure EBF2 exists
if "EBF2" not in adata_keloid.var_names:
    raise ValueError("EBF2 is not present in adata_keloid.var_names")

# get EBF2 expression
Xebf2 = adata_keloid[:, "EBF2"].X
if hasattr(Xebf2, "toarray"):
    Xebf2 = Xebf2.toarray().flatten()
else:
    Xebf2 = np.asarray(Xebf2).flatten()

# define EBF2-high cells, e.g. top 10%
threshold = np.quantile(Xebf2, 0.90)
adata_ebf2_hi = adata_keloid[Xebf2 >= threshold].copy()

print(adata_ebf2_hi)
print(adata_ebf2_hi.obs["leiden"].value_counts())


# In[59]:


top50_ebf2_genes = top50_ebf2.index.tolist()
print(top50_ebf2_genes)


# In[21]:


# 1. Check which samples actually made it into adata_custom
print("Samples in adata_custom:")
print(adata_fibro_no45.obs['orig.ident'].unique())
print("\nCell counts:")
print(adata_fibro_no45.obs['orig.ident'].value_counts())

# 2. Compare sizes
print(f"\nOriginal adata_keloid: {adata_fibro_no45.n_obs:,} cells")
print(f"New adata_fibro_no45:     {adata_fibro_no45.n_obs:,} cells")

# 3. Check if Kd samples are REALLY missing
kd_check = adata_fibro_no45.obs['orig.ident'].str.contains('Kd', na=False)
print(f"\nKd samples in custom: {kd_check.sum()} found")

# 4. Verify clusters preserved
print("\nCluster distribution preserved?")
print("adata_fibro_no45:", adata_fibro_no45.obs['leiden'].value_counts())
print("adata_fibro_no45:", adata_fibro_no45.obs['leiden'].value_counts())


# In[25]:


# Your desired samples (keloid + these Kd samples)
my_samples = ['KF1', 'KF2', 'KF3', 'K013CASE','K009CASE','K007CASE','K012CASE','Kd1','Kd2','Kd4','Kd3','K007CTRL', 'K009CTRL', 'K012CTRL', 'K013CTRL', 'NF1', 'NF2', 'NF3', 'Nsc1', 'Nsc2', 'Nsc3']

# From ALL fibroblasts (includes Kd normals)
adata_fibro_no45 = adata_fibro_no45[adata_fibro_no45.obs['orig.ident'].isin(my_samples)].copy()

print("Samples included:")
print(adata_fibro_no45.obs['orig.ident'].value_counts().sort_values(ascending=False))
print(f"Total cells: {adata_fibro_no45.n_obs:,}")


# In[27]:


# 1. Identify which samples are "Unknown"
unknown_mask = adata_fibro_no45.obs['sample_status'] == 'Unknown'
print("Unknown samples:")
print(adata_fibro_no45.obs.loc[unknown_mask, 'orig.ident'].value_counts())
print("Total unknown cells:", unknown_mask.sum())


# In[28]:


def classify_sample_fixed(sample):
    keloid_patterns = ['KF', 'CASE','Kd']  # KF1,2,3 + K00*CASE
    normal_patterns = ['Nsc', 'Nsk', 'NF','CTRL']  # Kd*, Nsc*, Nsk*, NF*

    if any(pattern in sample for pattern in keloid_patterns):
        return 'Keloid'
    elif any(pattern in sample for pattern in normal_patterns):
        return 'Normal'
    else:
        return 'Unknown'

# Re-run on ALL cells
adata_fibro_no45.obs['sample_statu\s'] = adata_fibro_no45.obs['orig.ident'].apply(classify_sample_fixed)

print("Fixed counts:")
print(adata_fibro_no45.obs['sample_status'].value_counts())


# In[29]:


def classify_sample_fixed(sample):
    sample = str(sample)

    keloid_patterns = ['KF', 'CASE', 'Kd']
    normal_patterns = ['Nsc', 'Nsk', 'NF', 'CTRL']

    if any(pattern in sample for pattern in keloid_patterns):
        return 'Keloid'
    elif any(pattern in sample for pattern in normal_patterns):
        return 'Normal'
    else:
        return 'Unknown'


# Classify all cells
adata_fibro_no45.obs['sample_status'] = (
    adata_fibro_no45.obs['orig.ident']
    .apply(classify_sample_fixed)
)

print("Fixed counts:")
print(adata_fibro_no45.obs['sample_status'].value_counts())

print("\nSample classifications:")
print(
    adata_fibro_no45.obs[
        ['orig.ident', 'sample_status']
    ]
    .drop_duplicates()
    .sort_values('orig.ident')
    .to_string(index=False)
)


# In[56]:


sc.pl.umap(adata_custom, color='sample_status')
print("\nSamples by status:")
print(adata_custom.obs.groupby('orig.ident')['sample_status'].value_counts())


# In[49]:


# Keloid-only subset (from your custom 41K cells)
keloid_only_mask = adata_fibro_no45.obs['sample_status'] == 'Keloid'
adata_keloid_only = adata_fibro_no45[keloid_only_mask].copy()

print("New keloid-only dataset:")
print(f"Total cells: {adata_keloid_only.n_obs:,}")
print("Samples:")
print(adata_keloid_only.obs['orig.ident'].value_counts())
print("Clusters preserved:")
print(adata_keloid_only.obs['leiden'].value_counts())


# In[118]:


import pandas as pd

df = pd.read_csv("bindetect_results_all peaks.csv")


# In[50]:


# Keloid-only subset (from your custom 41K cells)
normal_only_mask = adata_fibro_no45.obs['sample_status'] == 'Normal'
data_normal_only = adata_fibro_no45[normal_only_mask].copy()

print("New keloid-only dataset:")
print(f"Total cells: {data_normal_only.n_obs:,}")
print("Samples:")
print(data_normal_only.obs['orig.ident'].value_counts())
print("Clusters preserved:")
print(data_normal_only.obs['leiden'].value_counts())


# In[31]:


sc.pl.umap(data_normal_only, 
           color=['leiden'], 
           ncols=1,
           save='leiden.png')


# In[122]:


import pandas as pd
import numpy as np
import scipy.sparse as sp

# -----------------------------
# SETTINGS
# -----------------------------
adata = adata_custom.copy()

sample_col = "orig.ident"
condition_col = "sample_status"   # use this if your Keloid/Normal labels are here
cluster_col = "leiden"

min_cells = 20

# -----------------------------
# PREPARE OBS
# -----------------------------
adata.obs[sample_col] = adata.obs[sample_col].astype(str)
adata.obs[condition_col] = adata.obs[condition_col].astype(str)
adata.obs[cluster_col] = adata.obs[cluster_col].astype(str)

# Keep only Keloid and Normal
adata = adata[
    adata.obs[condition_col].isin(["Keloid", "Normal"])
].copy()

print(adata.obs[condition_col].value_counts())

# -----------------------------
# CREATE PSEUDOBULK
# -----------------------------
all_counts = []
all_metadata = []

for cluster in sorted(adata.obs[cluster_col].unique()):

    print(f"\nProcessing cluster {cluster}")

    adata_cluster = adata[
        adata.obs[cluster_col] == cluster
    ].copy()

    for sample in adata_cluster.obs[sample_col].unique():

        sample_cells = adata_cluster[
            adata_cluster.obs[sample_col] == sample
        ]

        n_cells = sample_cells.n_obs

        if n_cells < min_cells:
            continue

        X = sample_cells.X

        if sp.issparse(X):
            summed_counts = np.asarray(X.sum(axis=0)).flatten()
        else:
            summed_counts = np.asarray(X.sum(axis=0)).flatten()

        sample_name = f"{sample}_cluster{cluster}"

        all_counts.append(
            pd.Series(
                summed_counts,
                index=adata_cluster.var_names,
                name=sample_name
            )
        )

        all_metadata.append({
            "sample": sample_name,
            "orig_sample": sample,
            "cluster": cluster,
            "condition": sample_cells.obs[condition_col].iloc[0],
            "n_cells": n_cells
        })

# -----------------------------
# FINAL MATRICES
# -----------------------------
count_df = pd.concat(all_counts, axis=1)

meta_df = pd.DataFrame(all_metadata)
meta_df = meta_df.set_index("sample")

print("Counts matrix:", count_df.shape)
print("Metadata:", meta_df.shape)

print(meta_df.groupby(["cluster", "condition"]).size())

# -----------------------------
# SAVE
# -----------------------------
count_df.to_csv("adata_custom_pseudobulk_counts.csv")
meta_df.to_csv("adata_custom_pseudobulk_metadata.csv")

print("\nSaved:")
print("adata_custom_pseudobulk_counts.csv")
print("adata_custom_pseudobulk_metadata.csv")


# In[85]:


adata_no4 = adata_custom[adata_custom.obs["leiden"] != "4"].copy()


# In[86]:


adata_no4.obs["leiden"].value_counts()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[67]:


import scanpy as sc
import pandas as pd

# Set the fibroblast-type order
fibro_order = [
    'Mesenchymal',
    'Myofibroblasts',
    'Proinflammatory',
    'Secretory papillary'
]

# Keep only categories that are actually present
present_types = [
    x for x in fibro_order
    if x in adata_fibro_no45.obs['fibro_type'].astype(str).unique()
]

adata_fibro_no45.obs['fibro_type'] = pd.Categorical(
    adata_fibro_no45.obs['fibro_type'],
    categories=present_types,
    ordered=True
)

# Colours matching the example
fibro_palette = {
    'Mesenchymal': '#377eb8',
    'Myofibroblasts': '#ff7f00',
    'Proinflammatory': '#35a875',
    'Secretory papillary': '#ef4b4b'
}

colors = [fibro_palette[x] for x in present_types]

# Plot original UMAP coordinates
sc.pl.umap(
    adata_fibro_no45,
    color='fibro_type',
    palette=colors,
    title='fibro_type',
    legend_loc='right margin',
    size=5,
    frameon=True,
    edges=False,
    show=True
)


# In[ ]:





# In[ ]:





# In[ ]:





# In[64]:


# Clean keloid-only UMAP
sc.pl.umap(adata_keloid_only, color=['leiden', 'paper_source'], ncols=1)

# Keloid cluster proportions
print("Keloid-only cluster proportions:")
print(adata_keloid.obs['leiden'].value_counts(normalize=True))


# In[47]:


# Step 1: Add paper_source to adata_custom (same mapping as before)
paper_mapping = {
    'KF1': 'Deng et al', 'KF2': 'Deng et al', 'KF3': 'Deng et al','NF1': 'Deng et al', 'NF2': 'Deng et al', 'NF3': 'Deng et al',
    'K007CASE': 'Li et al', 'K007CTRL': 'Li et al',
    'K009CASE': 'Li et al', 'K009CTRL': 'Li et al', 
    'K012CASE': 'Li et al', 'K012CTRL': 'Li et al',
    'K013CASE': 'Li et al', 'K013CTRL': 'Li et al',
    'Kd1': 'Direder et al', 'Kd2': 'Direder et al', 'Kd3': 'Direder et al',
    'Kd4': 'Direder et al', 'Nsc1': 'Direder et al', 
    'Nsc2': 'Direder et al', 'Nsc3': 'Direder et al','Nsk1':'Direder et al'
}

adata_fibro_no45.obs['paper_source'] = adata_fibro_no45.obs['orig.ident'].map(paper_mapping).fillna('Unknown')

print("✅ Paper distribution in adata_custom:")
print(adata_fibro_no45.obs['paper_source'].value_counts())


# In[65]:


paper_samples = (adata_custom.obs.groupby('paper_source')['orig.ident']
                .value_counts()
                .loc[lambda x: x > 0]  # Keep only real samples
                .sort_values(ascending=False))
print(paper_samples)


# In[35]:


# Step 1: Add paper_source to adata_custom (same mapping as before)
paper_mapping = {
    'KF1': 'Deng et al', 'KF2': 'Deng et al', 'KF3': 'Deng et al',
    'K007CASE': 'Li et al', 
    'K009CASE': 'Li et al', 
    'K012CASE': 'Li et al',
    'K013CASE': 'Li et al',
    'Kd1': 'Direder et al', 'Kd2': 'Direder et al', 'Kd3': 'Direder et al',
    'Kd4': 'Direder et al',

adata_fibro_no45.obs['paper_source'] = adata_fibro_no45.obs['orig.ident'].map(paper_mapping).fillna('Unknown')

print("✅ Paper distribution in adata_custom:")
print(adata_fibro_no45.obs['paper_source'].value_counts())


# In[36]:


paper_samples = (adata_fibro_no45.obs.groupby('paper_source')['orig.ident']
                .value_counts()
                .loc[lambda x: x > 0]  # Keep only real samples
                .sort_values(ascending=False))
print(paper_samples)


# In[32]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

sample_col = 'orig.ident'
cluster_col = 'leiden'
paper_col = 'paper_source'
condition_col = 'sample_status'

# Cluster names
celltype_map = {
    '0': 'Mesenchymal',
    '1': 'Proinflammatory',
    '2': 'Myofibroblast',
    '3': 'Secretory Papillary',
    '4': 'Proliferating'
}

df = adata_fibro_no45.obs[
    [sample_col, cluster_col, paper_col, condition_col]
].dropna().copy()

# Convert categorical columns to strings
df[sample_col] = df[sample_col].astype(str)
df[cluster_col] = df[cluster_col].astype(str)
df[paper_col] = df[paper_col].astype(str)
df[condition_col] = df[condition_col].astype(str)

# Add cell-type names
df['cell_type'] = df[cluster_col].map(celltype_map)

# Count cells per sample / condition / cell type
counts = (
    df.groupby([sample_col, condition_col, 'cell_type'])
    .size()
    .reset_index(name='n_cells')
)

# Total cells per sample / condition
totals = (
    df.groupby([sample_col, condition_col])
    .size()
    .reset_index(name='total_cells')
)

# Calculate percentage
plot_df = counts.merge(totals, on=[sample_col, condition_col])
plot_df['percentage'] = plot_df['n_cells'] / plot_df['total_cells'] * 100

# Add paper source
sample_info = df[[sample_col, paper_col]].drop_duplicates()
plot_df = plot_df.merge(sample_info, on=sample_col)

# Convert again after merge
plot_df['cell_type'] = plot_df['cell_type'].astype(str)
plot_df[condition_col] = plot_df[condition_col].astype(str)
plot_df[paper_col] = plot_df[paper_col].astype(str)

# Colours
paper_palette = {
    'Deng et al': '#1f77b4',
    'Li et al': '#d62728',
    'Direder et al': '#2ca02c'
}

condition_palette = {
    'Keloid': '#1f77b4',
    'Normal': '#ff7f0e'
}

# Order of cell types
celltype_order = [
    'Mesenchymal',
    'Myofibroblast',
    'Proinflammatory',
    'Secretory Papillary',
    'Proliferating'
]

celltype_to_x = {ct: i for i, ct in enumerate(celltype_order)}

# Offsets for bars
condition_offsets = {
    'Keloid': -0.2,
    'Normal': 0.2
}

# Add numeric x positions for dots
plot_df['x_pos'] = plot_df['cell_type'].map(celltype_to_x).astype(float)
plot_df['offset'] = plot_df[condition_col].map(condition_offsets).astype(float)
plot_df['x_bar'] = plot_df['x_pos'] + plot_df['offset']

# Mean percentage for bars
bar_df = (
    plot_df.groupby(['cell_type', condition_col], as_index=False)['percentage']
    .mean()
)

bar_df['x_pos'] = bar_df['cell_type'].map(celltype_to_x).astype(float)
bar_df['offset'] = bar_df[condition_col].map(condition_offsets).astype(float)
bar_df['x_bar'] = bar_df['x_pos'] + bar_df['offset']

plt.figure(figsize=(18, 6))
ax = plt.gca()

# Plot bars
for _, row in bar_df.iterrows():
    ax.bar(
        row['x_bar'],
        row['percentage'],
        width=0.35,
        color=condition_palette.get(row[condition_col], 'gray'),
        alpha=0.35
    )

# Plot sample dots
for _, row in plot_df.iterrows():
    jitter = np.random.uniform(-0.04, 0.04)

    ax.scatter(
        row['x_bar'] + jitter,
        row['percentage'],
        color=paper_palette.get(row[paper_col], 'gray'),
        edgecolor='black',
        s=60,
        zorder=10
    )

# X-axis labels
ax.set_xticks(range(len(celltype_order)))
ax.set_xticklabels(celltype_order, rotation=45, ha='right')

ax.set_ylabel('Percentage of cells per sample (%)')
ax.set_xlabel('Cell type')
ax.set_title('Cell-type composition across samples')

# Legends
condition_legend = [
    Patch(facecolor=v, alpha=0.35, label=k)
    for k, v in condition_palette.items()
]

paper_legend = [
    Line2D(
        [0], [0],
        marker='o',
        color='w',
        label=k,
        markerfacecolor=v,
        markeredgecolor='black',
        markersize=8
    )
    for k, v in paper_palette.items()
]

legend1 = ax.legend(
    handles=condition_legend,
    title='Condition',
    bbox_to_anchor=(1.02, 1),
    loc='upper left'
)

ax.add_artist(legend1)

ax.legend(
    handles=paper_legend,
    title='Paper source',
    bbox_to_anchor=(1.02, 0.65),
    loc='upper left'
)

plt.tight_layout()
plt.show()


# In[48]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# =========================================================
# 1. Define metadata columns
# =========================================================

sample_col = 'orig.ident'
cluster_col = 'leiden'
paper_col = 'paper_source'
condition_col = 'sample_status'


# =========================================================
# 2. Map Leiden clusters to fibroblast cell types
# =========================================================

celltype_map = {
    '0': 'Mesenchymal',
    '1': 'Proinflammatory',
    '2': 'Myofibroblast',
    '3': 'Secretory Papillary',
}


# =========================================================
# 3. Check that all required columns exist
# =========================================================

required_columns = [
    sample_col,
    cluster_col,
    paper_col,
    condition_col
]

missing_columns = [
    col for col in required_columns
    if col not in adata_fibro_no45.obs.columns
]

if missing_columns:
    raise KeyError(
        f"These columns are missing from adata_custom.obs: "
        f"{missing_columns}\n\n"
        f"Available columns:\n"
        f"{adata_fibro_no45.obs.columns.tolist()}"
    )


# =========================================================
# 4. Create the working dataframe
# =========================================================

df = adata_fibro_no45.obs[required_columns].copy()

# Convert values to strings
for col in required_columns:
    df[col] = df[col].astype(str)

# Remove rows with missing-like string values
df = df[
    ~df[required_columns].isin(
        ['nan', 'None', '']
    ).any(axis=1)
].copy()

# Map Leiden clusters to fibroblast cell types
df['cell_type'] = df[cluster_col].map(celltype_map)

# Remove clusters that are not included in celltype_map
df = df.dropna(subset=['cell_type']).copy()


# =========================================================
# 5. Inspect the samples and classifications
# =========================================================

print("Samples and conditions:")
print(
    df[
        [sample_col, condition_col, paper_col]
    ]
    .drop_duplicates()
    .sort_values([condition_col, sample_col])
    .to_string(index=False)
)

print("\nCell counts by type:")
print(df['cell_type'].value_counts())


# =========================================================
# 6. Count cells per sample, condition and cell type
# =========================================================

counts = (
    df.groupby(
        [sample_col, condition_col, 'cell_type'],
        observed=True
    )
    .size()
    .reset_index(name='n_cells')
)


# =========================================================
# 7. Calculate the total number of cells per sample
# =========================================================

totals = (
    df.groupby(
        [sample_col, condition_col],
        observed=True
    )
    .size()
    .reset_index(name='total_cells')
)


# =========================================================
# 8. Calculate cell-type percentage per sample
# =========================================================

plot_df = counts.merge(
    totals,
    on=[sample_col, condition_col],
    how='left'
)

plot_df['percentage'] = (
    plot_df['n_cells']
    / plot_df['total_cells']
    * 100
)


# =========================================================
# 9. Add combinations where a sample has zero cells
#    of a particular fibroblast type
# =========================================================

celltype_order = [
    'Mesenchymal',
    'Myofibroblast',
    'Proinflammatory',
    'Secretory Papillary']

sample_metadata = (
    df[
        [sample_col, condition_col, paper_col]
    ]
    .drop_duplicates(subset=[sample_col])
)

complete_index = pd.MultiIndex.from_product(
    [
        sample_metadata[sample_col].unique(),
        celltype_order
    ],
    names=[sample_col, 'cell_type']
)

complete_df = (
    complete_index
    .to_frame(index=False)
    .merge(
        sample_metadata,
        on=sample_col,
        how='left'
    )
)

plot_df = complete_df.merge(
    plot_df[
        [
            sample_col,
            condition_col,
            'cell_type',
            'n_cells',
            'total_cells',
            'percentage'
        ]
    ],
    on=[
        sample_col,
        condition_col,
        'cell_type'
    ],
    how='left'
)

# Set missing cell-type counts and percentages to zero
plot_df['n_cells'] = plot_df['n_cells'].fillna(0)
plot_df['percentage'] = plot_df['percentage'].fillna(0)

# Add total cell counts where needed
sample_totals = totals[
    [sample_col, 'total_cells']
].drop_duplicates()

plot_df = plot_df.drop(
    columns=['total_cells'],
    errors='ignore'
).merge(
    sample_totals,
    on=sample_col,
    how='left'
)


# =========================================================
# 10. Keep only Keloid and Normal conditions
# =========================================================

plot_df = plot_df[
    plot_df[condition_col].isin(
        ['Keloid', 'Normal']
    )
].copy()


# =========================================================
# 11. Define colours and plotting positions
# =========================================================

paper_palette = {
    'Deng et al': '#1f77b4',
    'Li et al': '#d62728',
    'Direder et al': '#2ca02c'
}

condition_palette = {
    'Keloid': '#1f77b4',
    'Normal': '#ff7f0e'
}

condition_order = [
    'Keloid',
    'Normal'
]

celltype_to_x = {
    cell_type: index
    for index, cell_type in enumerate(celltype_order)
}

condition_offsets = {
    'Keloid': -0.2,
    'Normal': 0.2
}

plot_df['x_pos'] = (
    plot_df['cell_type']
    .map(celltype_to_x)
    .astype(float)
)

plot_df['offset'] = (
    plot_df[condition_col]
    .map(condition_offsets)
    .astype(float)
)

plot_df['x_bar'] = (
    plot_df['x_pos']
    + plot_df['offset']
)


# =========================================================
# 12. Calculate mean and standard deviation across samples
# =========================================================

bar_df = (
    plot_df.groupby(
        ['cell_type', condition_col],
        observed=True,
        as_index=False
    )
    .agg(
        mean_percentage=('percentage', 'mean'),
        sd_percentage=('percentage', 'std'),
        n_samples=(sample_col, 'nunique')
    )
)

# SD is undefined for groups containing only one sample
bar_df['sd_percentage'] = (
    bar_df['sd_percentage']
    .fillna(0)
)

bar_df['x_pos'] = (
    bar_df['cell_type']
    .map(celltype_to_x)
    .astype(float)
)

bar_df['offset'] = (
    bar_df[condition_col]
    .map(condition_offsets)
    .astype(float)
)

bar_df['x_bar'] = (
    bar_df['x_pos']
    + bar_df['offset']
)

print("\nSummary statistics:")
print(
    bar_df[
        [
            'cell_type',
            condition_col,
            'n_samples',
            'mean_percentage',
            'sd_percentage'
        ]
    ]
    .sort_values(['cell_type', condition_col])
    .to_string(index=False)
)


# =========================================================
# 13. Plot mean bars with standard deviation error bars
# =========================================================

# Fixed random seed gives the same jitter each time
rng = np.random.default_rng(42)

fig, ax = plt.subplots(figsize=(14, 7))

# Plot mean bars and SD error bars
for _, row in bar_df.iterrows():

    ax.bar(
        x=row['x_bar'],
        height=row['mean_percentage'],
        width=0.35,
        color=condition_palette.get(
            row[condition_col],
            'gray'
        ),
        alpha=0.35,
        edgecolor='black',
        linewidth=0.8,
        yerr=row['sd_percentage'],
        capsize=5,
        error_kw={
            'elinewidth': 1.5,
            'capthick': 1.5,
            'ecolor': 'black'
        },
        zorder=2
    )


# =========================================================
# 14. Add individual sample points
# =========================================================

for _, row in plot_df.iterrows():

    jitter = rng.uniform(-0.045, 0.045)

    ax.scatter(
        row['x_bar'] + jitter,
        row['percentage'],
        color=paper_palette.get(
            row[paper_col],
            'gray'
        ),
        edgecolor='black',
        linewidth=0.7,
        s=65,
        alpha=0.9,
        zorder=10
    )


# =========================================================
# 15. Format the axes
# =========================================================

ax.set_xticks(range(len(celltype_order)))

ax.set_xticklabels(
    celltype_order,
    rotation=45,
    ha='right'
)

ax.set_ylabel(
    'Percentage of cells per sample (%)',
    fontsize=12
)

ax.set_xlabel(
    'Fibroblast cell type',
    fontsize=12
)

ax.set_title(
    'Fibroblast cell-type composition across samples\n'
    'Mean ± standard deviation',
    fontsize=14
)

ax.set_ylim(
    bottom=0
)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.grid(
    axis='y',
    linestyle='--',
    alpha=0.3,
    zorder=0
)


# =========================================================
# 16. Create legends
# =========================================================

condition_legend = [
    Patch(
        facecolor=condition_palette[condition],
        edgecolor='black',
        alpha=0.35,
        label=condition
    )
    for condition in condition_order
    if condition in plot_df[condition_col].unique()
]

paper_legend = [
    Line2D(
        [0],
        [0],
        marker='o',
        linestyle='None',
        label=paper,
        markerfacecolor=colour,
        markeredgecolor='black',
        markersize=8
    )
    for paper, colour in paper_palette.items()
    if paper in plot_df[paper_col].unique()
]

legend1 = ax.legend(
    handles=condition_legend,
    title='Condition',
    bbox_to_anchor=(1.02, 1),
    loc='upper left',
    frameon=False
)

ax.add_artist(legend1)

ax.legend(
    handles=paper_legend,
    title='Paper source',
    bbox_to_anchor=(1.02, 0.68),
    loc='upper left',
    frameon=False
)


# =========================================================
# 17. Display the figure
# =========================================================



# In[75]:


status_mapping = {
    # Deng
    'KF1': 'keloid',
    'KF2': 'keloid',
    'KF3': 'keloid',

    # Li
    'K007CASE': 'keloid',
    'K009CASE': 'keloid',
    'K012CASE': 'keloid',
    'K013CASE': 'keloid',

    # Direder
    'Kd1': 'keloid',
    'Kd2': 'keloid',
    'Kd3': 'keloid',
    'Kd4': 'keloid',

    # Add NORMAL samples here
    'N1': 'Normal',
    'N2': 'Normal',
    'N3': 'Normal',

    'NF1': 'Normal', 
    'NF2': 'Normal', 
    'NF3': 'Normal',

    'K007CTRL': 'Normal',
    'K009CTRL': 'Normal',
    'K012CTRL': 'Normal',
    'K013CTRL': 'Normal',

    'Nsc1': 'Normal', 
    'Nsc2': 'Normal', 
    'Nsc3': 'Normal',
    'Nsk1':'Normal'
}


# In[38]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# =========================================================
# 1. Define metadata columns
# =========================================================

sample_col = 'orig.ident'
cluster_col = 'leiden'
paper_col = 'paper_source'
condition_col = 'sample_status'


# =========================================================
# 2. Map Leiden clusters to fibroblast cell types
# =========================================================

celltype_map = {
    '0': 'Mesenchymal',
    '1': 'Proinflammatory',
    '2': 'Myofibroblast',
    '3': 'Secretory Papillary',
}


# =========================================================
# 3. Check that all required columns exist
# =========================================================

required_columns = [
    sample_col,
    cluster_col,
    paper_col,
    condition_col
]

missing_columns = [
    col
    for col in required_columns
    if col not in adata_fibro_no45.obs.columns
]

if missing_columns:
    raise KeyError(
        f"These columns are missing from adata_fibro_no45.obs: "
        f"{missing_columns}\n\n"
        f"Available columns:\n"
        f"{adata_fibro_no45.obs.columns.tolist()}"
    )


# =========================================================
# 4. Create the working dataframe
# =========================================================

df = adata_fibro_no45.obs[required_columns].copy()


# Convert values to strings
for col in required_columns:
    df[col] = df[col].astype(str)


# Remove rows with missing-like string values
df = df[
    ~df[required_columns].isin(
        ['nan', 'None', '']
    ).any(axis=1)
].copy()


# Map Leiden clusters to fibroblast cell types
df['cell_type'] = df[cluster_col].map(celltype_map)


# Remove clusters that are not included in celltype_map
df = df.dropna(
    subset=['cell_type']
).copy()


# =========================================================
# 5. Inspect samples, paper sources and classifications
# =========================================================

print("\n========================================")
print("Samples and conditions")
print("========================================")

print(
    df[
        [sample_col, condition_col, paper_col]
    ]
    .drop_duplicates()
    .sort_values(
        [condition_col, sample_col]
    )
    .to_string(index=False)
)


print("\n========================================")
print("Paper sources")
print("========================================")

print(
    df[paper_col]
    .value_counts()
)


print("\n========================================")
print("Conditions")
print("========================================")

print(
    df[condition_col]
    .value_counts()
)


print("\n========================================")
print("Cell counts by type")
print("========================================")

print(
    df['cell_type']
    .value_counts()
)


# =========================================================
# 6. Count cells per sample, condition and cell type
# =========================================================

counts = (
    df.groupby(
        [
            sample_col,
            condition_col,
            'cell_type'
        ],
        observed=True
    )
    .size()
    .reset_index(name='n_cells')
)


# =========================================================
# 7. Calculate total number of cells per sample
# =========================================================

totals = (
    df.groupby(
        [
            sample_col,
            condition_col
        ],
        observed=True
    )
    .size()
    .reset_index(name='total_cells')
)


# =========================================================
# 8. Calculate cell-type percentage per sample
# =========================================================

plot_df = counts.merge(
    totals,
    on=[
        sample_col,
        condition_col
    ],
    how='left'
)


plot_df['percentage'] = (
    plot_df['n_cells']
    / plot_df['total_cells']
    * 100
)


# =========================================================
# 9. Add combinations where a sample has zero cells
#    of a particular fibroblast type
# =========================================================

celltype_order = [
    'Mesenchymal',
    'Myofibroblast',
    'Proinflammatory',
    'Secretory Papillary'
]


sample_metadata = (
    df[
        [
            sample_col,
            condition_col,
            paper_col
        ]
    ]
    .drop_duplicates(
        subset=[sample_col]
    )
)


complete_index = pd.MultiIndex.from_product(
    [
        sample_metadata[sample_col].unique(),
        celltype_order
    ],
    names=[
        sample_col,
        'cell_type'
    ]
)


complete_df = (
    complete_index
    .to_frame(index=False)
    .merge(
        sample_metadata,
        on=sample_col,
        how='left'
    )
)


plot_df = complete_df.merge(
    plot_df[
        [
            sample_col,
            condition_col,
            'cell_type',
            'n_cells',
            'total_cells',
            'percentage'
        ]
    ],
    on=[
        sample_col,
        condition_col,
        'cell_type'
    ],
    how='left'
)


# Set missing cell-type counts and percentages to zero
plot_df['n_cells'] = (
    plot_df['n_cells']
    .fillna(0)
)


plot_df['percentage'] = (
    plot_df['percentage']
    .fillna(0)
)


# Add total cell counts where needed
sample_totals = (
    totals[
        [
            sample_col,
            'total_cells'
        ]
    ]
    .drop_duplicates()
)


plot_df = (
    plot_df
    .drop(
        columns=['total_cells'],
        errors='ignore'
    )
    .merge(
        sample_totals,
        on=sample_col,
        how='left'
    )
)


# =========================================================
# 10. Keep only Keloid and Normal conditions
# =========================================================

plot_df = plot_df[
    plot_df[condition_col].isin(
        [
            'Keloid',
            'Normal'
        ]
    )
].copy()


# =========================================================
# 11. Define colours and plotting positions
# =========================================================

# ---------------------------------------------------------
# Paper-source colours
# ---------------------------------------------------------
# NDFs has been added here as PURPLE
# ---------------------------------------------------------

paper_palette = {
    'Deng et al': '#1f77b4',      # Blue
    'Li et al': '#d62728',        # Red
    'Direder et al': '#2ca02c',   # Green
    'NDFs': '#9467bd'             # Purple
}


# ---------------------------------------------------------
# Condition colours
# ---------------------------------------------------------

condition_palette = {
    'Keloid': '#1f77b4',
    'Normal': '#ff7f0e'
}


condition_order = [
    'Keloid',
    'Normal'
]


# ---------------------------------------------------------
# X-axis positions
# ---------------------------------------------------------

celltype_to_x = {
    cell_type: index
    for index, cell_type
    in enumerate(celltype_order)
}


# ---------------------------------------------------------
# Horizontal offsets for Keloid and Normal bars
# ---------------------------------------------------------

condition_offsets = {
    'Keloid': -0.2,
    'Normal': 0.2
}


plot_df['x_pos'] = (
    plot_df['cell_type']
    .map(celltype_to_x)
    .astype(float)
)


plot_df['offset'] = (
    plot_df[condition_col]
    .map(condition_offsets)
    .astype(float)
)


plot_df['x_bar'] = (
    plot_df['x_pos']
    + plot_df['offset']
)


# =========================================================
# 12. Calculate mean and standard deviation across samples
# =========================================================

bar_df = (
    plot_df.groupby(
        [
            'cell_type',
            condition_col
        ],
        observed=True,
        as_index=False
    )
    .agg(
        mean_percentage=(
            'percentage',
            'mean'
        ),
        sd_percentage=(
            'percentage',
            'std'
        ),
        n_samples=(
            sample_col,
            'nunique'
        )
    )
)


# SD is undefined for groups containing only one sample
bar_df['sd_percentage'] = (
    bar_df['sd_percentage']
    .fillna(0)
)


bar_df['x_pos'] = (
    bar_df['cell_type']
    .map(celltype_to_x)
    .astype(float)
)


bar_df['offset'] = (
    bar_df[condition_col]
    .map(condition_offsets)
    .astype(float)
)


bar_df['x_bar'] = (
    bar_df['x_pos']
    + bar_df['offset']
)


# =========================================================
# 13. Print summary statistics
# =========================================================

print("\n========================================")
print("Summary statistics")
print("========================================")

print(
    bar_df[
        [
            'cell_type',
            condition_col,
            'n_samples',
            'mean_percentage',
            'sd_percentage'
        ]
    ]
    .sort_values(
        [
            'cell_type',
            condition_col
        ]
    )
    .to_string(index=False)
)


# =========================================================
# 14. Create figure
# =========================================================

# Fixed random seed gives the same jitter every time
rng = np.random.default_rng(42)


fig, ax = plt.subplots(
    figsize=(14, 7)
)


# =========================================================
# 15. Plot mean bars with SD error bars
# =========================================================

for _, row in bar_df.iterrows():

    ax.bar(
        x=row['x_bar'],
        height=row['mean_percentage'],
        width=0.35,

        # Condition colour
        color=condition_palette.get(
            row[condition_col],
            'gray'
        ),

        alpha=0.35,

        edgecolor='black',
        linewidth=0.8,

        # Standard deviation
        yerr=row['sd_percentage'],
        capsize=5,

        error_kw={
            'elinewidth': 1.5,
            'capthick': 1.5,
            'ecolor': 'black'
        },

        zorder=2
    )


# =========================================================
# 16. Add individual sample points
# =========================================================

for _, row in plot_df.iterrows():

    # Random horizontal jitter
    jitter = rng.uniform(
        -0.045,
        0.045
    )

    # Get paper colour
    point_colour = paper_palette.get(
        row[paper_col],
        'gray'
    )

    ax.scatter(
        row['x_bar'] + jitter,
        row['percentage'],

        color=point_colour,

        edgecolor='black',
        linewidth=0.7,

        s=65,
        alpha=0.9,

        zorder=10
    )


# =========================================================
# 17. Format X-axis
# =========================================================

ax.set_xticks(
    range(
        len(celltype_order)
    )
)


ax.set_xticklabels(
    celltype_order,
    rotation=45,
    ha='right'
)


# =========================================================
# 18. Format Y-axis
# =========================================================

ax.set_ylabel(
    'Percentage of cells per sample (%)',
    fontsize=12
)


ax.set_xlabel(
    'Fibroblast cell type',
    fontsize=12
)


ax.set_title(
    'Fibroblast cell-type composition across samples\n'
    'Mean ± standard deviation',
    fontsize=14
)


ax.set_ylim(
    bottom=0
)


# =========================================================
# 19. Remove unnecessary spines
# =========================================================

ax.spines['top'].set_visible(False)

ax.spines['right'].set_visible(False)


# =========================================================
# 20. Add horizontal grid
# =========================================================

ax.grid(
    axis='y',
    linestyle='--',
    alpha=0.3,
    zorder=0
)


# =========================================================
# 21. Create condition legend
# =========================================================

condition_legend = [
    Patch(
        facecolor=condition_palette[condition],
        edgecolor='black',
        alpha=0.35,
        label=condition
    )

    for condition in condition_order

    if condition in plot_df[
        condition_col
    ].unique()
]


# =========================================================
# 22. Create paper-source legend
# =========================================================

paper_legend = [
    Line2D(
        [0],
        [0],

        marker='o',

        linestyle='None',

        label=paper,

        markerfacecolor=colour,

        markeredgecolor='black',

        markersize=8
    )

    for paper, colour
    in paper_palette.items()

    if paper in plot_df[
        paper_col
    ].unique()
]


# =========================================================
# 23. Add condition legend
# =========================================================

legend1 = ax.legend(
    handles=condition_legend,

    title='Condition',

    bbox_to_anchor=(
        1.02,
        1
    ),

    loc='upper left',

    frameon=False
)


# Keep condition legend while adding second legend
ax.add_artist(legend1)


# =========================================================
# 24. Add paper-source legend
# =========================================================

ax.legend(
    handles=paper_legend,

    title='Paper source',

    bbox_to_anchor=(
        1.02,
        0.68
    ),

    loc='upper left',

    frameon=False
)


# =========================================================
# 25. Adjust layout
# =========================================================

plt.tight_layout()


# =========================================================
# 26. Display figure
# =========================================================

plt.show()


# In[39]:


# Convert values to strings
for col in required_columns:
    df[col] = df[col].astype(str).str.strip()

# Clean paper source names
df[paper_col] = (
    df[paper_col]
    .str.strip()
    .replace({
        'NDF': 'NDFs',
        'NDFs ': 'NDFs',
        'NDFs': 'NDFs',
        'Deng': 'Deng et al',
        'Deng et al.': 'Deng et al',
        'Li': 'Li et al',
        'Li et al.': 'Li et al',
        'Direder': 'Direder et al',
        'Direder et al.': 'Direder et al'
    })
)


# In[44]:





# In[128]:


adata_custom.obs['sample_status'] = (
    adata_custom.obs['orig.ident']
    .map(status_mapping)
)


# In[33]:


# Rename Leiden clusters into fibroblast types
adata_custom.obs["fbs_type"] = adata_custom.obs["leiden"].astype(str).replace({
    "0": "Mesenchymal",
    "1": "Proinflammatory",
    "2": "Myofibroblasts",
    "3": "Secretory Papillary",
    "4": "Proliferating"
})

# Set order for the plot
cluster_order = [
    "Mesenchymal",
    "Proinflammatory",
    "Myofibroblasts",
    "Secretory Papillary",
    "Proliferating"
]

adata_custom.obs["fbs_type"] = (
    adata_custom.obs["fbs_type"]
    .astype("category")
    .cat.set_categories(cluster_order, ordered=True)
)

# Check
print(adata_custom.obs["fbs_type"].value_counts())


# In[60]:


import scanpy as sc

markers = {
    "Mesenchymal": [
        "POSTN", "ASPN", "ADAM12", "COL1A1", "OGN"
    ],

    "Myofibroblasts": [
        "ACTA2", "RGS5", "TAGLN", "MYL9"
    ],

    "Proinflammatory": [
        "APOD", "CXCL14", "C3", "TXNIP", "CCL19"
    ],

    "Secretory Papillary": [
        "ZFP36", "ZFP36L1", "CFD", "DCN"
    ],

    "Proliferating": [
        "NME2", "CRIP1", "EEF1G", "SNHG29", "SEPTIN7"
    ]
}

cluster_order = [
    "Mesenchymal",
    "Myofibroblasts",
    "Proinflammatory",
    "Secretory Papillary",
    "Proliferating"
]

adata_fibro.obs["fbs_type"] = adata_fibro.obs["leiden"].astype(str).replace({
    "0": "Mesenchymal",
    "1": "Proinflammatory",
    "2": "Myofibroblasts",
    "3": "Secretory Papillary",
    "4": "Proliferating"
})

adata_fibro.obs["fbs_type"] = (
    adata_fibro.obs["fbs_type"]
    .astype("category")
    .cat.set_categories(cluster_order, ordered=True)
)

sc.pl.dotplot(
    adata_fibro,
    var_names=markers,
    groupby="fbs_type",
    standard_scale="var",
    swap_axes=False,
    dendrogram=False,
    dot_max=0.8,
    dot_min=0.05,
    color_map="Reds",
    save="_fibroblast_subcluster_marker_dotplot.png"
)


# In[68]:


mesenchymal_genes = [
    'POSTN', 'COL1A1', 'COL3A1',
    'ASPN', 'DCN', 'LUM',
    'OGN', 'ADAM12'
]

adipogenic_genes = [
    'APOD', 'CFD', 'FABP4',
    'PPARG', 'CEBPA',
    'PLIN2', 'ADIPOQ'
]

chondrogenic_genes = [
    'COMP', 'SOX9', 'RUNX2',
    'BMP2', 'COL2A1',
    'ACAN', 'MATN3'
]


# In[97]:


import scanpy as sc

sc.tl.score_genes(
    adata_keloid_only,
    mesenchymal_genes,
    score_name='Mesenchymal_signature'
)

sc.tl.score_genes(
    adata_keloid_only,
    adipogenic_genes,
    score_name='Adipogenic_signature'
)

sc.tl.score_genes(
    adata_keloid_only,
    chondrogenic_genes,
    score_name='Chondrogenic_signature'
)


# In[98]:


sc.pl.umap(
    adata_keloid_only,
    color=[
        'Mesenchymal_signature',
        'Adipogenic_signature',
        'Chondrogenic_signature'
    ],
    cmap='RdYlBu_r',
    vmin='p1',
    vmax='p99'
)


# In[26]:


import pandas as pd
import numpy as np
import scipy.sparse as sp

adata = adata_custom.copy()

sample_col = "orig.ident"
condition_col = "Group_merged"
cluster_col = "leiden"

adata.obs[sample_col] = adata.obs[sample_col].astype(str)
adata.obs[condition_col] = adata.obs[condition_col].astype(str)
adata.obs[cluster_col] = adata.obs[cluster_col].astype(str)

pseudo_expr = []
metadata = []

for cluster in sorted(adata.obs[cluster_col].unique()):

    adata_cluster = adata[adata.obs[cluster_col] == cluster].copy()

    for sample in adata_cluster.obs[sample_col].unique():

        sample_cells = adata_cluster[
            adata_cluster.obs[sample_col] == sample
        ]

        X = sample_cells.X

        if sp.issparse(X):
            mean_expr = np.asarray(X.mean(axis=0)).flatten()
        else:
            mean_expr = np.asarray(X.mean(axis=0)).flatten()

        sample_name = f"{sample}_cluster{cluster}"

        pseudo_expr.append(
            pd.Series(
                mean_expr,
                index=adata.var_names,
                name=sample_name
            )
        )

        metadata.append({
            "sample": sample_name,
            "orig_sample": sample,
            "cluster": cluster,
            "condition": sample_cells.obs[condition_col].iloc[0],
            "n_cells": sample_cells.n_obs
        })

pseudo_expr_df = pd.concat(pseudo_expr, axis=1)
meta_df = pd.DataFrame(metadata).set_index("sample")

pseudo_expr_df.to_csv("scaled_pseudobulk_mean_expression.csv")
meta_df.to_csv("scaled_pseudobulk_metadata.csv")

print(pseudo_expr_df.shape)
print(meta_df.groupby(["cluster", "condition"]).size())


# In[24]:


print(adata_custom.layers.keys())
print(adata_custom.X.min(), adata_custom.X.max())


# In[19]:


adata_custom.obs["Group"].value_counts()


# In[20]:


adata_custom.obs["Group_merged"] = (
    adata_custom.obs["Group"]
    .replace({
        "CTRL": "Normal",
        "CASE": "Keloid"
    })
)


# In[21]:


adata_custom.obs["Group_merged"].value_counts()


# In[59]:


hox_signature = ["GLI2","COMP","RUNX2","MMP16","COL5A2","FBN1","EPHA2","COL6A1"]
# Keep only genes present in your dataset
hox_signature = [g for g in hox_signature if g in adata_keloid_only.var_names]

import scanpy as sc

sc.tl.score_genes(
    adata_keloid_only,
    gene_list=hox_signature,
    score_name="HOX_signature"
)




# In[60]:


sc.pl.umap(
    adata_keloid_only,
    color="HOX_signature",
    cmap="RdYlBu_r"
)


# In[83]:


genes = [
    "CEBPA", "ZNF423", "PPARG",
    "CEBPB", "EBF2","CDH11"]
genes = [g for g in genes if g in data_normal_only.var_names]

import scanpy as sc

sc.tl.score_genes(
    data_normal_only,
    gene_list=genes,
    score_name="genes"
)


# In[70]:


sc.pl.umap(
    adata_keloid_only,
    color="genes",
    cmap="RdYlBu_r"
)


# In[78]:


sc.pl.umap(
    adata_keloid_only,
    color="genes",
    cmap="RdYlBu_r",
    vmin=-1,
    vmax=2
)


# In[84]:


sc.pl.umap(
    data_normal_only,
    color="genes",
    cmap="RdYlBu_r",
    vmin=-1,
    vmax=2
)


# In[32]:


location_map = {
    # Liu
    "K007CASE": "chest",
    "K007CTRL": "chest",
    "K009CASE": "chest",
    "K009CTRL": "chest",
    "K012CASE": "chest",
    "K012CTRL": "chest",
    "K013CASE": "chest",
    "K013CTRL": "chest",

    # Deng
    "KF1": "chest",
    "KF2": "chest",
    "KF3": "back",
    "NF1": "chest",
    "NF2": "chest",
    "NF3": "back",

    # Direder
    "KD1": "chest",
    "KD2": "earlobe",
    "KD3": "earlobe",
    "KD4": "earlobe",
    "ND1": "abdomen"
}


# In[33]:


adata_fibro_no45.obs["anatomical_location"] = (
    adata_fibro_no45.obs["orig.ident"].map(location_map)
)

adata_fibro_no45.obs["anatomical_location"] = (
    adata_fibro_no45.obs["anatomical_location"].astype("category")
)


# In[34]:


adata_fibro_no45.obs[["orig.ident", "anatomical_location"]].drop_duplicates()


# In[35]:


adata_fibro_no45.obs["anatomical_location"].value_counts()


# In[36]:


location_map = {
    "K007CASE": "chest",
    "K007CTRL": "chest",
    "K009CASE": "chest",
    "K009CTRL": "chest",
    "K012CASE": "chest",
    "K012CTRL": "chest",
    "K013CASE": "chest",
    "K013CTRL": "chest",

    "KF1": "chest",
    "KF2": "chest",
    "KF3": "back",
    "NF1": "chest",
    "NF2": "chest",
    "NF3": "back",

    "Kd1": "chest",
    "Kd2": "earlobe",
    "Kd3": "earlobe",
    "Kd4": "earlobe",

    "ND1": "abdomen",
    "Nsc1": "abdomen",
    "Nsc2": "abdomen",
    "Nsc3": "abdomen",
    "Nsk1": "abdomen"
}

adata_fibro_no45.obs["anatomical_location"] = (
    adata_fibro_no45.obs["orig.ident"]
    .astype(str)
    .map(location_map)
)

adata_fibro_no45.obs["anatomical_location"].value_counts(dropna=False)


# In[37]:


adata_fibro_no45.obs["orig.ident"].unique()


# In[38]:


sc.pl.umap(
    adata_fibro_no45,
    color="anatomical_location",
    legend_loc="right margin"
)


# In[39]:


hox_genes = [
    "HOXA1", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXA10",
    "HOXB1", "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB9",
    "HOXC4", "HOXC5", "HOXC6", "HOXC8", "HOXC9", "HOXC10",
    "HOXD3", "HOXD4", "HOXD8", "HOXD9", "HOXD10", "HOXD11"
]

# keep only genes present in your dataset
hox_genes = [g for g in hox_genes if g in adata_fibro_no45.var_names]

print(hox_genes)


# In[40]:


import pandas as pd
import scanpy as sc

hox_avg = pd.DataFrame(
    adata_fibro_no45[:, hox_genes].X.toarray(),
    columns=hox_genes,
    index=adata_fibro_no45.obs["anatomical_location"]
)

hox_avg["location"] = hox_avg.index

hox_avg = hox_avg.groupby("location").mean()

hox_avg


# In[61]:


hox_genes = [
    "HOXA1", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXA10",
    "HOXB1", "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB9",
    "HOXC4", "HOXC5", "HOXC6", "HOXC8", "HOXC9", "HOXC10",
    "HOXD3", "HOXD4", "HOXD8", "HOXD9", "HOXD10", "HOXD11"
]

# keep genes that exist in your dataset
hox_genes = [g for g in hox_genes if g in adata_fibro_no45.var_names]

sc.pl.dotplot(
    adata_fibro_no45,
    var_names=hox_genes,
    groupby="anatomical_location",
    standard_scale="var",
    cmap="RdYlBu_r",
    dot_max=0.8,
    figsize=(12,6)
)


# In[62]:


adata_fibro_no45.obs["location_status"].value_counts()


# In[64]:


location_map = {
    "K007CASE": "chest",
    "K007CTRL": "chest",
    "K009CASE": "chest",
    "K009CTRL": "chest",
    "K012CASE": "chest",
    "K012CTRL": "chest",
    "K013CASE": "chest",
    "K013CTRL": "chest",

    "KF1": "chest",
    "KF2": "chest",
    "KF3": "back",
    "NF1": "chest",
    "NF2": "chest",
    "NF3": "back",

    "Kd1": "chest",
    "Kd2": "earlobe",
    "Kd3": "earlobe",
    "Kd4": "earlobe",

    "ND1": "abdomen",
    "Nsc1": "abdomen",
    "Nsc2": "abdomen",
    "Nsc3": "abdomen",
    "Nsk1": "abdomen"
}

adata_fibro_no45.obs["anatomical_location"] = (
    adata_fibro_no45.obs["orig.ident"]
    .astype(str)
    .map(location_map)
)

adata_fibro_no45.obs["anatomical_location"].value_counts(dropna=False)


# In[65]:


adata_fibro_no45.obs["location_status"] = (
    adata_fibro_no45.obs["sample_status"].astype(str)
    + "_"
    + adata_fibro_no45.obs["anatomical_location"].astype(str)
)

adata_fibro_no45.obs["location_status"] = (
    adata_fibro_no45.obs["location_status"].astype("category")
)


# In[67]:


import numpy as np
keloid_samples = [
    "K007CASE", "K009CASE", "K012CASE", "K013CASE",
    "KF1", "KF2", "KF3",
    "Kd1", "Kd2", "Kd3", "Kd4"
]

adata_fibro_no45.obs["sample_status"] = np.where(
    adata_fibro_no45.obs["orig.ident"].isin(keloid_samples),
    "Keloid",
    "Normal"
)


# In[68]:


keloid_samples = [
    "K007CASE", "K009CASE", "K012CASE", "K013CASE",
    "KF1", "KF2", "KF3",
    "Kd1", "Kd2", "Kd3", "Kd4"
]

adata_fibro_no45.obs["location_status"] = np.where(
    adata_fibro_no45.obs["orig.ident"].isin(keloid_samples),
    "Keloid_",
    "Normal_"
) + adata_fibro_no45.obs["anatomical_location"].astype(str)


# In[135]:


adata_fibro_no45.obs["location_status"].value_counts()


# In[69]:


sc.pl.dotplot(
    adata_fibro_no45,
    var_names=hox_genes,
    groupby="location_status",
    standard_scale="var",
    cmap="RdYlBu_r"
)


# In[45]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# =========================================================
# 1. Define metadata columns
# =========================================================

sample_col = 'orig.ident'
cluster_col = 'leiden'
paper_col = 'paper_source'
condition_col = 'sample_status'


# =========================================================
# 2. Map Leiden clusters to fibroblast cell types
# =========================================================

celltype_map = {
    '0': 'Mesenchymal',
    '1': 'Proinflammatory',
    '2': 'Myofibroblast',
    '3': 'Secretory Papillary',
}


# =========================================================
# 3. Check that all required columns exist
# =========================================================

required_columns = [
    sample_col,
    cluster_col,
    paper_col,
    condition_col
]

missing_columns = [
    col
    for col in required_columns
    if col not in adata_fibro_no45.obs.columns
]

if missing_columns:
    raise KeyError(
        f"These columns are missing from adata_fibro_no45.obs: "
        f"{missing_columns}\n\n"
        f"Available columns:\n"
        f"{adata_fibro_no45.obs.columns.tolist()}"
    )


# =========================================================
# 4. Create working dataframe
# =========================================================

df = adata_fibro_no45.obs[required_columns].copy()

for col in required_columns:
    df[col] = df[col].astype(str)

df = df[
    ~df[required_columns].isin(
        ['nan', 'None', '']
    ).any(axis=1)
].copy()


# =========================================================
# 5. Map Leiden clusters to fibroblast cell types
# =========================================================

df['cell_type'] = df[cluster_col].map(celltype_map)

df = df.dropna(
    subset=['cell_type']
).copy()


# =========================================================
# 6. Inspect samples and paper sources
# =========================================================

print("\n========================================")
print("Samples, conditions and papers")
print("========================================")

print(
    df[
        [
            sample_col,
            condition_col,
            paper_col
        ]
    ]
    .drop_duplicates()
    .sort_values(
        [
            condition_col,
            paper_col,
            sample_col
        ]
    )
    .to_string(index=False)
)


print("\n========================================")
print("Unique paper sources")
print("========================================")

print(
    df[
        [
            paper_col,
            condition_col
        ]
    ]
    .drop_duplicates()
    .sort_values(
        [
            paper_col,
            condition_col
        ]
    )
    .to_string(index=False)
)


# =========================================================
# 7. Count cells per sample, condition and cell type
# =========================================================

counts = (
    df.groupby(
        [
            sample_col,
            condition_col,
            'cell_type'
        ],
        observed=True
    )
    .size()
    .reset_index(name='n_cells')
)


# =========================================================
# 8. Calculate total cells per sample
# =========================================================

totals = (
    df.groupby(
        [
            sample_col,
            condition_col
        ],
        observed=True
    )
    .size()
    .reset_index(name='total_cells')
)


# =========================================================
# 9. Calculate percentage of each fibroblast type per sample
# =========================================================

plot_df = counts.merge(
    totals,
    on=[
        sample_col,
        condition_col
    ],
    how='left'
)


plot_df['percentage'] = (
    plot_df['n_cells']
    / plot_df['total_cells']
    * 100
)


# =========================================================
# 10. Add zero values for missing fibroblast types
# =========================================================

celltype_order = [
    'Mesenchymal',
    'Myofibroblast',
    'Proinflammatory',
    'Secretory Papillary'
]


sample_metadata = (
    df[
        [
            sample_col,
            condition_col,
            paper_col
        ]
    ]
    .drop_duplicates(
        subset=[sample_col]
    )
)


complete_index = pd.MultiIndex.from_product(
    [
        sample_metadata[sample_col].unique(),
        celltype_order
    ],
    names=[
        sample_col,
        'cell_type'
    ]
)


complete_df = (
    complete_index
    .to_frame(index=False)
    .merge(
        sample_metadata,
        on=sample_col,
        how='left'
    )
)


plot_df = complete_df.merge(
    plot_df[
        [
            sample_col,
            condition_col,
            'cell_type',
            'n_cells',
            'total_cells',
            'percentage'
        ]
    ],
    on=[
        sample_col,
        condition_col,
        'cell_type'
    ],
    how='left'
)


plot_df['n_cells'] = (
    plot_df['n_cells']
    .fillna(0)
)


plot_df['percentage'] = (
    plot_df['percentage']
    .fillna(0)
)


# =========================================================
# 11. Restore total cell numbers
# =========================================================

sample_totals = (
    totals[
        [
            sample_col,
            'total_cells'
        ]
    ]
    .drop_duplicates()
)


plot_df = (
    plot_df
    .drop(
        columns=['total_cells'],
        errors='ignore'
    )
    .merge(
        sample_totals,
        on=sample_col,
        how='left'
    )
)


# =========================================================
# 12. Keep only Keloid and Normal samples
# =========================================================

plot_df = plot_df[
    plot_df[condition_col].isin(
        [
            'Keloid',
            'Normal'
        ]
    )
].copy()


# =========================================================
# 13. Define bar colours
# =========================================================

condition_palette = {
    'Keloid': '#1f77b4',
    'Normal': '#ff7f0e'
}


condition_order = [
    'Keloid',
    'Normal'
]


# =========================================================
# 14. Define PAPER colours
# =========================================================
#
# IMPORTANT:
# Points are coloured according to paper_source,
# NOT according to individual sample.
#
# Add every paper name that appears in your dataset here.
# =========================================================

paper_palette = {
    'Deng et al': '#1f77b4',       # blue
    'Li et al': '#d62728',         # red
    'Direder et al': '#2ca02c',    # green
    'Liu et al': '#9467bd',        # purple
    'NDFs': '#ff7f0e'              # orange
}


# =========================================================
# 15. Automatically warn if a paper is missing from palette
# =========================================================

papers_in_data = (
    plot_df[paper_col]
    .dropna()
    .unique()
    .tolist()
)


missing_papers = [
    paper
    for paper in papers_in_data
    if paper not in paper_palette
]


if missing_papers:

    print("\n========================================")
    print("WARNING: papers missing from paper_palette")
    print("========================================")

    for paper in missing_papers:
        print(paper)

    print(
        "\nThese papers will appear grey until you add "
        "them to paper_palette."
    )


# =========================================================
# 16. Define X-axis positions
# =========================================================

celltype_to_x = {
    cell_type: index
    for index, cell_type
    in enumerate(celltype_order)
}


# =========================================================
# 17. Define offsets for Keloid and Normal
# =========================================================

condition_offsets = {
    'Keloid': -0.2,
    'Normal': 0.2
}


plot_df['x_pos'] = (
    plot_df['cell_type']
    .map(celltype_to_x)
    .astype(float)
)


plot_df['offset'] = (
    plot_df[condition_col]
    .map(condition_offsets)
    .astype(float)
)


plot_df['x_bar'] = (
    plot_df['x_pos']
    + plot_df['offset']
)


# =========================================================
# 18. Calculate mean and SD across samples
# =========================================================

bar_df = (
    plot_df.groupby(
        [
            'cell_type',
            condition_col
        ],
        observed=True,
        as_index=False
    )
    .agg(
        mean_percentage=(
            'percentage',
            'mean'
        ),
        sd_percentage=(
            'percentage',
            'std'
        ),
        n_samples=(
            sample_col,
            'nunique'
        )
    )
)


bar_df['sd_percentage'] = (
    bar_df['sd_percentage']
    .fillna(0)
)


bar_df['x_pos'] = (
    bar_df['cell_type']
    .map(celltype_to_x)
    .astype(float)
)


bar_df['offset'] = (
    bar_df[condition_col]
    .map(condition_offsets)
    .astype(float)
)


bar_df['x_bar'] = (
    bar_df['x_pos']
    + bar_df['offset']
)


# =========================================================
# 19. Print summary statistics
# =========================================================

print("\n========================================")
print("Summary statistics")
print("========================================")

print(
    bar_df[
        [
            'cell_type',
            condition_col,
            'n_samples',
            'mean_percentage',
            'sd_percentage'
        ]
    ]
    .sort_values(
        [
            'cell_type',
            condition_col
        ]
    )
    .to_string(index=False)
)


# =========================================================
# 20. Create figure
# =========================================================

rng = np.random.default_rng(42)

fig, ax = plt.subplots(
    figsize=(14, 7)
)


# =========================================================
# 21. Plot mean bars with SD error bars
# =========================================================

for _, row in bar_df.iterrows():

    ax.bar(
        x=row['x_bar'],
        height=row['mean_percentage'],
        width=0.35,

        color=condition_palette.get(
            row[condition_col],
            'gray'
        ),

        alpha=0.35,

        edgecolor='black',
        linewidth=0.8,

        yerr=row['sd_percentage'],
        capsize=5,

        error_kw={
            'elinewidth': 1.5,
            'capthick': 1.5,
            'ecolor': 'black'
        },

        zorder=2
    )


# =========================================================
# 22. Add individual sample points
#     COLOUR = PAPER SOURCE
# =========================================================

for _, row in plot_df.iterrows():

    jitter = rng.uniform(
        -0.045,
        0.045
    )


    point_colour = paper_palette.get(
        row[paper_col],
        'gray'
    )


    ax.scatter(
        row['x_bar'] + jitter,
        row['percentage'],

        color=point_colour,

        edgecolor='black',
        linewidth=0.7,

        s=65,
        alpha=0.9,

        zorder=10
    )


# =========================================================
# 23. Format X-axis
# =========================================================

ax.set_xticks(
    range(
        len(celltype_order)
    )
)


ax.set_xticklabels(
    celltype_order,
    rotation=45,
    ha='right'
)


# =========================================================
# 24. Format Y-axis
# =========================================================

ax.set_ylabel(
    'Percentage of cells per sample (%)',
    fontsize=12
)


ax.set_xlabel(
    'Fibroblast cell type',
    fontsize=12
)


ax.set_title(
    'Fibroblast cell-type composition across samples\n'
    'Mean ± standard deviation',
    fontsize=14
)


ax.set_ylim(
    bottom=0
)


# =========================================================
# 25. Remove unnecessary spines
# =========================================================

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


# =========================================================
# 26. Add horizontal grid
# =========================================================

ax.grid(
    axis='y',
    linestyle='--',
    alpha=0.3,
    zorder=0
)


# =========================================================
# 27. Create condition legend
# =========================================================

condition_legend = [
    Patch(
        facecolor=condition_palette[condition],
        edgecolor='black',
        alpha=0.35,
        label=condition
    )

    for condition in condition_order

    if condition in plot_df[
        condition_col
    ].unique()
]


# =========================================================
# 28. Create PAPER SOURCE legend
# =========================================================

paper_legend = [
    Line2D(
        [0],
        [0],

        marker='o',

        linestyle='None',

        label=paper,

        markerfacecolor=colour,

        markeredgecolor='black',

        markersize=8
    )

    for paper, colour
    in paper_palette.items()

    if paper in plot_df[
        paper_col
    ].unique()
]


# =========================================================
# 29. Add condition legend
# =========================================================

legend1 = ax.legend(
    handles=condition_legend,

    title='Condition',

    bbox_to_anchor=(
        1.02,
        1
    ),

    loc='upper left',

    frameon=False
)


ax.add_artist(legend1)


# =========================================================
# 30. Add paper-source legend
# =========================================================

ax.legend(
    handles=paper_legend,

    title='Paper source',

    bbox_to_anchor=(
        1.02,
        0.68
    ),

    loc='upper left',

    frameon=False
)


# =========================================================
# 31. Adjust layout
# =========================================================

plt.tight_layout()


# =========================================================
# 32. Display figure
# =========================================================

plt.show()


# In[70]:


adata_keloid = adata_fibro_no45[
    adata_fibro_no45.obs["location_status"].str.startswith("Keloid_")
].copy()

sc.pl.dotplot(
    adata_keloid,
    var_names=hox_genes,
    groupby="location_status",
    standard_scale="var",
    cmap="RdYlBu_r")


# In[75]:


adata_normal = adata_fibro_no45[
    adata_fibro_no45.obs["location_status"].str.startswith("Normal_")
].copy()

sc.pl.dotplot(
    adata_normal,
    var_names=hox_genes,
    groupby="location_status",
    standard_scale="var",
    cmap="RdYlBu_r")


# In[76]:


adata_keloid = adata_fibro_no45[
    adata_fibro_no45.obs["location_status"].str.startswith("Keloid_")
].copy()

sc.pl.dotplot(
    adata_keloid,
    var_names=hox_genes,
    groupby="location_status",
    standard_scale="var",
    cmap="RdYlBu_r",
    dot_min=0,
    dot_max=0.4,
    smallest_dot=0
)


# In[77]:


# Normal samples
adata_normal = adata_fibro_no45[
    adata_fibro_no45.obs["location_status"].str.startswith("Normal_")
].copy()

sc.pl.dotplot(
    adata_normal,
    var_names=hox_genes,
    groupby="location_status",
    standard_scale="var",
    cmap="RdYlBu_r",
    dot_min=0,
    dot_max=0.4,
    smallest_dot=0
)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[32]:


import pandas as pd
import scanpy as sc

# Read WITHOUT header - gets the actual gene list
keloid_genes = pd.read_csv("keloid_abdomen_DEG.csv", header=None)[0].dropna().tolist()
ndf_genes = pd.read_csv("NDF_abdomen_DEGlist.csv", header=None)[0].dropna().tolist()

print(f"Keloid genes: {len(keloid_genes)}")
print(f"NDF genes: {len(ndf_genes)}")
print("Keloid first 5:", keloid_genes[:5])
print("NDF first 5:", ndf_genes[:5])


# In[ ]:





# In[ ]:





# In[ ]:





# In[33]:


# Filter genes that exist in your scRNA-seq data
keloid_genes_in_data = [g for g in keloid_genes if g in adata.var_names]
ndf_genes_in_data = [g for g in ndf_genes if g in adata.var_names]

print(f"Keloid genes found: {len(keloid_genes_in_data)}/{len(keloid_genes)}")
print(f"NDF genes found: {len(ndf_genes_in_data)}/{len(ndf_genes)}")

# Test keloid signature across clusters
sc.tl.rank_genes_groups(adata_custom, groupby='leiden', genes=keloid_genes_in_data, method='wilcoxon')

# KEY VISUALIZATION - shows which clusters express your signature!
sc.pl.dotplot(adata_custom, keloid_genes_in_data[:25], groupby='leiden', 
              standard_scale='var', save='_keloid_signature.png')


# In[58]:


# Score cells based on your gene lists
sc.tl.score_genes(adata_custom, keloid_genes_in_data, score_name='Keloid_Score')
sc.tl.score_genes(adata_custom, ndf_genes_in_data, score_name='NDF_Score')

# Plot on UMAP - red/orange regions = high signature expression!
sc.pl.umap(adata_custom, color=['Keloid_Score', 'NDF_Score'], 
           cmap='viridis', size=50)


# In[35]:


# Score cells based on your gene lists
sc.tl.score_genes(adata_custom, keloid_genes_in_data, score_name='Keloid_Score')
sc.tl.score_genes(adata_custom, ndf_genes_in_data, score_name='NDF_Score')

# Plot on UMAP - red/orange regions = high signature expression!
sc.pl.umap(adata_custom, color=['Keloid_Score', 'NDF_Score'], 
           cmap='viridis', size=50)


# In[36]:


# Violin plot for signature SCORES across all clusters
sc.pl.violin(adata_custom, ['Keloid_Score', 'NDF_Score'], 
             groupby='leiden', 
             rotation=90,  # Rotate cluster labels
             stripplot=False)  # Cleaner without dots


# In[ ]:


# Define your fibroblast names (must match leiden categories exactly)
new_names = ['Mesenchymal', 'Proinflammatory', 'SMC/Myofibroblasts', 
             'Secretory Papillary', 'Proliferating']

# Apply the names (assumes leiden has categories 0,1,2,3,4)
adata_fibro_no2.obs['leiden'] = adata_fibro_no2.obs['leiden'].astype('category')
adata_fibro_no2.rename_categories('leiden', new_names)

# Verify it worked
print(adata_fibro_no2.obs['leiden'].cat.categories)


# In[37]:


# Calculate % cells expressing your gene list per paper
def genes_expressed_prop(adata_subset, genes):
    expressed_mask = adata_subset[:, genes].X > 0
    return expressed_mask.sum(1) / len(genes)  # Proportion of genes ON per cell

# Per paper proportions
papers = ['Deng et al', 'Li et al', 'Direder et al']
proportions = {}

for paper in papers:
    paper_mask = adata_custom.obs['paper_source'] == paper
    if paper_mask.sum() > 0:
        prop = genes_expressed_prop(adata_custom[paper_mask], keloid_genes_in_data).mean()
        proportions[paper] = prop * 100
        print(f"{paper}: {prop*100:.1f}% cells express avg {len(keloid_genes_in_data)} genes")

# Plot
plt.figure(figsize=(8, 6))
pd.Series(proportions).plot(kind='bar')
plt.title('Keloid Signature Genes Expressed - % Cells per Paper')
plt.ylabel('% Cells Expressing Genes')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('gene_list_proportion_per_paper.png', dpi=300)
plt.show()


# In[101]:


# % MEAN expression level (more realistic)
mean_expression = adata_custom[:, keloid_genes_in_data].X.mean(axis=1)
adata_custom.obs['keloid_mean_expr'] = mean_expression

# Per paper_cluster MEAN expression
results_mean = adata_custom.obs.groupby('paper_cluster')['keloid_mean_expr'].mean() * 100
print("MEAN keloid gene expression per cluster per paper:")
print(results_mean.round(2))

# Heatmap of mean expression
results_df = results_mean.reset_index()
results_df['Paper'] = results_df['paper_cluster'].str.split('_').str[0]
results_df['Cluster'] = results_df['paper_cluster'].str.split('_').str[1]

pivot_mean = results_df.pivot(index='Paper', columns='Cluster', values='keloid_mean_expr')
plt.figure(figsize=(10, 6))
sns.heatmap(pivot_mean, annot=True, fmt='.2f', cmap='Reds', 
            cbar_kws={'label': 'Mean Expression'})
plt.title('Mean Keloid Signature Expression in Fibroblasts - Per Paper Per Cluster')
plt.savefig('keloid_mean_expression_per_paper_cluster.png', dpi=300)
plt.show()


# In[94]:


# % MEAN expression level (more realistic)
mean_expression = adata_keloid_only[:, keloid_genes_in_data].X.mean(axis=1)
adata_keloid_only.obs['keloid_mean_expr'] = mean_expression

# Per paper_cluster MEAN expression
results_mean = adata_keloid_only.obs.groupby('paper_cluster')['keloid_mean_expr'].mean() * 100
print("MEAN keloid gene expression per cluster per paper:")
print(results_mean.round(2))

# Heatmap of mean expression
results_df = results_mean.reset_index()
results_df['Paper'] = results_df['paper_cluster'].str.split('_').str[0]
results_df['Cluster'] = results_df['paper_cluster'].str.split('_').str[1]

pivot_mean = results_df.pivot(index='Paper', columns='Cluster', values='keloid_mean_expr')
plt.figure(figsize=(10, 6))
sns.heatmap(pivot_mean, annot=True, fmt='.2f', cmap='Reds', 
            cbar_kws={'label': 'Mean Expression'})
plt.title('Mean Keloid Signature Expression in Keloid cells - Per Paper Per Cluster')
plt.savefig('keloid_mean_expression_per_paper_cluster.png', dpi=300)
plt.show()


# In[72]:


# 1. Find top 5 markers PER cluster (paper_cluster level)
sc.tl.rank_genes_groups(adata_keloid_only, groupby='paper_cluster', method='wilcoxon')

# 2. Get marker table
markers_df = sc.get.rank_genes_groups_df(adata_keloid_only, group=None)
print("Top 5 markers per paper_cluster:")
print(markers_df.groupby('group').head(15)[['group', 'names', 'logfoldchanges']].head(100))


# In[ ]:


orig.ident
K009CTRL    877
K012CTRL    835
K007CTRL    688
K013CTRL    349
K009CASE      0
K007CASE      0
K012CASE      0
K013CASE      0
KF1           0
KF2           0
KF3           0
Kd1           0
Kd2           0
Kd3           0
Kd4           0
NF1           0
NF2           0
NF3           0
Nsc1          0
Nsc2          0
Nsc3          0
Nsk1          0


# In[86]:


# Keloid-only subset (from your custom 41K cells)
Normal_only_mask = adata_custom.obs['sample_status'] == 'Normal'
adata_Normal_only = adata_custom[Normal_only_mask].copy()

print("New keloid-only dataset:")
print(f"Total cells: {adata_Normal_only.n_obs:,}")
print("Samples:")
print(adata_Normal_only.obs['orig.ident'].value_counts())
print("Clusters preserved:")
print(adata_Normal_only.obs['leiden'].value_counts())


# In[92]:


# % MEAN expression level (more realistic)
mean_expression = adata_keloid_only[:, ndf_genes_in_data].X.mean(axis=1)
adata_keloid_only.obs['keloid_mean_expr'] = mean_expression

# Per paper_cluster MEAN expression
results_mean = adata_keloid_only.obs.groupby('paper_cluster')['keloid_mean_expr'].mean() * 100
print("MEAN keloid gene expression per cluster per paper:")
print(results_mean.round(2))

# Heatmap of mean expression
results_df = results_mean.reset_index()
results_df['Paper'] = results_df['paper_cluster'].str.split('_').str[0]
results_df['Cluster'] = results_df['paper_cluster'].str.split('_').str[1]

pivot_mean = results_df.pivot(index='Paper', columns='Cluster', values='keloid_mean_expr')
plt.figure(figsize=(10, 6))
sns.heatmap(pivot_mean, annot=True, fmt='.2f', cmap='Reds', 
            cbar_kws={'label': 'Mean Expression'})
plt.title('Mean Normal Signature Expression in Keloid cells- Per Paper Per Cluster')
plt.savefig('keloid_mean_expression_per_paper_cluster.png', dpi=300)
plt.show()


# In[97]:


# % MEAN expression level (more realistic)
mean_expression = adata_Normal_only[:, ndf_genes_in_data].X.mean(axis=1)
adata_Normal_only.obs['keloid_mean_expr'] = mean_expression

# Per paper_cluster MEAN expression
results_mean = adata_Normal_only.obs.groupby('paper_cluster')['keloid_mean_expr'].mean() * 100
print("MEAN keloid gene expression per cluster per paper:")
print(results_mean.round(2))

# Heatmap of mean expression
results_df = results_mean.reset_index()
results_df['Paper'] = results_df['paper_cluster'].str.split('_').str[0]
results_df['Cluster'] = results_df['paper_cluster'].str.split('_').str[1]

pivot_mean = results_df.pivot(index='Paper', columns='Cluster', values='keloid_mean_expr')
plt.figure(figsize=(10, 6))
sns.heatmap(pivot_mean, annot=True, fmt='.2f', cmap='Reds', 
            cbar_kws={'label': 'Mean Expression'})
plt.title('Mean Normal Signature Expression in NDF- Per Paper Per Cluster')
plt.savefig('keloid_mean_expression_per_paper_cluster.png', dpi=300)
plt.show()


# In[26]:


# Your ATAC-seq results might look like:
atac_results = pd.read_csv("atacseqsignature.csv")
signature_genes = atac_results.head(100)['Gene List'].tolist()
# ['FOXO1', 'COL1A1', 'TGFB1', 'ACTA2', ...]  # Your 100 most significant genes


# In[27]:


sig_genes = [g for g in signature_genes if g in adata_custom.var_names]


# In[28]:


print("Actual column names:")
print(atac_results.columns.tolist())
print("\nFirst few rows:")
print(atac_results.head())


# In[29]:


print(f"Found {len(sig_genes)}/{len(signature_genes)} genes in scRNA-seq data")


# In[30]:


import numpy as np
from scipy import sparse

cluster_hits = {}
for cluster in adata_custom.obs['leiden'].unique():
    # Integer indices, NOT boolean mask
    cluster_cells = np.where(adata_custom.obs['leiden'] == cluster)[0]

    # Use adata.X directly (skip .raw entirely)
    cluster_data = adata_custom[cluster_cells, sig_genes].X

    if sparse.issparse(cluster_data):
        cluster_expr = (cluster_data > 0).mean(axis=0).A1
    else:
        cluster_expr = (cluster_data > 0).mean(axis=0)

    num_hits = np.sum(cluster_expr > 0)
    cluster_hits[cluster] = num_hits

top_cluster = max(cluster_hits, key=cluster_hits.get)
print(f"🏆 Top cluster: {top_cluster} with {cluster_hits[top_cluster]}/{len(sig_genes)} genes")


# In[31]:


# Highlight top cluster
adata_custom.obs['atac_highlight'] = adata_custom.obs['leiden'] == top_cluster

# Signature score + clusters
sc.tl.score_genes(adata_custom, sig_genes, score_name='ATAC_signature')
sc.pl.umap(adata_custom, color=['ATAC_signature', 'leiden', 'atac_highlight'], 
           ncols=3, color_map='viridis', size=8, frameon=False,
           title=['ATAC Score', 'Clusters', f'Top: {top_cluster}'])


# In[32]:


sc.pl.dotplot(adata_custom, sig_genes[:20], groupby='leiden')


# In[33]:


sc.pl.violin(adata_custom, 'ATAC_signature', groupby='leiden')


# In[36]:


import pandas as pd

atac_df = pd.read_csv("/data/home/hmz255/ATACGENES.csv")

# adjust column name if needed
atac_genes = atac_df.iloc[:, 2].dropna().unique().tolist()

print(len(atac_genes))


# In[77]:


import pandas as pd

atac_df = pd.read_csv("/data/home/hmz255/ATACGENES.csv")

# KEEP ONLY SIGNIFICANT ONES (adjust column name!)
sig_atac = atac_df[atac_df['padj'] < 0.05]

atac_genes = (
    sig_atac.iloc[:, 2]
    .dropna()
    .astype(str)
    .str.upper()
    .str.strip()
    .unique()
    .tolist()
)

print("Signature size:", len(atac_genes))


# In[78]:


import scanpy as sc

sc.tl.score_genes(
    adata_custom,
    gene_list=atac_genes,
    score_name="ATAC_score"
)


# In[79]:


cluster_scores = adata_custom.obs.groupby('leiden')['ATAC_score'].mean().sort_values(ascending=False)

print(cluster_scores)


# In[80]:


sc.pl.umap(
    adata_custom,
    color=["ATAC_score", "leiden"],
    cmap="viridis",
    size=8
)


# In[81]:


import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(6,4))
sns.barplot(x=cluster_scores.values, y=cluster_scores.index, palette="viridis")
plt.title("ATAC signature per cluster")
plt.xlabel("Average ATAC score")
plt.ylabel("Cluster")
plt.show()


# In[82]:


import scanpy as sc

sc.pl.violin(
    adata_custom,
    keys="ATAC_score",
    groupby="leiden",
    stripplot=False,
    jitter=0.2,
    scale="width",
    inner="quartile",
    rotation=90
)


# In[83]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# compute mean ATAC score per cluster
plot_df = adata_custom.obs.groupby('leiden')['ATAC_score'].mean().reset_index()

# sort for nicer plot
plot_df = plot_df.sort_values('ATAC_score', ascending=False)

plt.figure(figsize=(7,4))
sns.barplot(
    data=plot_df,
    x='ATAC_score',
    y='leiden',
    palette='viridis'
)

plt.axvline(0, color='black', linestyle='--')
plt.title("ATAC signature per cluster")
plt.xlabel("Mean ATAC score")
plt.ylabel("Cluster")
plt.show()


# In[146]:


atac_genes = [g for g in atac_genes if g in adata_custom.var_names]


# In[148]:


sc.tl.score_genes(
    adata_custom,
    gene_list=atac_genes,
    score_name='ATAC_DA_score'
)


# In[150]:


sc.pl.umap(adata_custom, color='ATAC_DA_score')


# In[151]:


sc.pl.violin(adata_custom, keys='ATAC_DA_score', groupby='leiden')


# In[232]:


import gseapy as gp

enr = gp.enrichr(
    gene_list=atac_genes,
    gene_sets=[
        'KEGG_2021_Human'
    ],
    organism='h. sapiens'
)


# In[236]:


# adjust thresholds if needed
pval_cutoff = 0.05
min_genes = 5

filtered = enr.results[
    (enr.results['P-value'] < pval_cutoff) &
    (enr.results['Genes'].str.split(';').apply(len) >= min_genes)
]


# In[254]:


filtered.to_csv("atac_enriched_pathways.csv", index=False)


# In[243]:


filtered  # dict: {pathway_name: [genes]}


# In[241]:


for name, genes in gene_sets.items():
    if len(genes) >= 3:  # safety check
        sc.tl.score_genes(adata_custom, genes, score_name=name)


# In[188]:


sc.pl.violin(
    adata_custom,
    keys=list(gene_sets.keys()),
    groupby='leiden',
    rotation=90
)


# In[245]:


sc.pl.dotplot(
    adata_custom,
    var_names=list(filtered.keys()),
    groupby='leiden'
)


# In[231]:


sc.pl.dotplot(
    adata_keloid_only,
    var_names=list(gene_sets.keys()),
    groupby='leiden'
)


# In[247]:


filtered = enr.results[
    (enr.results['P-value'] < pval_cutoff) &
    (enr.results['Genes'].str.split(';').apply(len) >= min_genes)
]


# In[252]:


cluster_scores = adata_keloid_only.obs.groupby("leiden")[list(gene_sets.keys())].mean()


# In[253]:


go_terms = list(gene_sets.keys())

sc.pl.dotplot(
    adata_keloid_only,
    var_names=go_terms,
    groupby="leiden",
    standard_scale="var"
)


# In[289]:


import pandas as pd

tf_df = pd.read_csv("diffbinddetect.csv")

# adjust column name if needed
tf_genes = tf_df["motif"].astype(str).str.strip()


# In[290]:


print(filtered_tfs.columns)
print(filtered_tfs.head())


# In[291]:


import scanpy as sc
import pandas as pd

sig_tfs = (
    filtered_tfs.loc[filtered_tfs["KDF_NDF_highlighted"] == True, "motif"]
    .astype(str)
    .unique()
)

sig_tfs_in_scrna = [tf for tf in sig_tfs if tf in adata.var_names]

print("Significant TFs in CSV:", len(sig_tfs))
print("Significant TFs found in scRNA-seq:", len(sig_tfs_in_scrna))
print(sig_tfs_in_scrna)


# In[292]:


sc.pl.umap(
    adata_keloid_only,
    color=sig_tfs_in_scrna,
    ncols=4,
    frameon=False,
    cmap="viridis"
)


# In[276]:


cluster_col = "leiden"   # change if your cluster column is different

sc.pl.dotplot(
    adata_keloid_only,
    var_names=sig_tfs_in_scrna,
    groupby=cluster_col,
    standard_scale="var"
)


# In[ ]:


sc.tl.score_genes(
    adata_keloid_only,
    gene_list=sig_tfs_in_scrna,
    score_name="TF_avg_score"
)

sc.pl.umap(
    adata_keloid_only,
    color="TF_avg_score",
    cmap="viridis",
    frameon=False
)


# In[ ]:





# In[ ]:





# In[277]:


sig_tfs = (
    filtered_tfs.loc[filtered_tfs["KDF_NDF_highlighted"] == True, "motif"]
    .astype(str)
    .unique()
)

sig_tfs_in_scrna = [tf for tf in sig_tfs if tf in adata.var_names]

print(len(sig_tfs_in_scrna))


# In[278]:


sc.tl.score_genes(
    adata_custom,
    gene_list=sig_tfs_in_scrna,
    score_name="TF_signature"
)


# In[279]:


sc.pl.umap(
    adata_custom,
    color=["leiden", "TF_signature"],
    cmap="viridis",
    frameon=False
)


# In[280]:


cluster_col = "leiden"

cluster_signature = (
    adata_custom.obs
    .groupby(cluster_col)["TF_signature"]
    .mean()
    .sort_values(ascending=False)
)

cluster_signature


# In[281]:


sc.tl.rank_genes_groups(
    adata_custom,
    groupby="leiden",
    method="wilcoxon"
)


# In[288]:


cluster_col = "leiden"   # change if needed

sc.tl.rank_genes_groups(
    adata_custom,
    groupby=cluster_col,
    groups=["2"],
    reference="rest",
    method="wilcoxon"
)

sc.pl.rank_genes_groups(
    adata_custom,
    n_genes=40,
    sharey=False
)


# In[287]:


cluster_col = "leiden"   # change if needed

sc.tl.rank_genes_groups(
    adata_keloid_only,
    groupby=cluster_col,
    groups=["2"],
    reference="rest",
    method="wilcoxon"
)

sc.pl.rank_genes_groups(
    adata_keloid_only,
    n_genes=40,
    sharey=False
)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[23]:


# 1. Extract fibroblast clusters 0,3,4,10
fibro_clusters = [0, 3, 4, 10]
fibro_mask = adata_filtered.obs['seurat_clusters'].isin(fibro_clusters)
adata_fibro = adata_filtered[fibro_mask].copy()

print(f"Fibroblast clusters: {adata_fibro.n_obs} cells")


# In[25]:


# 2. Use TOP 10% CompFibro score within fibroblasts (not fixed cutoff)
comp_scores = adata_fibro.obs['CompFibro']
cutoff = comp_scores.quantile(0.9)  # Top 10%
comp_pos_mask = adata_fibro.obs['CompFibro'] > cutoff

comp_pos_cells = adata_fibro[comp_pos_mask].copy()
print(f"✅ COMP+ fibroblasts: {comp_pos_cells.n_obs} cells (top {100*(1-0.9):.0f}%)")
print(f"Cutoff: {cutoff:.2f}")


# In[26]:


# ONE LINER - bulletproof
mean_expr = np.asarray(comp_pos_cells.X.mean(axis=0)).flatten()
top50_comp_genes = pd.Series(mean_expr, index=comp_pos_cells.var_names).nlargest(100).index.tolist()


# In[27]:


# Remove .tolist() - it's already a list!
print("✅ Top 50 COMP+ fibroblast genes:")
print(top50_comp_genes[:10])  # Just slice the list!

# Score + plot
sc.tl.score_genes(adata_fibro, gene_list=top50_comp_genes, score_name='CompFibro_hi')
sc.pl.dotplot(adata_fibro, ['CompFibro_hi'], groupby='seurat_clusters', swap_axes=True)
sc.pl.umap(adata_fibro, color=['CompFibro_hi', 'seurat_clusters'])


# In[ ]:


# 1. Save as CSV (recommended)
top50_df = pd.DataFrame({'CompFibro_top50_genes': top50_comp_genes})
top50_df.to_csv('comp_fibroblast_top50_genes.csv', index=False)
print("✅ Saved to comp_fibroblast_top50_genes.csv")


# In[55]:


# Subset to POSTN-hi cells (top 25% POSTN expression)
postn_expr = adata_fibro_no2[:, 'POSTN'].X.mean(axis=1).A1
postn_hi_threshold = np.quantile(postn_expr, 0.75)
postn_hi = adata_fibro_no2[postn_expr > postn_hi_threshold].copy()

print(f"POSTN-hi cells: {len(postn_hi)} ({len(postn_hi)/len(adata_fibro_no2)*100:.1f}%)")

# YOUR exact Col1hi method - top 50 signature genes
mean_expr = postn_hi.X.mean(axis=0).A1
postn_hi_sig = pd.Series(mean_expr, index=postn_hi.var_names).nlargest(100).index.tolist()

print("✅ POSTN-hi signature (top 10):", postn_hi_sig[:10])
pd.Series(postn_hi_sig).to_csv('POSTN_hi_fibroblast_signature.csv')


# In[49]:


import numpy as np

# FIXED - check type first, no .A1 needed
postn_expr = adata_fibro_no2[:, 'POSTN'].X.mean(axis=1)

# If numpy array (no .A1), use .squeeze()
if hasattr(postn_expr, 'A1'):
    postn_expr = postn_expr.A1
postn_expr = postn_expr.squeeze()  # Remove extra dimensions

postn_hi_threshold = np.quantile(postn_expr, 0.75)
postn_hi = adata_fibro_no2[postn_expr > postn_hi_threshold].copy()


# In[56]:


# Print the COMPLETE top 50 POSTN-hi signature
print("POSTN-hi Keloid Fibroblast Signature (Top 50):")
print(postn_hi_sig)

# Also save as numbered list for paper
for i, gene in enumerate(postn_hi_sig, 1):
    print(f"{i:2d}. {gene}")


# In[58]:


# FIXED - safe POSTN-high extraction
postn_expr = keloid_fibro[:, 'POSTN'].X.mean(axis=1)

# Handle both sparse & dense automatically
if scipy.sparse.issparse(postn_expr):
    postn_expr = postn_expr.A1
postn_expr = postn_expr.squeeze()  # Remove dimensions

postn_hi_threshold = np.quantile(postn_expr, 0.75)
postn_hi_keloid = keloid_fibro[postn_expr > postn_hi_threshold].copy()

print(f"POSTN-hi keloid fibroblasts: {len(postn_hi_keloid)} cells")


# In[52]:


# 1. SUBSET TO KELOID FIBROBLASTS ONLY first
keloid_fibro = adata_fibro_no2[adata_fibro_no2.obs['Group'] == 'Keloid'].copy()
print(f"Keloid fibroblasts: {len(keloid_fibro)} cells")

# 2. POSTN-hi cells from keloid fibroblasts (top 25%)
postn_expr = keloid_fibro[:, 'POSTN'].X.mean(axis=1).A1
postn_hi_threshold = np.quantile(postn_expr, 0.75)
postn_hi_keloid = keloid_fibro[postn_expr > postn_hi_threshold].copy()

print(f"POSTN-hi keloid fibroblasts: {len(postn_hi_keloid)} cells")

# 3. YOUR exact Col1hi method - top 50 signature
mean_expr = postn_hi_keloid.X.mean(axis=0).A1
postn_hi_sig_keloid = pd.Series(mean_expr, index=postn_hi_keloid.var_names).nlargest(50).index.tolist()

print("✅ POSTN-hi KELOID fibroblast signature (top 10):", postn_hi_sig_keloid[:10])
pd.Series(postn_hi_sig_keloid).to_csv('POSTN_hi_keloid_signature.csv')


# In[ ]:


# Print COMPLETE top 50 POSTN-hi signature (run after signature extraction)
print("POSTN-hi Keloid Fibroblast Signature (Top 50):")
for i, gene in enumerate(postn_hi_sig_keloid, 1):
    print(f"{i:2d}. {gene}")

# Save formatted for paper (Table S2)
top50_df = pd.DataFrame({
    'Rank': range(1, 51),
    'Gene': postn_hi_sig_keloid
})
top50_df.to_csv('POSTN_hi_keloid_top50.csv', index=False)
print(f"\n✅ Saved: POSTN_hi_keloid_top50.csv ({len(postn_hi_sig_keloid)} genes)")


# In[ ]:


# Check unique values in Group column
print("Available groups:")
print(adata_fibroblasts.obs['Group'].value_counts())

# Also check all .obs columns
print("\nAll .obs columns:")
print(adata_fibroblasts.obs.columns.tolist())


# In[ ]:


# Merge CASE into Keloid, CTRL into Normal
adata_fibroblasts.obs['Condition_merged'] = adata_fibroblasts.obs['Group'].replace({
    'CASE': 'Keloid',
    'CTRL': 'Normal'
})

print("Merged groups:")
print(adata_fibroblasts.obs['Condition_merged'].value_counts())


# In[60]:


# 1. POSTN-hi fibroblasts (high POSTN expressors)
postn_hi = adata_fibro_no2[adata_fibro_no2[:, 'POSTN'].X.toarray().flatten() > 0].copy()

# 2. Top 50 POSTN-hi signature (mean expression)
mean_expr_postn = postn_hi.X.mean(axis=0)
if hasattr(mean_expr_postn, 'A1'):
    mean_expr_postn = mean_expr_postn.A1

postn_top50 = pd.Series(mean_expr_postn, index=postn_hi.var_names).nlargest(50).index.tolist()

# 3. Score POSTN-hi signature across ALL fibroblasts
sc.tl.score_genes(adata_fibro_no2, gene_list=postn_top50, score_name='POSTNhi')

# 4. Visualize in fibroblast clusters only
sc.pl.dotplot(adata_fibro_no2, ['POSTNhi'], 
              groupby='seurat_clusters', 
              swap_axes=True, save='_postnhi_dot.png')

sc.pl.umap(adata_fibro_no2, 
           color=['POSTNhi', 'seurat_clusters'], 
           legend_loc='on data', save='_postnhi_umap.png')


# In[ ]:


print("✅ POSTN-hi Top 50 genes:")
print(postn_top50)

print("\nTop 10 POSTN-hi signature genes:")
for i, gene in enumerate(postn_top50[:100], 1):
    print(f"{i:2d}. {gene}")


# In[ ]:


# Use your keloid fibroblast subset (clusters 0,3,4,10 from keloid data)
mean_expr_keloid = adata_fibroblasts.X.mean(axis=0)
if hasattr(mean_expr_keloid, 'A1'):
    mean_expr_keloid = mean_expr_keloid.A1

# Top 50 MOST EXPRESSED in keloid fibroblasts (clean data)
keloid_comp_top50 = pd.Series(mean_expr_keloid, index=adata_fibroblasts.var_names).nlargest(50).index.tolist()

print("✅ Keloid fibroblast Comp Top 50 genes:")
print(keloid_comp_top50[:50])  # Top 10

# Save it
adata_fibroblasts.uns['Keloid_CompFibro_signature'] = keloid_comp_top50
pd.Series(keloid_comp_top50).to_csv('Keloid_CompFibro_top50.csv', index=False)


# In[ ]:


sc.tl.score_genes(adata_fibroblasts, gene_list=keloid_comp_top50, score_name='KeloidComp')
sc.pl.umap(adata_fibroblasts, color='KeloidComp', save='_keloidcomp_umap.png')


# In[ ]:


# You need normal fibroblasts for comparison
# Assuming you have normal skin data or normal fibroblast clusters

# 1. Define normal fibroblasts (if you have them)
normal_fibro_clusters = ['1', '2']  # UPDATE with your normal fibroblast clusters
normal_fibros = adata_filtered[
    adata_filtered.obs['seurat_clusters'].isin(normal_fibro_clusters)
].copy()

# 2. Keloid fibroblasts (your clusters 0,3,4,10)
keloid_fibros = adata_fibroblasts

# 3. Differential expression: Keloid vs Normal fibroblasts
sc.tl.rank_genes_groups(adata_filtered, 
                       groupby='seurat_clusters',
                       groups=['0','3','4','10'],     # Keloid fibros
                       reference_groups=normal_fibro_clusters,  # Normal fibros
                       method='wilcoxon')

# 4. Top keloid-specific genes
keloid_de_genes = sc.get.rank_genes_groups_df(adata_filtered, group=['0','3','4','10'])
keloid_top50 = keloid_de_genes.nsmallest(50, 'logfoldchanges')['names'].tolist()

print("✅ TOP 50 KELIOD-SPECIFIC fibroblast genes:")
print(keloid_top50[:20])

# 5. Score keloid-specific signature
sc.tl.score_genes(adata_fibroblasts, gene_list=keloid_top50, score_name='KeloidSpecific')
sc.pl.umap(adata_fibroblasts, color='KeloidSpecific')


# In[ ]:


# Check if you have sample metadata
print("Available .obs columns:", adata_filtered.obs.columns.tolist())
print("\nSample Normal/Keloid column?", any('sample' in col.lower() or 'condition' in col.lower() or 'kd' in col.lower() for col in adata_filtered.obs.columns))

# Look for sample IDs (usually 'Kd1', 'Normal1', etc.)
print("\nUnique sample values:")
for col in adata_filtered.obs.columns:
    if adata_filtered.obs[col].nunique() < 20:  # Likely sample/condition
        print(f"{col}: {adata_filtered.obs[col].unique()[:5]}")


# In[ ]:


# Keloid fibroblasts only
keloid_fibro_mask = (
    (adata_filtered.obs['Group'] == 'Keloid') & 
    (adata_filtered.obs['seurat_clusters'].isin([0,3,4,10]))
)
keloid_fibros = adata_filtered[keloid_fibro_mask].copy()

# Normal fibroblasts only  
normal_fibro_mask = (
    (adata_filtered.obs['Group'] == 'Normal') & 
    (adata_filtered.obs['seurat_clusters'].isin([0,3,4,10]))
)
normal_fibros = adata_filtered[normal_fibro_mask].copy()

print(f"Keloid fibroblasts: {keloid_fibros.n_obs}")
print(f"Normal fibroblasts: {normal_fibros.n_obs}")


# In[ ]:


# Add condition to full dataset
adata_filtered.obs['fibro_condition'] = 'other'
adata_filtered.obs.loc[keloid_fibro_mask, 'fibro_condition'] = 'keloid_fibro'
adata_filtered.obs.loc[normal_fibro_mask, 'fibro_condition'] = 'normal_fibro'
adata_filtered.obs['fibro_condition'] = pd.Categorical(adata_filtered.obs['fibro_condition'])

# DE: Keloid fibroblasts vs Normal fibroblasts
sc.tl.rank_genes_groups(adata_filtered, 
                       groupby='fibro_condition',
                       groups=['keloid_fibro'], 
                       reference='normal_fibro',
                       method='wilcoxon')


# In[ ]:


keloid_de = sc.get.rank_genes_groups_df(adata_filtered, group='keloid_fibro')
keloid_top5 = keloid_de.head(100)['names'].tolist()

print("✅ TOP 5 KELIOD vs NORMAL FIBROBLAST genes:")
for i, gene in enumerate(keloid_top5, 1):
    print(f"{i}. {gene} (logFC: {keloid_de.iloc[i-1]['logfoldchanges']:.2f})")


# In[ ]:





# In[ ]:


# Get keloid fibroblast DE genes with logFC > 0.5
keloid_de_filtered = keloid_de[
    (keloid_de['logfoldchanges'] > 0.5) & 
    (keloid_de['pvals_adj'] < 0.05)
].head(20)  # Top 20 with logFC > 0.5

print("✅ TOP KELIOD FIBROBLAST GENES (logFC > 0.5, padj < 0.05):")
for i, (_, row) in enumerate(keloid_de_filtered.iterrows(), 1):
    print(f"{i:2d}. {row['names']:12s} logFC: {row['logfoldchanges']:6.2f} p-adj: {row['pvals_adj']:.2e}")

# Save all logFC > 0.5 genes
keloid_sig_genes = keloid_de[
    (keloid_de['logfoldchanges'] > 0.5) & 
    (keloid_de['pvals_adj'] < 0.05)
]['names'].tolist()

print(f"\nTotal keloid fibroblast genes (logFC > 0.5): {len(keloid_sig_genes)}")
print("Top 5:", keloid_sig_genes[:5])


# In[ ]:


# Built-in DE visualization
sc.pl.rank_genes_groups(adata_filtered, 
                       groups=['keloid_fibro'], 
                       n_genes=20, 
                       sharey=False, 
                       save='_keloid_fibro_de.png')


# In[ ]:


top_de_genes = keloid_de_filtered.head(20)['names'].tolist()
sc.pl.dotplot(adata_filtered, 
              top_de_genes, 
              groupby='fibro_condition', 
              swap_axes=True, 
              save='_keloid_de_dotplot.png')


# In[ ]:


all_sig_genes = keloid_de_filtered['names'].tolist()  # All your filtered genes
sc.pl.dotplot(adata_filtered, 
              all_sig_genes, 
              groupby='fibro_condition', 
              swap_axes=True, 
              standard_scale='var',
              save='_keloid_all_sig_dotplot.png')


# In[ ]:


# CASE = keloid fibroblasts
case_fibro_mask = (
    (adata_filtered.obs['Group'] == 'CASE') & 
    (adata_filtered.obs['seurat_clusters'].isin([0,3,4,10]))
)

# CTRL = normal fibroblasts  
ctrl_fibro_mask = (
    (adata_filtered.obs['Group'] == 'CTRL') & 
    (adata_filtered.obs['seurat_clusters'].isin([0,3,4,10]))
)

case_fibros = adata_filtered[case_fibro_mask].copy()
ctrl_fibros = adata_filtered[ctrl_fibro_mask].copy()

print(f"CASE (keloid) fibroblasts: {case_fibros.n_obs}")
print(f"CTRL (normal) fibroblasts: {ctrl_fibros.n_obs}")


# In[ ]:


# Keloid fibroblasts = CASE + fibro clusters
keloid_fibro_mask = (
    (adata_filtered.obs['Group'] == 'CASE') & 
    (adata_filtered.obs['seurat_clusters'].isin([0,3,4,10]))
)

# Normal fibroblasts = CTRL + fibro clusters
normal_fibro_mask = (
    (adata_filtered.obs['Group'] == 'CTRL') & 
    (adata_filtered.obs['seurat_clusters'].isin([0,3,4,10]))
)

print(f"Keloid fibroblasts: {adata_filtered[keloid_fibro_mask].n_obs}")
print(f"Normal fibroblasts: {adata_filtered[normal_fibro_mask].n_obs}")


# In[ ]:


# STEP 2: CREATE THE COLUMN (you missed this!)
keloid_fibro_mask = (
    (adata_filtered.obs['Group'] == 'CASE') & 
    (adata_filtered.obs['seurat_clusters'].isin([0,3,4,10]))
)
normal_fibro_mask = (
    (adata_filtered.obs['Group'] == 'CTRL') & 
    (adata_filtered.obs['seurat_clusters'].isin([0,3,4,10]))
)

# CREATE THE COLUMN
adata_filtered.obs['keloid_normal_fibro'] = 'non_fibro'
adata_filtered.obs.loc[keloid_fibro_mask, 'keloid_normal_fibro'] = 'keloid'
adata_filtered.obs.loc[normal_fibro_mask, 'keloid_normal_fibro'] = 'normal'
adata_filtered.obs['keloid_normal_fibro'] = pd.Categorical(adata_filtered.obs['keloid_normal_fibro'])

print("Column created!")
print(adata_filtered.obs['keloid_normal_fibro'].value_counts())


# In[ ]:


# CREATE 'keloid_normal_fibro' column (you skipped this!)
adata_filtered.obs['keloid_normal_fibro'] = 'non_fibro'
adata_filtered.obs.loc[keloid_fibro_mask, 'keloid_normal_fibro'] = 'keloid'
adata_filtered.obs.loc[normal_fibro_mask, 'keloid_normal_fibro'] = 'normal'
adata_filtered.obs['keloid_normal_fibro'] = pd.Categorical(adata_filtered.obs['keloid_normal_fibro'])

print("✅ Column created!")
print(adata_filtered.obs['keloid_normal_fibro'].value_counts())


# In[ ]:


# Save ALL DE genes
keloid_de.to_csv('keloid_vs_normal_fibroblast_DE_all.csv', index=False)

# Save ONLY significant genes (logFC > 0.5, padj < 0.05)
keloid_sig = keloid_de[
    (keloid_de['logfoldchanges'] > 0.5) & 
    (keloid_de['pvals_adj'] < 0.05)
]
keloid_sig.to_csv('keloid_vs_normal_fibroblast_DE_significant.csv', index=False)

print(f"✅ ALL DE genes saved: {len(keloid_de)} rows")
print(f"✅ Significant genes saved: {len(keloid_sig)} rows")


# In[ ]:


# Check what's in Group first
print(adata_fibroblasts.obs['Group'].value_counts())

# Subset to keloid fibroblasts only
adata_keloid = adata_fibroblasts[adata_fibroblasts.obs['Group'] == 'keloid'].copy()
# OR if it's in Condition_merged:
# adata_keloid = adata_fibroblasts[adata_fibroblasts.obs['Condition_merged'] == 'keloid'].copy()


# In[ ]:


# Remap CASE -> keloid, CTRL -> normal in Group column
adata_fibroblasts.obs['Group'] = adata_fibroblasts.obs['Group'].replace({
    'CASE': 'keloid',
    'CTRL': 'normal'
})

# Verify the new grouping
print(adata_fibroblasts.obs['Group'].value_counts())
# Now: Normal(10303) + CTRL(2758) = normal group
#      Keloid(9756) + CASE(2713) = keloid group


# In[ ]:


# Merge all keloid variants into one, all normal variants into one
adata_fibroblasts.obs['Group_clean'] = adata_fibroblasts.obs['Group'].replace({
    'keloid': 'keloid',    # lowercase -> keloid  
    'normal': 'Normal',    # lowercase -> normal
    'Keloid': 'keloid',    # titlecase -> keloid
    'Normal': 'Normal'     # titlecase -> normal
})

# Check clean groups
print(adata_fibroblasts.obs['Group_clean'].value_counts())
# Should show: normal (10303+2713=13016), keloid (9756+2758=12514)


# In[ ]:


# Keloid fibroblasts only (all keloid variants combined)
adata_keloid = adata_fibroblasts[adata_fibroblasts.obs['Group_clean'] == 'keloid'].copy()
print(f"Keloid fibroblasts: {adata_keloid.n_obs}")


# In[ ]:


postn_hi = adata_keloid[adata_keloid[:, 'POSTN'].X.toarray().flatten() > 0].copy()
mean_expr_postn = postn_hi.X.mean(axis=0)
if hasattr(mean_expr_postn, 'A1'):
    mean_expr_postn = mean_expr_postn.A1
postn_top50 = pd.Series(mean_expr_postn, index=postn_hi.var_names).nlargest(100).index.tolist()
sc.tl.score_genes(adata_keloid, gene_list=postn_top50, score_name='POSTNhi_new')

sc.pl.dotplot(adata_keloid, ['POSTNhi_new'], groupby='seurat_clusters', 
              swap_axes=True, save='_keloid_postnhi_dot.png')
sc.pl.umap(adata_keloid, color=['POSTNhi_new', 'seurat_clusters'], 
           legend_loc='on data', save='_keloid_postnhi_umap.png')


# In[ ]:


col1_expr = (adata_keloid[:, ['COL1A1', 'COL1A2']].X.toarray().sum(axis=1).flatten())
comp_hi = adata_keloid[col1_expr > np.quantile(col1_expr, 0.9)].copy()
mean_expr_comp = comp_hi.X.mean(axis=0)
if hasattr(mean_expr_comp, 'A1'):
    mean_expr_comp = mean_expr_comp.A1
comp_top50 = pd.Series(mean_expr_comp, index=comp_hi.var_names).nlargest(100).index.tolist()
sc.tl.score_genes(adata_keloid, gene_list=comp_top50, score_name='COMP_new')

sc.pl.dotplot(adata_keloid, ['COMP_new'], groupby='seurat_clusters', 
              swap_axes=True, save='_keloid_comp_dot.png')
sc.pl.umap(adata_keloid, color=['COMP_new', 'seurat_clusters'], 
           legend_loc='on data', save='_keloid_comp_umap.png')


# In[ ]:


# Print the top 50 genes for both signatures

# POSTN-hi top 50 genes
print("=== POSTN-hi Top 50 Genes ===")
print(postn_top50)

# COMP top 50 genes  
print("\n=== COMP Top 50 Genes ===")
print(comp_top50)

# Save as CSV for easy access
import pandas as pd
pd.DataFrame({'POSTNhi_genes': postn_top50}).to_csv('keloid_postnhi_top50.csv')
pd.DataFrame({'COMP_genes': comp_top50}).to_csv('keloid_comp_top50.csv')


# In[ ]:


print(adata_fibroblasts.obs['Group'].value_counts())
# Now: Normal(10303) + CTRL(2758) = normal group
#      Keloid(9756) + CASE(2713) = keloid group


# In[ ]:


# First, clean up the Group column to have ONLY 'normal' and 'keloid'
adata_fibroblasts.obs['keloid_normal_fibro'] = adata_fibroblasts.obs['Group'].replace({
    'Normal': 'normal',
    'Keloid': 'keloid', 
    'keloid': 'keloid',
    'normal': 'normal'
}).astype('category')

print(adata_fibroblasts.obs['keloid_normal_fibro'].value_counts())
# Should show clean: normal (~13k), keloid (~12k)


# In[ ]:


# DE: keloid vs normal fibroblasts
sc.tl.rank_genes_groups(adata_fibroblasts, 
                       groupby='keloid_normal_fibro',
                       groups=['keloid'], 
                       reference='normal',
                       method='wilcoxon')

keloid_de = sc.get.rank_genes_groups_df(adata_fibroblasts, group='keloid')
keloid_sig = keloid_de[
    (keloid_de['logfoldchanges'] > 0.5) & 
    (keloid_de['pvals_adj'] < 0.05)
]['names'].tolist()

print(f"✅ Keloid vs Normal fibroblast genes (logFC > 0.5): {len(keloid_sig)}")
print("Top 10:", keloid_sig[:10])

# Save
pd.DataFrame({'keloid_signature_genes': keloid_sig}).to_csv('keloid_vs_normal_sig.csv')


# In[ ]:


# 1. Save full DE results
keloid_de.to_csv('keloid_vs_normal_fibro_DE_full.csv', index=False)
keloid_sig_top50 = keloid_de.head(50)['names'].tolist()
pd.DataFrame({'top50_keloid_genes': keloid_sig_top50}).to_csv('keloid_top50_genes.csv')

print(f"✅ Saved full DE results and top 50 genes")
print("Top 50 keloid genes:", keloid_sig_top50)


# In[ ]:


# Perfect Scanpy built-in dotplot for top 50 keloid DE genes
sc.pl.dotplot(adata_fibroblasts, 
              keloid_sig_top50, 
              groupby='keloid_normal_fibro',
              swap_axes=False,  # genes on Y, groups on X
              standard_scale='var',  # scale expression 0-1 per gene
              save='_top50_keloid_vs_normal.png')

# Alternative: show logFC as color + % expressed as size
sc.pl.rank_genes_groups_dotplot(adata_fibroblasts, 
                               groupby='keloid_normal_fibro',
                               n_genes=50,
                               gene_symbols='names',
                               save='_top50_rank_dotplot.png')


# In[ ]:


# EBF2 Signature (Keloid Only)
# 1. EBF2-hi fibroblasts (high EBF2 expressors) 
ebf2_hi = adata_keloid[adata_keloid[:, 'EBF2'].X.toarray().flatten() > 0].copy()

# 2. Top 50 EBF2-hi signature (mean expression)
mean_expr_ebf2 = ebf2_hi.X.mean(axis=0)
if hasattr(mean_expr_ebf2, 'A1'):
    mean_expr_ebf2 = mean_expr_ebf2.A1
ebf2_top50 = pd.Series(mean_expr_ebf2, index=ebf2_hi.var_names).nlargest(50).index.tolist()

# 3. Score EBF2-hi signature across ALL keloid fibroblasts
sc.tl.score_genes(adata_keloid, gene_list=ebf2_top50, score_name='EBF2hi')

# 4. Visualize
sc.pl.dotplot(adata_keloid, ['EBF2hi'], groupby='seurat_clusters', 
              swap_axes=True, save='_keloid_ebf2_dot.png')
sc.pl.umap(adata_keloid, color=['EBF2hi', 'seurat_clusters'], 
           legend_loc='on data', save='_keloid_ebf2_umap.png')


# In[ ]:


# Print and save ASPN + EBF2 top 50 gene lists



print("\n=== EBF2-hi Top 50 Genes ===")
print(ebf2_top50)

# Save all top 50 signatures together
all_top50 = pd.DataFrame({
    'POSTNhi_top50': postn_top50,
    'COMP_top50': comp_top50, 
    'EBF2hi_top50': ebf2_top50
})

all_top50.to_csv('keloid_fibroblast_all_top50_signatures.csv', index=False)
print("\n✅ Saved all 4 top 50 signatures to 'keloid_fibroblast_all_top50_signatures.csv'")


# In[ ]:


# EBF2 Signature (ALL fibroblasts - normal + keloid)
# 1. EBF2-hi fibroblasts (high EBF2 expressors) across ALL fibroblasts
ebf2_hi_all = adata_fibroblasts[adata_fibroblasts[:, 'EBF2'].X.toarray().flatten() > 0].copy()

# 2. Top 50 EBF2-hi signature (mean expression)
mean_expr_ebf2_all = ebf2_hi_all.X.mean(axis=0)
if hasattr(mean_expr_ebf2_all, 'A1'):
    mean_expr_ebf2_all = mean_expr_ebf2_all.A1
ebf2_top50_all = pd.Series(mean_expr_ebf2_all, index=ebf2_hi_all.var_names).nlargest(50).index.tolist()

# 3. Score EBF2-hi signature across ALL fibroblasts
sc.tl.score_genes(adata_fibroblasts, gene_list=ebf2_top50_all, score_name='EBF2hi_all')

# 4. Visualize across all fibroblast clusters
sc.pl.dotplot(adata_fibroblasts, ['EBF2hi_all'], groupby='seurat_clusters', 
              swap_axes=True, save='_allfibro_ebf2_dot.png')
sc.pl.umap(adata_fibroblasts, color=['EBF2hi_all', 'seurat_clusters', 'keloid_normal_fibro'], 
           legend_loc='on data', save='_allfibro_ebf2_umap.png')

# Compare EBF2 across conditions
sc.pl.violin(adata_fibroblasts, ['EBF2hi_all'], groupby='keloid_normal_fibro', 
             save='_ebf2_by_condition.png')


# In[ ]:


print("=== EBF2-hi Top 50 Genes (ALL fibroblasts) ===")
print(ebf2_top50_all)

pd.DataFrame({'EBF2hi_all_top50': ebf2_top50_all}).to_csv('all_fibroblasts_ebf2_top50.csv')


# In[ ]:


# PPARGC1B Signature (Keloid Only) - same workflow
# 1. PPARGC1B-hi fibroblasts (high PPARGC1B expressors)
ppargc1b_hi = adata_keloid[adata_keloid[:, 'PPARGC1B'].X.toarray().flatten() > 0].copy()

# 2. Top 50 PPARGC1B-hi signature (mean expression)
mean_expr_ppargc1b = ppargc1b_hi.X.mean(axis=0)
if hasattr(mean_expr_ppargc1b, 'A1'):
    mean_expr_ppargc1b = mean_expr_ppargc1b.A1
ppargc1b_top50 = pd.Series(mean_expr_ppargc1b, index=ppargc1b_hi.var_names).nlargest(50).index.tolist()

# 3. Score PPARGC1B-hi signature across ALL keloid fibroblasts
sc.tl.score_genes(adata_keloid, gene_list=ppargc1b_top50, score_name='PPARGC1Bhi')

# 4. Visualize
sc.pl.dotplot(adata_keloid, ['PPARGC1Bhi'], groupby='seurat_clusters', 
              swap_axes=True, save='_keloid_ppargc1b_dot.png')
sc.pl.umap(adata_keloid, color=['PPARGC1Bhi', 'seurat_clusters'], 
           legend_loc='on data', save='_keloid_ppargc1b_umap.png')


# In[ ]:


print("=== PPARGC1B-hi Top 50 Genes ===")
print(ppargc1b_top50)

pd.DataFrame({'PPARGC1Bhi_top50': ppargc1b_top50}).to_csv('keloid_ppargc1b_top50.csv')


# In[ ]:


# Fixed mesenchymal adipocyte genes (core set - all should exist)
mesenchymal_adipo_sig = ['PPARG', 'CEBPA', 'FABP4', 'ADIPOQ', 'PLIN1', 'LPL', 'LEP','EBF2']

# Check which actually exist first
available_genes = [g for g in mesenchymal_adipo_sig if g in adata_keloid.var_names]
print("Available adipocyte genes:", available_genes)

# Dotplot with available genes only
sc.pl.dotplot(adata_keloid, available_genes, groupby='seurat_clusters', 
              save='_mesenchymal_adipo_genes.png')

# Score with available genes
sc.tl.score_genes(adata_keloid, gene_list=available_genes, score_name='Adipogenic')

sc.pl.umap(adata_keloid, color=['Adipogenic', 'seurat_clusters'], 
           save='_adipogenic_umap.png')

# Compare to your other signatures
sc.pl.dotplot(adata_keloid, ['POSTNhi_new', 'COMP_new', 'Adipogenic'], 
              groupby='seurat_clusters', swap_axes=True, 
              save='_adipogenic_vs_fibrotic.png')


# In[ ]:


# FABP4-hi Signature (Keloid Only) - same workflow as POSTN/COMP/ASPN
# 1. FABP4-hi fibroblasts (high FABP4 expressors)
fabp4_hi = adata_keloid[adata_keloid[:, 'FABP4'].X.toarray().flatten() > 0].copy()

# 2. Top 50 FABP4-hi signature (mean expression)
mean_expr_fabp4 = fabp4_hi.X.mean(axis=0)
if hasattr(mean_expr_fabp4, 'A1'):
    mean_expr_fabp4 = mean_expr_fabp4.A1
fabp4_top50 = pd.Series(mean_expr_fabp4, index=fabp4_hi.var_names).nlargest(50).index.tolist()

# 3. Score FABP4-hi signature across ALL keloid fibroblasts
sc.tl.score_genes(adata_keloid, gene_list=fabp4_top50, score_name='FABP4hi')

# 4. Visualize
sc.pl.dotplot(adata_keloid, ['FABP4hi'], groupby='seurat_clusters', 
              swap_axes=True, save='_keloid_fabp4_dot.png')
sc.pl.umap(adata_keloid, color=['FABP4hi', 'seurat_clusters'], 
           legend_loc='on data', save='_keloid_fabp4_umap.png')


# In[ ]:


print("=== FABP4-hi Top 50 Genes ===")
print(fabp4_top50)

pd.DataFrame({'FABP4hi_top50': fabp4_top50}).to_csv('keloid_fabp4_top50.csv')


# In[ ]:


# Mesenchymal Adipocyte-hi Signature (using available core genes)
# First check which adipocyte genes exist
adipo_genes = ['PPARG', 'CEBPA', 'FABP4', 'ADIPOQ', 'PLIN1', 'LPL', 'LEP']
available_adipo = [g for g in adipo_genes if g in adata_keloid.var_names]
print("Available mesenchymal adipocyte genes:", available_adipo)

# 1. Adipocyte-hi fibroblasts (high expression of available adipo genes)
adipo_expr = adata_keloid[:, available_adipo].X.mean(axis=1).toarray().flatten()
adipo_hi = adata_keloid[adipo_expr > np.quantile(adipo_expr, 0.9)].copy()  # Top 10%

# 2. Top 50 mesenchymal adipocyte-hi signature
mean_expr_adipo = adipo_hi.X.mean(axis=0)
if hasattr(mean_expr_adipo, 'A1'):
    mean_expr_adipo = mean_expr_adipo.A1
adipo_top50 = pd.Series(mean_expr_adipo, index=adipo_hi.var_names).nlargest(50).index.tolist()

# 3. Score across ALL keloid fibroblasts
sc.tl.score_genes(adata_keloid, gene_list=adipo_top50, score_name='Adipocyte_hi')

# 4. Visualize
sc.pl.dotplot(adata_keloid, ['Adipocyte_hi'], groupby='seurat_clusters', 
              swap_axes=True, save='_keloid_adipocyte_hi_dot.png')
sc.pl.umap(adata_keloid, color=['Adipocyte_hi', 'seurat_clusters'], 
           legend_loc='on data', save='_keloid_adipocyte_hi_umap.png')


# In[ ]:


# Mesenchymal Adipocyte-hi Signature (FIXED)
adipo_genes = ['PPARG', 'CEBPA', 'FABP4', 'ADIPOQ', 'PLIN1', 'LPL', 'LEP']
available_adipo = [g for g in adipo_genes if g in adata_keloid.var_names]
print("Available mesenchymal adipocyte genes:", available_adipo)

# 1. FIXED: Adipocyte-hi fibroblasts (top 10% expression)
adipo_expr = adata_keloid[:, available_adipo].X.mean(axis=1).A1  # Use .A1 for sparse!
adipo_hi = adata_keloid[adipo_expr > np.quantile(adipo_expr, 0.9)].copy()

# 2. Top 50 mesenchymal adipocyte-hi signature
mean_expr_adipo = adipo_hi.X.mean(axis=0).A1  # .A1 again!
adipo_top50 = pd.Series(mean_expr_adipo, index=adipo_hi.var_names).nlargest(100).index.tolist()

# 3. Score across ALL keloid fibroblasts
sc.tl.score_genes(adata_keloid, gene_list=adipo_top50, score_name='Adipocyte_hi')

# 4. Visualize
sc.pl.dotplot(adata_keloid, ['Adipocyte_hi'], groupby='seurat_clusters', 
              swap_axes=True, save='_keloid_adipocyte_hi_dot.png')
sc.pl.umap(adata_keloid, color=['Adipocyte_hi', 'seurat_clusters'], 
           save='_keloid_adipocyte_hi_umap.png')


# In[ ]:


print("=== Adipocyte-hi Top 50 ===")
print(adipo_top50)
pd.DataFrame({'Adipocyte_hi_top50': adipo_top50}).to_csv('keloid_adipocyte_hi_top50.csv')


# In[ ]:


# Just check chondrogenic gene expression (no signature needed)
chondro_genes = ['RUNX1', 'ACAN', 'OGN', 'BMP2', 'COMP']
available_chondro = [g for g in chondro_genes if g in adata_keloid.var_names]
print("Available chondrogenic genes:", available_chondro)

# Simple dotplot to see expression across clusters
sc.pl.dotplot(adata_keloid, available_chondro, groupby='seurat_clusters', 
              save='_chondro_expression.png')

# UMAP colored by each gene
sc.pl.umap(adata_keloid, color=available_chondro, 
           save='_chondro_umap.png')

# Violin plot by condition (if you want keloid vs normal comparison)
sc.pl.violin(adata_keloid, available_chondro, groupby='keloid_normal_fibro', 
             save='_chondro_violin.png')

