## Introduction

I’ll walk you through a step-by-step single-cell RNA sequencing (scRNA-seq) analysis pipeline using Scanpy in Python. I’ll explain each step, share my thought process, and provide code snippets that you can run yourself.

## Step 1: Setting Up the Environment

Before we start, we need to set up a clean environment to ensure everything runs smoothly. Here’s how I did it:

```python
# Create a virtual environment
python -m venv scRNA_env
source scRNA_env/bin/activate  # On Windows: scRNA_env\Scripts\activate

# Install required libraries
pip install scanpy numpy pandas matplotlib seaborn scikit-learn

```

---


## Step 2: Loading the Data

The first step is to load the data. For 10X Genomics data, we use the `read_10x_mtx` function from Scanpy. Here’s how I did it:

```python
import scanpy as sc

# Load 10X data
adata = sc.read_10x_mtx(
    "path/to/matrix/folder",  # Replace with your data path
    var_names='gene_symbols',  # Use gene symbols (or 'gene_ids' for Ensembl IDs)
    cache=True                # Cache for faster reloading
)
print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes.")
```

---

## Step 3: Quality Control (QC)

Single-cell data can be noisy, so we need to perform quality control to filter out low-quality cells and genes. Here’s what I did:

```python
# Filter cells and genes
sc.pp.filter_cells(adata, min_genes=200)  # Keep cells with ≥200 genes
sc.pp.filter_genes(adata, min_cells=3)    # Keep genes detected in ≥3 cells

# Annotate mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Adjust prefix if needed
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Visualize QC metrics
sc.pl.violin(adata, ['n_genes', 'n_counts', 'pct_counts_mt'], save='_qc.png')
```

---

## Step 4: Normalization and Log Transformation

Next, we normalize the data to account for differences in sequencing depth and apply a log transformation to stabilize variance:

```python
# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)  # Normalize to 10,000 counts per cell
sc.pp.log1p(adata)  # Log-transform the data
```

---

## Step 5: Feature Selection

To reduce noise, we focus on the **highly variable genes (HVGs)** that show the most variation across cells:

```python
# Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]  # Subset to HVGs
```

---

## Step 6: Dimensionality Reduction

To visualize and cluster the data, we reduce its dimensionality using PCA and UMAP:

```python
# Scale the data
sc.pp.scale(adata, max_value=10)

# Perform PCA
sc.tl.pca(adata, svd_solver='arpack')

# Compute UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
```

---

## Step 7: Clustering

We cluster the cells using the Leiden algorithm to identify distinct cell populations:

```python
# Cluster cells
sc.tl.leiden(adata)
```

---

## Step 8: Visualization

Finally, we visualize the results to interpret the clusters and identify marker genes:

```python
# Plot UMAP
sc.pl.umap(adata, color=['leiden'], save='_clusters.png')

# Identify marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_markers.png')
```

---

## Step 9: Saving the Results

To make sure we can revisit the analysis later, we save the processed data and results:

```python
# Save the AnnData object
adata.write('processed_data.h5ad')

# Export marker genes to CSV
markers = sc.get.rank_genes_groups_df(adata, group='0')
markers.to_csv('marker_genes_cluster_0.csv')
```

---

## Conclusion

That’s it! This pipeline takes raw scRNA-seq data and processes it to identify cell types, visualize clusters, and find marker genes.

