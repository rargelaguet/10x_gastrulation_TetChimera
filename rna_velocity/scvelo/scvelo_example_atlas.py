###
### scvelo example with the atlas

import numpy as np
import scvelo as scv
import scanpy as sc
import matplotlib.pyplot as plt

# Assuming the object has batch corrected pcs and atlas umap
adata = sc.read('adata_umap_pca.h5',cache=True)

scv.pp.filter_and_normalize(adata, n_top_genes=5000, min_shared_counts=20)

sc.pp.neighbors(adata, n_neighbors = 30, n_pcs = 50)
sc.tl.louvain(adata, flavor='igraph')

scv.pp.moments(adata, n_pcs=50, n_neighbors=30)

scv.tl.recover_dynamics(adata) # required for dynamical mode (it uses both, PCA and UMAP)

scv.tl.velocity(adata, mode="dynamical") #add time covariate in the future

scv.tl.velocity_graph(adata)

# set approx True to use batch corrected PCs
scv.tl.velocity_graph(adata, approx = True, n_neighbors = 30)

# Scale colours found in the atlas github
colPalette_celltype = [celltype_colours[i] for i in sorted(np.unique(adata_proc.obs['celltype']))]
adata_proc.uns['celltype_colors'] = colPalette_celltype

colPalette_stage = [stage_colours_wot[i] for i in sorted(np.unique(adata_proc.obs['stage']))]
adata_proc.uns['stage_colors'] = colPalette_stage

fig, axs = plt.subplots(1,1, figsize = (20,20))

scv.pl.velocity_embedding_stream(
    adata, basis="umap", color="celltype", size=50, ax = axs)

fig, axs = plt.subplots(1,1, figsize = (20,20))

# Latent time recovery
scv.tl.velocity_pseudotime(adata)

scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')

# Confidence stats (see Speed and coherence)
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])

# Get top likelihood genes 
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False)

# Cluster-specific top-likelihood genes
scv.tl.rank_dynamical_genes(adata, groupby='louvain') # could be celltypes

# Differential kinetics
var_names = ['Gata1', 'Etv2']
scv.tl.differential_kinetic_test(adata, var_names=var_names, groupby='louvain')
# Here the example in the documentation is more complete

# You can use PAGA with scVelo too but I have not used it in the atlas

