import anndata as anndata
import scanpy as sc

def load_adata(adata_file, metadata_file = None, normalise = False, cells = None, features = None, filter_lowly_expressed_genes = False, set_colors = True):

	adata = sc.read(adata_file)

	if cells is not None:
		# np.in1d(cells, adata.obs.index) # super slow
		adata = adata[cells,:]

	if features is not None:
		adata = adata[:,features]

	if metadata_file is not None:
		metadata = pd.read_table(metadata_file, delimiter="\t", header=0).set_index("cell", drop=False)
		metadata = metadata.loc[cells]
		assert np.all(adata.obs.index.isin(metadata.cell))
		# assert np.all(metadata.cell.isin(adata.obs.index))
		assert metadata.shape[0] == adata.shape[0]
		adata.obs = metadata#.reindex(adata.obs.index)

	if filter_lowly_expressed_genes:
		sc.pp.filter_genes(adata, min_counts=10)

	if normalise:
		sc.pp.normalize_total(adata, target_sum=None, exclude_highly_expressed=False)
		sc.pp.log1p(adata)

	if set_colors:
		colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_").replace("/","_")] for i in sorted(np.unique(adata.obs['celltype.mapped']))]
		adata.uns['celltype.mapped_colors'] = colPalette_celltypes
		colPalette_stages = [opts["stages_colors"][i.replace(" ","_").replace("/","_")] for i in sorted(np.unique(adata.obs['stage']))]
		adata.uns['stage_colors'] = colPalette_stages

	return adata
