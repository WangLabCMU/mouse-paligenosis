#Load packages
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
scv.settings.verbosity = 3  
scv.settings.presenter_view = True  
scv.set_figure_params('scvelo')  

# Prepare input data
adata_h5ad = sc.read("./adata.h5ad", cache=True)
adata_loom = scv.read("./adata.loom", cache=True)
adata_h5ad.obs.index = [element[0] for element in adata_h5ad.obs.index.str.split("-")]
adata_loom.obs.index = [element[1] for element in adata_loom.obs.index.str.split(":")]
adata_m = scv.utils.merge(adata_h5ad, adata_loom, copy=True)

# Preprocesses, normalizes, and analyzes data, inferring cell dynamics and trajectories, and computing velocity-based pseudotime and latent time
scv.pp.filter_and_normalize(adata_m, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_m, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata_m, n_jobs=30)                                            
scv.tl.velocity(adata_m, mode='dynamical')
scv.tl.velocity_graph(adata_m, n_jobs=30)
scv.tl.terminal_states(adata_m)
scv.tl.velocity_pseudotime(adata_m)
scv.tl.latent_time(adata_m)

# Visualize RNA velocity plot
scv.pl.proportions(adata_m, groupby="leiden", figsize=(6, 2))
scv.pl.velocity_embedding_stream(adata_m, basis='umap', color="leiden", dpi=100)
scv.pl.velocity_embedding(adata_m, basis='umap', arrow_length=3, arrow_size=2, color="leiden", dpi=105)
scv.pl.scatter(adata_m, color=[ 'root_cells', 'end_points'], basis='umap') 
scv.pl.scatter(adata_m, color='velocity_pseudotime', cmap='gnuplot', basis="umap")
scv.pl.scatter(adata_m, color='latent_time', color_map='gnuplot', size=80)
