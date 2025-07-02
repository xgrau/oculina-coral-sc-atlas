# libraries
import sys
import numpy as np
import pandas as pd
import samap
from samap.utils import save_samap, load_samap, df_to_dict, substr
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from samalg import SAM
import scipy as sp

# output
out_fn = "results_samap/"

# list of multispecies comparisons: Oculina spp, all scleractinians, or all anthozoans
com_list = [ 
	{ "alpoculi" : ["Ocupat","Ocuarb"] },
	{ "alpscler" : ["Ocupat","Ocuarb","Spin","Amil"] },
	{ "alpantho" : ["Ocupat","Ocuarb","Spin","Amil","Nvec","Xesp"] },
	{ "subsym_Nv" : ["opasym","oarsym","Nvec"] },
	{ "subapo_Nv" : ["opaapo","oarapo","Nvec"] }
]

# list of datasets with sps id
dat_dict = {
	"Ocupat" : "Ocupat",
	"Ocuarb" : "Ocuarb",
	"Spin" : "Spin",
	"Amil" : "Amil",
	"Nvec" : "Nvec",
	"Xesp" : "Xesp",
	"opasym" : "Ocupat",
	"opaapo" : "Ocupat",
	"oarsym" : "Ocuarb",
	"oarapo" : "Ocuarb"
}


# functions
def q(x):
	return np.array(list(x))

def smod_get_mapping_scores_top_f(sm, keys, f_top = 0):
	"""Calculate mapping scores
	Parameters
	----------
	sm: SAMAP object

	keys: dict, annotation vector keys for at least two species with species identifiers as the keys
		e.g. {'pl':'tissue','sc':'tissue'}

	n_top: int, optional, default 0
		If `n_top` is 0, average the alignment scores for all cells in a pair of clusters.
		Otherwise, average the alignment scores of the top `n_top` cells in a pair of clusters.
		Set this to non-zero if you suspect there to be subpopulations of your cell types mapping
		to distinct cell types in the other species.
	Returns
	-------
	D - table of highest mapping scores for cell types 
	A - pairwise table of mapping scores between cell types across species
	"""


	if len(list(keys.keys()))<len(list(sm.sams.keys())):
		samap = SAM(counts = sm.samap.adata[np.in1d(sm.samap.adata.obs['species'],list(keys.keys()))])
	else:
		samap=sm.samap

	clusters = []
	ix = np.unique(samap.adata.obs['species'],return_index=True)[1]
	skeys = q(samap.adata.obs['species'])[np.sort(ix)]

	for sid in skeys:
		clusters.append(q([sid+'_'+str(x) for x in sm.sams[sid].adata.obs[keys[sid]]]))

	cl = np.concatenate(clusters)
	l = "{}_mapping_scores".format(';'.join([keys[sid] for sid in skeys]))
	samap.adata.obs[l] = pd.Categorical(cl)

	CSIMth, clu = smod_compute_csim_top_f(samap, l, f_top = f_top, prepend = False)

	A = pd.DataFrame(data=CSIMth, index=clu, columns=clu)
	i = np.argsort(-A.values.max(0).flatten())
	H = []
	C = []
	for I in range(A.shape[1]):
		x = A.iloc[:, i[I]].sort_values(ascending=False)
		H.append(np.vstack((x.index, x.values)).T)
		C.append(A.columns[i[I]])
		C.append(A.columns[i[I]])
	H = np.hstack(H)
	D = pd.DataFrame(data=H, columns=[C, ["Cluster","Alignment score"]*(H.shape[1]//2)])
	return D, A

def smod_compute_csim_top_f(samap, key, X=None, prepend=True, f_top = 0):
	splabels = q(samap.adata.obs['species'])
	skeys = splabels[np.sort(np.unique(splabels,return_index=True)[1])]

	cl = []
	clu = []
	for sid in skeys:
		if prepend:
			cl.append(sid+'_'+q(samap.adata.obs[key])[samap.adata.obs['species']==sid].astype('str').astype('object'))
		else:
			cl.append(q(samap.adata.obs[key])[samap.adata.obs['species']==sid])            
		clu.append(np.unique(cl[-1]))

	clu = np.concatenate(clu)
	cl = np.concatenate(cl)

	CSIM = np.zeros((clu.size, clu.size))
	if X is None:
		X = samap.adata.obsp["connectivities"].copy()

	xi,yi = X.nonzero()
	spxi = splabels[xi]
	spyi = splabels[yi]

	filt = spxi!=spyi
	di = X.data[filt]
	xi = xi[filt]
	yi = yi[filt]

	px,py = xi,cl[yi]
	p = px.astype('str').astype('object')+';'+py.astype('object')

	A = pd.DataFrame(data=np.vstack((p, di)).T, columns=["x", "y"])
	valdict = df_to_dict(A, key_key="x", val_key="y")   
	cell_scores = [valdict[k].sum() for k in valdict.keys()]
	ixer = pd.Series(data=np.arange(clu.size),index=clu)
	if len(valdict.keys())>0:
		xc,yc = substr(list(valdict.keys()),';')
		xc = xc.astype('int')
		yc=ixer[yc].values
		cell_cluster_scores = sp.sparse.coo_matrix((cell_scores,(xc,yc)),shape=(X.shape[0],clu.size)).A

		for i, c in enumerate(clu):
			if f_top > 0:
				n_top = int(len(np.sort(cell_cluster_scores[cl==c],axis=0)) * f_top)
				CSIM[i, :] = np.sort(cell_cluster_scores[cl==c],axis=0)[-n_top:].mean(0)
			else:
				CSIM[i, :] = cell_cluster_scores[cl==c].mean(0)

		CSIM = np.stack((CSIM,CSIM.T),axis=2).max(2)
		CSIMth = CSIM / samap.adata.uns['mapping_K']
		return CSIMth,clu
	else:
		return np.zeros((clu.size, clu.size)), clu


# loop
for n,com in enumerate(com_list):
	
	# query species data
	cid = [ c for c in com.keys() ][0]
	ddi_list = [ c for c in com.values() ][0]
	
	# get list for each sample id
	sps_list = []
	for ddi in ddi_list:
		spi = dat_dict[ddi]
		sps_list.append(spi)
		print("samap %s | data for %s will be mapped to %s" % (cid, ddi, spi))
	
	# cluster level dictionary (using sps id as keys)
	cts_list = [ "celltype" for s in ddi_list ]
	cts_dict = dict(zip(sps_list, cts_list))
	
	# load samap
	print("samap %s | load object" % cid)
	samm = load_samap("%s/samap.multi.%s.object" % (out_fn, cid))

	# broad cell types
	print("samap %s | alignment scores, cell types" % cid)
	map_topct,map_score = smod_get_mapping_scores_top_f(samm, keys = cts_dict, f_top = 0.9)
	map_score.to_csv("%s/samap.multi.%s.scores.cts.csv" % (out_fn, cid), sep = "\t")
	map_topct.to_csv("%s/samap.multi.%s.topcts.cts.csv" % (out_fn, cid), sep = "\t")
	
	# omit this, very slow?
	print("samap %s | alignment scores, gene pairs" % cid)
	gpf = samap.analysis.GenePairFinder(samm, keys = cts_dict)
	gene_pairs = gpf.find_all(align_thr = 0.05)
	gene_pairs.to_csv("%s/samap.multi.%s.scores.cts.shared_genes.csv" % (out_fn, cid), index = False)

	# end		
	print("samap %s | done" % cid)


print("all done!")
