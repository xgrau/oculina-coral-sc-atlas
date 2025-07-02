# libraries
import pandas as pd
import numpy as np
import sys
import samalg
import scanpy as sc

# output
out_fn = "results_samap/"

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

# query species
for ddi in dat_dict.keys():
	
	# read
	spi = dat_dict[ddi]
	print("# loading %s" % (ddi))
	mat = sc.read_mtx("results_samap/data/sam.%s.counts.mtx" % ddi)
	sct = pd.read_csv("results_samap/data/sam.%s.cell_to_cts.csv" % ddi, sep="\t", header = None, names = ["cell","celltype"])
	gen = pd.read_csv("results_samap/data/sam.%s.genes.txt" % ddi, sep="\t", names = ["gene"])

	# create SAM objects for concatenated dataset
	print("# run SAM %s" % (ddi))
	sam = samalg.SAM(counts=[mat.X.transpose(), gen["gene"].values,sct["cell"].values ])
	sam.preprocess_data()
	for i in range(sct.shape[1]):
		sam.adata.obs[sct.columns[i]] = sct.iloc[:,i].values
	sam.run()

	# save precomputed sam object
	sam.save_anndata("%s/sam.%s.object.h5ad" % (out_fn, ddi))

print("all done!")
