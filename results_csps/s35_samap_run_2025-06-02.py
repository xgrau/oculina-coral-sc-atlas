# libraries
import sys
import numpy as np
import pandas as pd
import samap
import scanpy as sc
import samalg
from samap.mapping import SAMAP
from samap.utils import save_samap, load_samap

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
	
	# sample dictionary (using sps id as keys)
	sam_dict = dict(zip(sps_list, [None] * len(sps_list)))

	# cluster level dictionary (using sps id as keys)
	cts_list = [ "celltype" for s in ddi_list ]
	cts_dict = dict(zip(sps_list, cts_list))

	# load data into sample dictionary
	for ddi in ddi_list:
		spi = dat_dict[ddi]
		print("samap %s | loading %s data into %s slot" % (cid, ddi, spi))
		sam_dict[spi] = samalg.SAM()
		sam_dict[spi].load_data("%s/sam.%s.object.h5ad" % (out_fn, ddi))
		sc.pp.highly_variable_genes(sam_dict[spi].adata)

	# run samap, all homlogs with cell type clusters
	print("samap %s | init SAMAP..." % (cid))
	samm = SAMAP(sam_dict, f_maps = "results_samap/blasts/", keys = cts_dict)
	print("samap %s | run SAMAP..." % (cid))
	samm.run()
	print("samap %s | save SAMAP..." % (cid))
	save_samap(samm, "%s/samap.multi.%s.object" % (out_fn, cid))
	print("samap %s | done!" % (cid))

print("all done!")
