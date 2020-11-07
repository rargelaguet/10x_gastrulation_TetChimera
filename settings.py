import numpy as np
import AnnData as anndata
import scipy as s
import scanpy as sc
import matplotlib.pyplot as pl
import seaborn
import pandas as pd
import csv
import os
from re import search

#########
## I/O ##
#########

io = {}

host = os.uname()[1]
if search("ricard", host):
	io["basedir"] = "/Users/ricard/data/10x_gastrulation_TetChimera"
	io["atlas.basedir"] = "/Users/ricard/data/gastrulation10x"
	io["gene_metadata"] = "/Users/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
elif search("ebi", host):
	io["basedir"] = "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera"
	io["atlas.basedir"] = "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
	io["gene_metadata"] = "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
else:
	print("Computer not recognised"); exit()

io["metadata"] = io["basedir"] + "/sample_metadata.txt.gz"

# Atlas information
io["atlas.metadata"] = io["atlas.basedir"] + "/sample_metadata.txt.gz"


#############
## Options ##
#############

opts = {}

opts["celltypes"] = [
	"Epiblast",
	"Primitive_Streak",
	"Caudal_epiblast",
	"PGC",
	"Anterior_Primitive_Streak",
	"Notochord",
	"Def._endoderm",
	"Gut",
	"Nascent_mesoderm",
	"Mixed_mesoderm",
	"Intermediate_mesoderm",
	"Caudal_Mesoderm",
	"Paraxial_mesoderm",
	"Somitic_mesoderm",
	"Pharyngeal_mesoderm",
	"Cardiomyocytes",
	"Allantois",
	"ExE_mesoderm",
	"Mesenchyme",
	"Haematoendothelial_progenitors",
	"Endothelium",
	"Blood_progenitors_1",
	"Blood_progenitors_2",
	"Erythroid1",
	"Erythroid2",
	"Erythroid3",
	"NMP",
	"Rostral_neurectoderm",
	"Caudal_neurectoderm",
	"Neural_crest",
	"Forebrain_Midbrain_Hindbrain",
	"Spinal_cord",
	"Surface_ectoderm",
	"Visceral_endoderm",
	"ExE_endoderm",
	"ExE_ectoderm",
	"Parietal_endoderm"
]

opts["celltype.colors"] = {
	"Epiblast" : "#635547",
	"Primitive_Streak" : "#DABE99",
	"Caudal_epiblast" : "#9e6762",
	"PGC" : "#FACB12",
	"Anterior_Primitive_Streak" : "#c19f70",
	"Notochord" : "#0F4A9C",
	"Def._endoderm" : "#F397C0",
	"Gut" : "#EF5A9D",
	"Nascent_mesoderm" : "#C594BF",
	"Mixed_mesoderm" : "#DFCDE4",
	"Intermediate_mesoderm" : "#139992",
	"Caudal_Mesoderm" : "#3F84AA",
	"Paraxial_mesoderm" : "#8DB5CE",
	"Somitic_mesoderm" : "#005579",
	"Pharyngeal_mesoderm" : "#C9EBFB",
	"Cardiomyocytes" : "#B51D8D",
	"Allantois" : "#532C8A",
	"ExE_mesoderm" : "#8870ad",
	"Mesenchyme" : "#cc7818",
	"Haematoendothelial_progenitors" : "#FBBE92",
	"Endothelium" : "#ff891c",
	"Blood_progenitors" : "#c9a997",
	"Blood_progenitors_1" : "#f9decf",
	"Blood_progenitors_2" : "#c9a997",
	"Erythroid" : "#EF4E22",
	"Erythroid1" : "#C72228",
	"Erythroid2" : "#f79083",
	"Erythroid3" : "#EF4E22",
	"NMP" : "#8EC792",
	"Neurectoderm" : "#65A83E",
	"Rostral_neurectoderm" : "#65A83E",
	"Caudal_neurectoderm" : "#354E23",
	"Neural_crest" : "#C3C388",
	"Forebrain_Midbrain_Hindbrain" : "#647a4f",
	"Spinal_cord" : "#CDE088",
	"Surface_ectoderm" : "#f7f79e",
	"Visceral_endoderm" : "#F6BFCB",
	"ExE_endoderm" : "#7F6874",
	"ExE_ectoderm" : "#989898",
	"Parietal_endoderm" : "#1A1A1A"
}

opts["stages_colors"] = {
	'E6.5':"#D53E4F",
	'E6.75':"#F46D43",
	'E7.0':"#FDAE61",
	'E7.5':"#FFFFBF",
	'E7.25':"#FEE08B",
	'E7.75':"#E6F598",
	'E8.0':"#ABDDA4",
	'E8.5':"#3288BD",
	'E8.25':"#66C2A5",
	'mixed_gastrulation': "#A9A9A9"  
	
}


opts["batches"] = [

	# First batches (all failed QC)
	# SIGAA3_E8.5_pool1_Host-WT_L001
	# SIGAB3_E8.5_pool1_TET-TKO_L002
	# SIGAC3_E8.5_pool2_Host-WT_L003
	# SIGAD3_E8.5_pool2_TET-TKO_L004
	# SIGAE3_E7.5_pool1_Host-WT_L005
	# SIGAF3_E7.5_pool1_TET-TKO_L006
	# SIGAG3_E8.5_hashing_Host-WT_L007
	# SIGAH3_E8.5_hasting_TET-TKO_L00

	# Second batch
	"E75_TET_TKO_L002",
	"E75_WT_Host_L001",
	"E85_Rep1_TET_TKO_L004",
	"E85_Rep1_WT_Host_L003",
	"E85_Rep2_TET_TKO_L006",
	"E85_Rep2_WT_Host_L005"

	# Third batch (failed QC)
	# "SIGAE4_E105_3_TET123_Chimera_Host_L005", 
	# "SIGAF4_E105_3_TET123_Chimera_TKO_L006", 
	# "SIGAG4_E105_5_TET123_Chimera_Host_L007", 
	# "SIGAH4_E105_5_TET123_Chimera_TKO_L008"

	# Fourth batch
	# "SIGAC2_TET_TKO_E9_5_Head1_L002",
	# "SIGAD2_TET_TKO_E9_5_Trunk1_L002",
	# "SIGAE2_TET_TKO_E9_5_Tail1_L002",
	# "SIGAE6_TET_TKO_E9_5_Head2_L003",
	# "SIGAF2_TET_TKO_E9_5_YS1_L002",
	# "SIGAF6_TET_TKO_E9_5_Trunk2_L003",
	# "SIGAG6_TET_TKO_E9_5_Tail2_L003",
	# "SIGAH6_TET_TKO_E9_5_YS2_L003"
]

opts["classes"] = [
	"E7.5_Host", 
	"E7.5_TET_TKO", 
	"E8.5_Host", 
	"E8.5_TET_TKO"
	# "E9.5_TET_TKO",
	# "E10.5_Host", 
	# "E10.5_TET_TKO"
]

##########################
## Load sample metadata ##
##########################