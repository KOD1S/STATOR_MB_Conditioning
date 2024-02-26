import pandas as pd
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
from matplotlib import cm
import argparse
import seaborn as sns
viridis = cm.get_cmap('viridis', 12)
sns.set_style("white")
sns.set_palette('colorblind')
from utilities import *
import collections
from collections import OrderedDict
import itertools
import random
import hypernetx as hnx
import scanpy as sc
from scipy.stats import binom_test
from statsmodels.stats.multitest import multipletests
import upsetplot as up
from upsetplot import generate_counts
from upsetplot import from_indicators
import io
from PIL import Image
print('Modules imported \n')

dataPath = "./trainingData_10000Cells_1000Genes.csv"
MCMCgraphPath = "./MCMCgraph_10000Cells_1000Genes.csv"
CPDAGgraphPath = "./Data/CPDAGgraph_10000Cells_1000Genes.csv"

## Set filepaths of interaction files (STATOR outputs these into the NF_TL_pipeline/coupling_output folder)
pathTo3pts = "./Outputs/interactions_withinMB_3pts_MCMCgraph_10000Cells_1000Genes_expectations.npy"
pathTo4pts = "./Outputs/interactions_withinMB_4pts_MCMCgraph_10000Cells_1000Genes_expectations.npy"
pathTo5pts = "./Outputs/interactions_withinMB_5pts_MCMCgraph_10000Cells_1000Genes_expectations.npy"

# Choose an identifier for this MB condition (for example the gene name you are setting to 1)
File_ID = 'X'

minStateDeviation = 0.001
stateDevAlpha = 0.05
sigHOIthreshold = 0.05
trainDat = pd.read_csv(dataPath)
DSname = MCMCgraphPath.split('.')[0]
MCMCadjMat = pd.read_csv(MCMCgraphPath, index_col=0)
MCMCgraph = ig.Graph.Adjacency(MCMCadjMat.values.tolist())
MCMCgraph.vs['label'] = trainDat.columns.values
CPDAGadjMat = pd.read_csv(CPDAGgraphPath, index_col=0)
CPDAGgraph = ig.Graph.Adjacency(CPDAGadjMat.values.tolist())
CPDAGgraph.vs['label'] = trainDat.columns.values
genes = trainDat.columns.values
scObj = sc.AnnData(trainDat) # Data is loaded as scanpy Anndata objects for easy embedding plots. 

# HHOIs stores all the significant (n>1)-point interactions present among interacting triplets, quadruplets, and pentuplets
HHOIs = {}
alpha = sigHOIthreshold
for order, intPath in enumerate([pathTo3pts, pathTo4pts, pathTo5pts]):
	ints = np.load(intPath, allow_pickle=True)
	if len(ints)>0:
		perfectSigEsts = list(map(lambda x: (((x[[4, 5, 6]]==0).all()) & (x[3]<alpha)), ints))
		HHOIs[f'n{order+3}'] = ints[perfectSigEsts][:, [0, -1]]
	else: HHOIs[f'n{order+3}'] = []

def findsubsets(s, n):
	return list(itertools.combinations(s, n))

f=2
kwargs = {'with_node_counts': False, 'with_node_labels':False, 'with_edge_labels':False}
concatInts = lambda x: ''.join(map(str, x))
ordersToPlot = [3, 4, 5]
enrichments = {}

for order in ordersToPlot:
	nStates = 2**order
	# Binary representation of each possible state:
	binStates = [np.array(list(format(x, f"0{order}b"))).astype(bool) for x in range(nStates)]
	# Stores the log 2-fold change of each d-tuple
	log2Fs = []
	if len(HHOIs[f'n{order}'])>0:
		# loop over the HOIs, use their value (weight w) and gene tuple:
		for w, geneTuple in HHOIs[f'n{order}'][:, [0, -1]]:
			ID = '_'.join(genes[geneTuple])
			# The local CPDAG structure will be plotted, and its layout will be used for the hypergraph of interactions. 
			g = findLocalGraph(geneTuple, CPDAGgraph, order=0)
			layout_c = g.layout('circle')
			# Layout needs to be manipulated a bit so that hypernetx can use it for the hypergraph:
			# tmp = dict(zip([x.replace('.', '-') for x in g.vs['label']], np.array(layout_c)*np.array([1, -1])))
			# NOTE: I have removed the string replacement, as it no longer seems to be necessary and was causing problems. 
			tmp = dict(zip([x for x in g.vs['label']], np.array(layout_c)*np.array([1, -1])))
			# Constructing the hypergraph
			edges = {'maxOrder': genes[geneTuple]}
			weights = [w] # Stores the actual interactions -- used to colour the edges in the hypergraph. ########### Changed This to 3, order ##############
			for i in range(3, order):
				for subset in findsubsets(geneTuple, i):
					for tup in HHOIs[f'n{i}']:
						if sorted(tup[-1]) == sorted(subset):
							edges[len(edges)+1] =  [genes[g] for g in subset]
							weights.append(tup[0])
							break
			#  ************************ Set Conditioned Genes ************************ 
			unConditionedGenes = trainDat.iloc[:, geneTuple]
			conditionedGenes = conditionOnMB(geneTuple, MCMCgraph, trainDat, mode='Min') # required
			#  ************************ Calculate deviations ************************ 
			# Count the occurence of each of the possible binary states in the conditioned data:
			binCounts = np.bincount(list(map(lambda x: int(x, 2), list(map(concatInts, conditionedGenes.values)))), minlength=nStates) # requuired
			# The means determine the null hypothesis
			means = conditionedGenes.mean(axis=0)
			pStates = np.array([np.prod([m if state[i] else 1-m for i, m in enumerate(means)]) for state in binStates])
			expected = pStates*len(conditionedGenes)
			# Calculate the log 2-fold enrichment, and p-values based on the binomial null
			# Vectorised over the 2^n states
			log2FoldEnrichment = np.log2(binCounts/expected)
			enrichment_pval = np.array([binom_test(binCounts[i], p = pStates[i], n = len(conditionedGenes), alternative='greater') for i in range(len(binCounts))])
			log2Fs.append([log2FoldEnrichment, enrichment_pval, geneTuple])
		log2Fs  = np.array(log2Fs, dtype=object)
		enrichments[f'n{order}'] = log2Fs
	else: enrichments[f'n{order}'] = []

# Construct a dataframe with the deviations of each positively enriched state:
devDict = []
for order in [3, 4, 5]:
	for devs, pvals, interactors in enrichments[f'n{order}']:
		ID = '_'.join(genes[interactors])
		for devStateInd in np.where((devs>=0))[0]:
			# Format the binary state
			maxDevState = format(devStateInd, f"0{order}b")
			devDict.append([ID, maxDevState, devs[devStateInd], pvals[devStateInd]])

if len(devDict)>0:
	deviators = pd.DataFrame(devDict, columns=['genes', 'state', 'enrichment', 'pval'])
	deviators = deviators.sort_values(by='pval', ascending=True)
	# Add the Benjamini-Yekutieli corrected p-values:
	deviators['pval_corrected'] = multipletests(deviators['pval'], method='fdr_by')[1]
	# Binreps is the binary represenation of the interactions: binReps[i] is 1 if cell i is in the deviating state, 0 otherwise.  
	binReps = np.array(deviators.apply(lambda x: (trainDat[x['genes'].rsplit('_')]==[int(g) for g in list(str(x['state']))]).all(axis=1), axis=1))*1
	# Add the cell IDs to each of the d-tuples:
	# (has to be made into a list to avoid the ellipsis in the output)
	deviators['cellIDs'] = [np.where(binRep)[0].astype(int).tolist() for binRep in binReps]
	# Filter out the d-tuples that are significantly enriched:
	strongDeviators = deviators.query(f'pval_corrected <= {stateDevAlpha} and enrichment >= {minStateDeviation}').copy()
	deviators.to_csv(f'all_DTuples_' + File_ID + '.csv')
	strongDeviators.to_csv(f'top_DTuples_' + File_ID + '.csv')
	pd.DataFrame(binReps).to_csv(f'DTuples_binaryReps_' + File_ID + '.csv')
else:
	pd.DataFrame(data=[]).to_csv(f'all_DTuples_' + File_ID + '.csv')
	pd.DataFrame(data=[]).to_csv(f'top_DTuples_' + File_ID + '.csv')
