# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run Scrublet doublet scoring for primary snRNA-seq libraries.
# Steps: 1) Read per-library matrices; 2) calculate scores and embeddings; 3) write doublet diagnostics for QC filtering.
# =======================================
# Source reformatted from snRNA_AgeEx/MarkdownRNAQC5.Rmd.
# Set the variables in config/.env before running this script.
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.paths import publication_path

import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import gzip
import scrublet as scr
import pandas as pd 
#Function to get gene list 
def load_genes(filename, delimiter='\t', column=0, skip_rows=0):
    gene_list = []
    gene_dict = {}

    with gzip.open(filename, 'rb') as f:
        for iL in range(skip_rows):
            f.readline()
        for l in f:
            gene = l.decode().strip('\n').split(delimiter)[column]
            if gene in gene_dict:
                gene_dict[gene] += 1
                gene_list.append(gene + '__' + str(gene_dict[gene]))
                if gene_dict[gene] == 2:
                    i = gene_list.index(gene)
                    gene_list[i] = gene + '__1'
            else: 
                gene_dict[gene] = 1
                gene_list.append(gene)
    return gene_list

#Set plot variables 
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

#Define Runs 
Runs = ["MAE2_1","MAE4_3","MAS1_4","MAS3_2","YE2_1","YE4_3","YS1_4","YS3_2"]


#Loop Thorugh Runs running scrublet 
for Run in Runs:
    #read in 
    print("run ="+Run)
    # publication_path() returns str(Path(...)), which drops the trailing slash; the paths below are built by
    # string concatenation (out_dir + Run + ..., input_dir + folder_dir), so the separator is added explicitly.
    out_dir = publication_path("HEART_SNRNA_ROOT", "data/2022-07-04.Scrublet") + "/"
    input_dir = publication_path("HEART_SNRNA_ROOT", "data/CellRanger") + "/"
    folder_dir = Run +'/outs/filtered_feature_bc_matrix/'
    counts_matrix = scipy.io.mmread(input_dir + folder_dir + 'matrix.mtx.gz').T.tocsc()
    counts_matrix = counts_matrix[:,:-9]
    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    genes = np.array(load_genes(input_dir + folder_dir + 'features.tsv.gz', delimiter='\t', column=1))
    genes = genes[:-9]
    print('Number of genes in gene list: {}'.format(len(genes)))
    #Initialize Scrublet object
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    #

    #Run the default pipeline
    #Doublet simulation
    #Normalization, gene filtering, rescaling, PCA
    #Doublet score calculation
    #Doublet score threshold detection and doublet calling

    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
    scrub.plot_histogram()
    #
    #
    scrub.plot_histogram()[0].savefig(out_dir + Run +'_scrublet_scores.pdf')
    #
    #
    print(doublet_scores[0:5])
    print(doublet_scores.shape)
    pd.DataFrame(doublet_scores).to_csv(out_dir + Run +'_doublet_scores.txt',sep='\t')
    #
    #

    #Get 2-D embedding to visualize the results
    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    print('Done.')

    #Plot doublet predictions on 2-D embedding
    scrub.plot_embedding('UMAP', order_points=True)[0].savefig(out_dir + Run + '_doublet_UMAP.pdf')
    print(Run + " completed")

