#!/usr/bin/env python3
# Copyright (c) [2021] [Monica T. Hannani]
# mhannani@ukaachen.de


''' Compositional data analysis with scCODA'''


# Initialize
import scanpy as sc
import numpy as numpy
import pandas as pd
import matplotlib.pyplot as plt

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz


# Read cell counts
indir = '/Users/monica/Dropbox/UKA/marmoset/'
cell_counts = pd.read_csv(indir + 'sample_cell_type_counts.csv')


# Convert data to anndata object and add condition
data = dat.from_pandas(cell_counts, covariate_columns = ['sample_condition'])
data.obs['condition'] = data.obs['sample_condition'].str.replace(r'\d+_', '')


# Compare compositional changes between young and old marmoset kidneys
case = 'old'
control = 'young'
reference_cell_type = 'DCT'


# Subset
data_cond = data[data.obs['condition'].isin([case, control])]
data_cond.obs['condition'] = data_cond.obs['condition'].astype(str)
data_cond.obs['condition'] = data_cond.obs['condition'].astype('category').cat.reorder_categories([control, case])


# Visualize cell type proportions
# viz.boxplots(data_cond, feature_name = 'condition')


# Setup model and inference
# Use condition as covariate in model
# Use DCT as reference cell type (they vary very little in all samples and show low difference in condition)
reference_cell_type = 'DCT'
model_cond = mod.CompositionalAnalysis(data_cond, formula = 'condition', reference_cell_type = reference_cell_type)


# Run Hamiltonian Monte Carlo (HMC) sampling for parameter inference (Markov chain Monte Carlo, MCMC)
sim_results = model_cond.sample_hmc()


# Save results
sim_results.save(indir + 'CoDa_young_vs_old')


