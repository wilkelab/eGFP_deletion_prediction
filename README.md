This is the main README file for article "Computational prediction of the tolerance to amino-acid deletion in green-fluorescent protein" authored by Eleisha L. Jackson, Stephanie J. Spielman and Claus O. Wilke

The directory contains all of the associated data and scripts for this paper. 

Contents:

data/ 
    This directory contains all features for each deletion analyzed in the paper.
    Contents:
    egfp_functional_data.csv - Contains all functional data for the mutants
    egfp_relax_model_scores.csv - Contains all of the energy scores for all protein mutant models
    egfp_structural_data.csv - Contains all of the structure (WCN, RSA, SS) information for all mutants

plotting_scripts/
    This directory contains a script (pca_plots.R) to plot the PCA of the data
    
r_scripts/
    relax_regression_analysis.R - This a script that runs the logistic regression analysis analysis

rsa/
    Contains scripts (mac_calc_rsa.py, mac_calc_sa.py) that calculate the Solvent Accessibility (SA) and Relative Solvent Accessibility (RSA) for the pdb 4EUL. It also contains the raw calculated RSA and SA values.

scripts/
    This directory contains python scripts used in the analysis
    renumber_pdb.py - A script that renumbers pdb files
    summarize_scores.py - A script that calculates the mean score for each of the designed models

wcn/
    Contains calc_wcn.py, the script used to calculate the Weighted Contact Numbers (WCN) for each residue and 4EUL_A_wcn.csv, a file containing the weighted contact number values (WCN)
