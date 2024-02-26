# STATOR_MB_Conditioning
Scripts which utilize the STATOR pipeline to allow runs to use pre-existing adjacency matrices for testing different Markov blanket conditions.

In order to estimate interactions between a set of genes from sc-RNAseq data, STATOR (https://doi.org/10.1101/2023.12.18.572232) must condition on all other genes being zero.
In practice this would leave very few or no samples with which to estimate from, so only genes within the Markov blanket (Pearl, Judea 1988 ISBN 0-934613-73-7) are conditioned upon.
By default, STATOR conditions on these genes being zero, but some interactions may only be estimable under certain Markov blanket conditions. To avoid running the entire pipeline
for each condition, these additional files have been created. These scripts utilize pre-existing binarized & filtered scRNA-seq data, and precalculated adjacency matrix. For this
sub-pipeline give valid results, the same dataset must be used as was used for the initial STATOR run which generated these files.
Note that these scripts require singularity installed.

1. Run the STATOR pipeline (available at https://github.com/AJnsm/Stator)
    This will output a trainingData.csv file and MCMCgraph.csv file, available inside the 'output' folder.



2. Copy and edit the file (in this repo) 'calcHOIsWithinMB_InteractionsOnly.py' with filepaths of trainingData & MCMCgraph, and the gene names you wish to condition on 1.
  
3. Run this file in the environment used to run STATOR (load singularity

5. In the same environment used to run STATOR (with singularity module loaded), run the file 
    This will output a set of .npy files containing the interactions to the folder the file was run from.

