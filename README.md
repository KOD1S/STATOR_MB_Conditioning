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

2. Once the pipeline has finished running, make sure singularity is loaded (example command: module load singularity).

3. Choose a folder in which to run from, and place the trainingData.csv & MCMCgraph.csv in this folder.

4. Pull the container docker://ajnsm/py_nf_container_new

5. Run the command: singularity shell --bind /XXXX:/mnt py_nf_container_new.sif
   where XXXX is the path to the folder you are running from.

6. Copy and edit the file (in this repo) 'calcHOIsWithinMB_InteractionsOnly.py' with filepaths of trainingData & MCMCgraph, and the gene names you wish to condition on 1.
  
7. Run 'calcHOIsWithinMB_InteractionsOnly.py' in the containerised environment shell you initialised in step 5.
   This will generate .npy files containing all estimated interactions.

9. Copy and edit the file (in this repo) 'createHOIsummaries_DTuplesOnly.py' with filepaths of the 3, 4 and 5 point interaction .npy files and an output file identifier.

10. Run 'createHOIsummaries_DTuplesOnly.py' in the containerised environment shell, which will generate all dtuples. This script will not output any plots.
