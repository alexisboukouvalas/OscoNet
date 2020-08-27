# OscoNet
Bootstrap-based OscoNet method: 
Method to infer sinusoidal oscillations in single cell data. This software reproduces the approach presented in
'OscoNet: inferring oscillatory gene networks'
by Luisa Cutillo, Alexis Boukouvalas, Elli Marinopoulou, Nancy Papalopulu & Magnus Rattray 
BMC Bioinformatics volume 21, Article number: 351 (2020)
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03561-y
contact: l.cutillo@leeds.ac.uk

## Installation
1. Create new environment `conda create --name fullosconet python=3` and to activate it with 
`conda activate fullosconet` before proceeding.
1. Install required packages using `pip install -r requirements.txt`
1. Install package `pip install -e .`
1. Install numba `conda install numba`
1. Verify your installation by running `pytest` from the project root directory `FullOscoNet`.
Note this can take around 1-2 minutes.


## Synthetic data
1. Run ` python OscopeBootstrap/oscope_tf.py --test_mode
` for a simple demonstration of the method on synthetic data. This will run under a quick configuration
to demonstrate the capabilities of the method. This should take 10-20 seconds.
1. Remove the `--test_mode` flag for a 1000-sample bootstrap test on the exact synthetic run configuraiton used in the paper
(1000 genes, 100 cells with 3 clusters of co-oscillating genes).

## Notebooks
1. `notebooks/OscoNet introduction.ipynb`: provides an introduction to the hypothesis test on a simple synthetic example.
2. `notebooks/Reproduce_figures_5_7.ipynb` : pseudotime on Whitfield microarray data. To see how the spectral embedding
pseudotime method can be applied.
3. `notebooks/Reproduce_Table1.ipynb`: Reproduce table 1 from OscoNet paper
4. `notebooks/Reproduce_Table5.ipynb`: Reproduce table 5 from OscoNet paper
5. `notebooks/Reproduce_Figure4.ipynb`: Reproduce Figure 4 from OscoNet paper



