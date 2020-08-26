# OscoNet
Bootstrap-based OscoNet method: 
Method to infer sinusoidal oscillations in single cell data.

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



