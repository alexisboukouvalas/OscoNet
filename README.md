# OscoNet
OscoNet code.

## A simple example

Just run OscoNetHypothesis.py in ipython.

Assuming all packages are installed, this will demonstrate the hypothesis test on a simple
toy example with 100 genes and 100 cells.

The code should take approximately 1 minute to run.
This is the output you should get:

```
G=10, N=100
Removing genes Index([], dtype='object')
postprocessed data size (10, 100)
DoTheBootstrap:: Running on 1  cores with 10 target genes.
DoTheBootstrapFast: Initial search Using NP
DoTheBootstrapFast: multiple passes per gene
Bootstrap parallel time with 1 workers took 0.5 secs. Bootstrap samples 100. Use TF 0.
Completed in 0.545126 seconds
qvalues pi0=1.091, estimated proportion of null features
got pi0 > 1 (1.091) while estimating qvalues, setting it to 1
Is adjacency matrix symmetric? True
Adjacency matrix now symmetric after applying AND transformation
Sparseness 0.266667
Writing edge network to toy.csv number of edges 12
False positive rate 0.0 True positive rate 1.0 False discrovery rate 0.0
True positives= 24 True negatives= 76 False positives= 0 False negatives 0
```

The program also generate a csv file ```toy.csv``` containing an edge network representation of
adjacency graph. The latter contains edges only on significant pairs
of co-oscillating genes.


To get a bigger network just change line
```
data, phaseG, angularSpeed = GetSimISyntheticData(fPlot=False, NG=10, G=100, N=100, noiseLevel=0, ngroups=2, fReturnTruth=True)
``` in ```OscoNetHypothesis.py```

```G``` is the number of genes, ```NG``` the number of genes in a group, ```ngroups``` the number of groups, ```N``` the number of cells.

The other parameter in the code is ```NumPerm=100```
which is the number of bootstrap samples.


## Dependencies

You need to install joblib
```pip install joblib```

