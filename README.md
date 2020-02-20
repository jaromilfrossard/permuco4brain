# `permuco4brain`


WARNING: permuco4brain currently works with the `flip` branch of permuco. Download it using: 

`devtools::install_github(jaromilfrossard/permuco, ref = "flip")`



This package is an add-on to permuco to compute cluster-mass test on spatio-temporal data, like full-scalp EEG data.

The package is still in developpement and will change.


The main function is brainperm() and it needs:

1. formula: a formula object definig the design (right side) and a 3D array containing the signals. Its dimension is the design X times X channel. Each "rows" of the 3D array should be related to the corresponding row of the data.
2. data: a dataframe containing the variables of the design.
3. graph: an igraph object defining spatial adjacency of the channels (third dimension of the signal). The names of the vertices should correspond to the names of the 3rd dimension of signal.




