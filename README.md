# `permuco4brain`


`permuco4brain` currently works with permuco 1.1.1 (github) . Download it using: 

`devtools::install_github("jaromilfrossard/permuco")`

This package is an add-on to permuco for computing clustermass tests on spatio-temporal data (eg: full-scalp EEG data).

The package is still under developpement and you can expect changes.

# Installation

Make sure to install `permuco4brain` with its full documentation:

`devtools::install_github("jaromilfrossard/permuco4brain", build_vignettes = TRUE)`


# Functionality

The main function is `brainperm()` and it needs:

1. `formula`: a formula object definig the design (right side) and a 3D array containing the signals. The signal is a 3D array stored in the global environment. Its dimension is the design X times X channel. Each "rows" of the 3D array should be related to the corresponding row of the data.
2. `data`: a dataframe containing the variables of the design.
3. `graph`: an igraph object defining spatial adjacency of the channels (third dimension of the signal). The names of the vertices should correspond to the names of the 3rd dimension of signal.

You can inspect result using the `summary()` method and vizualize them using `image()`.


# Documentation

1. The following vignette presents how to use `permuco4brain` in combinaison with `eeguana` :

`vignette("permuco4brain-with-eeguana", package = "permuco4brain")`


You will find information how to extract the 3D array, the design dataframe and the graph from an `eeg_lst` object of the `eeguana` package. Check the `eeguana` package for preprocessing EEG data within `R` (https://github.com/bnicenboim/eeguana).


2. The next vignette presents how to use `permuco4brain` with signals stored in edf files:

`vignette("download-example-cheval", package = "permuco4brain")`

By this tutorial you download EEG data from zenodo and analyze them.


3. The next vignette presents how to use `ggplot2` to produce graphical representation for publication:

`vignette("figure-ggplot2", package = "permuco4brain")`


In this tutorial, you will find 3 differents type of plot that can be used to present the results of a clustermass test. All of them use the `ggplot2` package and may easly be custom for publication.
