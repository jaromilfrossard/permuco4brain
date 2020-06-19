
<!-- README.md is generated from README.Rmd. Please edit that file -->

# permuco4brain

<!-- badges: start -->

<!-- badges: end -->

`permuco4brain` provides functions to compute permutation test in brain
imagery data. It is specially designed for M-EEG/ERP data.
`permuco4brain` is an add-on to `permuco` to computing clustermass tests
on spatio-temporal data (eg: full-scalp EEG data).

## Installation

`permuco4brain` currently works with `permuco` 1.1.1 (github). Download
it using:

``` r
devtools::install_github("jaromilfrossard/permuco", build_vignettes = TRUE)
```

The package is still under development and you can expect changes.

Make sure to install `permuco4brain` with its full documentation:

``` r
devtools::install_github("jaromilfrossard/permuco4brain", build_vignettes = TRUE)
```

Functionality

The main function is `brainperm()` and it needs:

1.  `formula`: a formula object defining the design (right side) and a
    3D array containing the signals (left side). The signal is a 3D
    array visible in the global environment. Its dimension is the design
    X samples X channels. Each “rows” of the 3D array should be related
    to the corresponding row of the design (`data` argument).
2.  `data`: a data-frame containing the variables of the design.
3.  `graph`: an `igraph` object defining spatial adjacency of the
    channels (third dimension of the signal). The names of the vertices
    should correspond to the names of the 3rd dimension of signal.

You can inspect result using the `summary()` method and visualize them
using `image()`.

## Documentation

1.  The following vignette presents how to use `permuco4brain` in
    combination with `eeguana` :

<!-- end list -->

``` r
vignette("permuco4brain-with-eeguana", package = "permuco4brain")
```

You will find information how to extract the 3D array, the design
data-frame and the graph from an `eeg_lst` object of the `eeguana`
package. Check the `eeguana` package for pre-processing EEG data within
`R` (<https://github.com/bnicenboim/eeguana>).

2.  The next vignette presents how to use `permuco4brain` with signals
    stored in edf files:

<!-- end list -->

``` r
vignette("download-example-cheval", package = "permuco4brain")
```

In this tutorial, you learn how to download EEG data from zenodo and
analyze them using `permuco4brain`.

3.  The next vignette presents how to use `ggplot2` to produce graphical
    representation for publication:

<!-- end list -->

``` r
vignette("figure-ggplot2", package = "permuco4brain")
```

In this tutorial, you find 3 different types of plot that present the
results of a clustermass test. All of them use the `ggplot2` package and
may easily be customized for publication.
