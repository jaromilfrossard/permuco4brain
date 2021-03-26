
<!-- README.md is generated from README.Rmd. Please edit that file -->

# permuco4brain

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/jaromilfrossard/permuco4brain/branch/master/graph/badge.svg)](https://codecov.io/gh/jaromilfrossard/permuco4brain?branch=master)
[![R-CMD-check](https://github.com/jaromilfrossard/permuco4brain/workflows/R-CMD-check/badge.svg)](https://github.com/jaromilfrossard/permuco4brain/actions)
<!-- badges: end -->

`permuco4brain` provides functions to compute permutation test in brain
imagery data. It is designed for M-EEG/ERP data. `permuco4brain` is an
add-on to `permuco` to computing cluster-mass tests, the “threshold-free
clusters-enhancement” and the Troendle’s procedure for tests
distribution in space and time (e.g.: full-scalp EEG data).

## Installation

`permuco4brain` currently works jointly with `permuco` 1.1.1 (github).
Download it using:

``` r
devtools::install_github("jaromilfrossard/permuco", build_vignettes = TRUE)
```

`permuco4brain` is still under development and you can expect changes.
Make sure to install `permuco4brain` with its full documentation:

``` r
devtools::install_github("jaromilfrossard/permuco4brain", build_vignettes = TRUE)
```

## Functionality

The main function of `permuco4brain` is `brainperm()` and it needs:

1.  `formula`: a formula object defining the design (right side) and a
    3D array containing the signals (left side). The signal is a 3D
    array visible in the global environment. Its dimension is the design
    X samples X channels. Each “rows” of the 3D array should be related
    to the corresponding row of the design (`data` argument).
2.  `data`: a dataframe containing the variables of the design.
3.  `graph`: an `igraph` object defining spatial adjacency of the
    channels (third dimension of the signal). The names of the vertices
    should correspond to the names of the 3rd dimension of signal.

By default, the `brainperm()` function produces the cluster-mass test,
but you can choose the TFCE with the argument `multcomp = "tfce"` or the
Troendle’s procedure `multcomp = "troendle"`.

`permuco4brain` uses the
[`future`](https://CRAN.R-project.org/package=future) package to handle
multi-cores computing.

You can inspect result using the `summary()` method and visualize them
using `image()`.

## Documentation

Visit <https://jaromilfrossard.github.io/permuco4brain> or check the
vignette:

1.  The
    [permuco4brain-with-eeguana](https://jaromilfrossard.github.io/permuco4brain/articles/permuco4brain-with-eeguana.html)
    vignette presents how to use `permuco4brain` in combination with
    `eeguana` :

<!-- end list -->

``` r
vignette("permuco4brain-with-eeguana", package = "permuco4brain")
```

You will find information how to extract the 3D array, the design
data-frame and the graph from an `eeg_lst` object of the `eeguana`
package. Check the `eeguana` package for pre-processing EEG data within
`R` (<https://github.com/bnicenboim/eeguana>).

2.  The
    [download-example-cheval](https://jaromilfrossard.github.io/permuco4brain/articles/download-example-cheval.html)
    vignette presents how to use `permuco4brain` with signals stored in
    `edf` files:

<!-- end list -->

``` r
vignette("download-example-cheval", package = "permuco4brain")
```

In this tutorial, you learn how to download EEG data from `zenodo` and
analyze them using `permuco4brain`.

3.  The
    [tfce](https://jaromilfrossard.github.io/permuco4brain/articles/tfce.html)
    vignette presents how to use `permuco4brain` and `future` to run the
    TFCE:

<!-- end list -->

``` r
vignette("tfce", package = "permuco4brain")
```

4.  The
    [figure-ggplot2](https://jaromilfrossard.github.io/permuco4brain/articles/figure-ggplot2.html)
    vignette presents how to use `ggplot2` to produce graphical
    representation of the results for publication:

<!-- end list -->

``` r
vignette("figure-ggplot2", package = "permuco4brain")
```

In this tutorial, you find 3 different types of plot that present the
results of a cluster-mass test or TFCE. All of them use the `ggplot2`
package and may easily be customized for publication.
