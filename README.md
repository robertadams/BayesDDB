# BayesDDB

Bayesian Dynamic Borrowing with Double Adjustment for augmenting clinical trial data with external data.

## Installation

```r
# From CRAN (when available)
install.packages("BayesDDB")

# Development version (public)
pak::pak("robertadams/BayesDDB")

# Development version (internal)
remotes::install_github("bayer-int/BayesDDB", auth_token =  <PAT>)
```

## Usage

The package provides four main analysis functions and a utility:

- `bdb_single_arm_bin()` — single-arm trial, binary outcome
- `bdb_single_arm_cont()` — single-arm trial, continuous outcome (requires JAGS)
- `bdb_two_arms_bin()` — two-arm trial, binary outcome
- `bdb_two_arms_cont()` — two-arm trial, continuous outcome (requires JAGS)
- `elastic_function()` — arctangent elastic function for outcome discounting

Minimal example with the elastic function:

```r
library(BayesDDB)
elastic_function(c(0.1, 0.5, 0.9), a = 1)
```

For full methodology and runnable examples, see the vignette:

```r
vignette("introduction", package = "BayesDDB")
```

To install from source with the vignette built: `devtools::install(build_vignettes = TRUE)`.

## License

GPL-3. See [LICENSE.md](LICENSE) for details.
