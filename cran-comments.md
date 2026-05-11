# CRAN submission comments

## Package

BayesDDB implements Bayesian dynamic borrowing with double adjustment for clinical trials augmented by external data.

## Dependencies

- **Imports:** dplyr, overlapping, ggplot2, ggpubr (required at runtime).
- **Suggests:** rjags — required only for continuous-outcome functions (`bdb_single_arm_cont`, `bdb_two_arms_cont`). The system dependency JAGS (https://mcmc-jags.sourceforge.io/) must be installed separately for those functions. Binary-outcome functions and `elastic_function` do not require rjags or JAGS. **Check without rjags:** Because rjags/JAGS may be unavailable on some systems (e.g. CRAN build servers), a complete check can be run with `_R_CHECK_FORCE_SUGGESTS_=false` so that missing Suggested packages do not cause an error.

## Testing

- `testthat` tests cover `elastic_function`, internal `mean_ci`, and the four main `bdb_*` functions.
- Tests that require rjags/JAGS use `skip_if_not_installed("rjags")` and `skip_on_cran()` so CRAN checks pass without JAGS.

## Notes

- No compiler or system libraries required beyond R and the listed R packages.
- Vignette builds with knitr/rmarkdown (Suggests).
