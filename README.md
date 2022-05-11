# robcovsel

- This package provides the functions to compute the ALRP (adaptive lasso based on robust pariwise correlations) and ALGR (adaptive lasso based on Gaussrank correlations) proposed by Peng Su, Samuel Muller and Garth Tarr. The manuscript could be found here (https://arxiv.org/abs/2110.12406).

- We added a demonstration (demo), a simulation workflow (simu) and two real data applications (boston and breastcancer) in vignettes

To get started, you can install the package using:

```r
remotes::install_github("PengSU517/robcovsel")
```

You may also consider installing the `shootings` package. The functions in this package were taken from https://github.com/ineswilms/sparse-shooting-S and compiled into a basic package to make them more convenient to use:

```r
remotes::install_github("PengSU517/shootings")
```


