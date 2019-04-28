# OptM

Estimating the optimal number of migration edges from Treemix

## Description
This package uses results from the population software 'Treemix' by [Pickrell and Pritchard (2012) DOI:10.1371/journal.pgen.1002967](https://doi.org/10.1371/journal.pgen.1002967) to estimate the optimal number of migrations edges to add to the tree. Previously, it was customary to stop adding migration edges when 99.8\% of variation in the data was explained, but optM automates this process using an _ad hoc_ statistic based on the second order rate of change in the log likelihood.  OptM has added functionality for various threshold modeling to compare with the ad hoc statistic.  The various methods are:

- "Evanno" - calculates an _ad hoc_ statistic we call deltaM based on the Evanno method, or second-order rate of change in likelihood weighted by the standard deviation.
- "linear" - estimates of the optimal M based on a piecewise linear (change point), bent cable (alpha), simple exponential (threshold, default 5\%), or non-linear least squares (threshold, default 5\%) models
- "SiZer" - a method to map and analyze derivatives for change point estimation for ecological modeling.


## Install OptM (from an R console)
- To install from CRAN
  * First install the R package 'SiZer' from CRAN using the command `install.packages("SiZer")`
  * Then install the OptM package using `install.packages("OptM")`
  * Load the package into your working R environment using `library(OptM)`

## Preparing the input files
To run OptM, you will need a folder of output files produced by Treemix v1.13.  The function optM will automatically search the folder for the _stem.llik_, _stem.modelcov.gz_, and _stem.cov.gz_ files; where "_stem_" is that provided to the _-o_ parameter of _treemix_.  It is recommended, but not required, to use _stem_ in the format _stem_\._i_\._M_; where

- _stem_ is any name you prefer
- _i_ is the iteration number for that value of _M_
- _M_ is the number of migration edges used for the treemix run (_-m_ parameter)

In order for optM to function properly, you must run:

- At least two iterations at each value of _M_ (the number of migration edges)
- _M_>2.
- The range for _M_ must be sequential integers (e.g., 1, 2, 3, etc)
- You do not need to run _M_=0 because _treemix_ automatically includes this as the null model in each run.

**NOTE:  There will be an error check to see if there is variation across iterations for each _M_.  In other words, if the data are very robust, you may get the same likelihood across all runs, thus the standard deviation across runs is zero and the _ad hoc_ statistic is undefined.  In this case, try making larger variations in the dataset (subsetting the SNPs, varying _-k_ in _treemix_, or other method of permutation/bootstrap).**

### Below is an example run of _treemix_ from a UNIX terminal for _M_={1-10} and 5 iterations per _M_:
```bash
for m in {1..10}
   do
   for i in {1..5}
      do
      treemix \
         -i test.treemix.gz \
         -o test.${i}.${m} \
         -global \
         -m ${m} \
         -k 1000
      done 
done
```

## To run OptM in R:
- First load the provided example data for a simulated dataset with 3 migration edges; and 10 iterations for _M_={1-10}
  * `folder <- system.file("extdata", package = "OptM")`
- Next, run _optM_ using the default "Evanno"-like method:
  * `test.optM = optM(folder)`
- Finally, plot the results:
  * `plot_optM(test.optM, method = "Evanno")`
  
  
- Alternatively, run using various linear modeling estimates rather than the _ad hoc_ statistic:
  * `folder <- system.file("extdata", package = "OptM")`
  * `test.linear = optM(folder, method = "linear")`
  * `plot_optM(test.linear, method = "linear")`

- OR using _SiZer_:
  * `folder <- system.file("extdata", package = "OptM")`
  * `test.sizer = optM(folder, method = "SiZer")`
  * `plot_optM(test.sizer, method = "SiZer")`
  

# Version History
- Version 0.1.3, 2019/4/23
  * The read.treemix function now searches for all treemix input files, and the specially formatted _stem_ is no longer required.  Thanks Jie Zhong!!!
- Version 0.1.2, 2019/3/1
  * Fixed typo in DESCRIPTION - Pickrell and Pritchard 2012, not 2002
  * In the `plot_optM`, changed the plotting color to have an alpha (semi-transparent) fill and Y-axis labels for Δm
- Version 0.1.1, 2019/1/2
  * Released the first version


## Citation
Fitak, R. R. (2018) optM: an R package to optimize the number of migration edges using threshold models. Journal of Heredity [in prep]

- Or enter the command `citation("OptM")` into your R console

## Contact
Robert Fitak  
Department of Biology  
Duke University  
USA  
rfitak9@gmail.com  