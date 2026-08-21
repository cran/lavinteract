# lavinteract NEWS

## Version 0.1.0

* Initial release.
* `lav_slopes()`: computes conditional (simple) slopes from a fitted 'lavaan' model.
  
## Version 0.1.1

* CRAN-requested correction: quote software names in Title and Description ('lavaan').

## Version 0.2.1

* `lav_vif()`: variance inflation factors for structural predictors.
* Examples for `lav_slopes()` and `lav_vif()` now fully runnable and fast.
* Package help cleaned.
* Added helper `%||%`.

## Version 0.2.2

* CRAN compliance: documented S3 method args, fixed Rd usage, and made R code ASCII-only.

## Version 0.3.2

* `lav_cv()`: repeated holdout (Monte Carlo) cross-validation of R^2 for 'lavaan' models.

## Version 0.3.3

* `lav_cv()`: CRAN-requested fixes for ASCII-only code, explicit external function calls, and S3 method documentation.

## Version 0.3.4

* `lav_cv()`: further CRAN-requested fixes for ASCII-only code and explicit external function calls.

## Version 0.3.5

* Fixed CITATION document.  

## Version 0.4.5

* `lav_fdr()`: false discovery rate (Benjamini-Yekutieli by default) correction for selected 'lavaan' parameter p-values. 

## Version 0.4.6

* `lav_fdr()`: fixed table appearance for long predictor or outcome names.   

## Version 0.5.0

* `lav_vif()`: added support for higher-order latent predictors, fixed multigroup handling, improved automatic data recovery.
* `lav_slopes()`: improved automatic data recovery, added automatic probe values for latent moderators in single-group models, added automatic plotting ranges.
* `lav_jn()`: Johnson-Neyman regions of significance for continuous moderators.
* `lav_deltaR2()`: incremental effect sizes (part R^2 and Cohen's f^2) for structural predictors via reduced-model comparisons.
* `lav_localfit()`: residual-based local fit diagnostics.

## Version 0.5.1

* `lav_localfit()`: fixed import of `utils::head()`.

## Version 0.5.2

* `lav_cv()`: factor-score method is now selectable via `scores`, and Bartlett scoring is the new default (previously regression).
* `lav_cv()`: added `aggregation` to choose the median (default) or mean split-wise summary.
* `lav_cv()`: the proportional overfitting index (percent drop from in-sample to cross-validated R^2) is now the primary reported metric, shown next to the absolute index (in-sample minus cross-validated R^2).
* `lav_cv()`: splits are retained only if the training-set solution converges and passes the post-fit admissibility check.
* `lav_cv()`: added an internal cap on split attempts, with a warning when the requested number of usable splits cannot be reached; the number of attempted splits is now returned.

## Version 0.5.3

* `lav_cv()`: wrapped the example in `\donttest{}` to keep check time within CRAN limits.
