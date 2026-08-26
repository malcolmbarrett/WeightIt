## To Do
* Create function for estimating treatment effects to remove marginaleffects as dependency (possibly in a new package)
* Implement RieszBoost for GBM weighting
* Drop `"cobalt.treat"` from `.treat_classes` in `R/treat.R`, once *cobalt* registers its `[` method on `treat` rather than on `cobalt.treat`. It is on `cobalt.treat` only because *WeightIt* 2.0.0 registered a competing `[.treat`; *cobalt* cannot move it while that version is the one on CRAN. See `cobalt/_dev/cens-transition.md`.

## Larger future directions

### Architecture

* **Unify the "weight one sample to fixed targets" problem.** IPCW (`.cens()`), the ATT, survey calibration, and transportability to an external population are all the same shape: weight one group so its covariate moments match a supplied target vector. `weightit2optweight.cens()` already goes through `optweight::optweight.svy.fit()`, and `weightit2ebal.cens()` is the ATE branch restricted to one group. Factoring this into a single internal solver interface (targets + tolerances + a divergence) would let every optimization-based method support all four cases from one implementation, instead of each gaining a bespoke `.cens` branch.
* **Make `Mparts` a real class rather than a list of closures.** The M-estimation parts are passed around as a bare list and manipulated by helpers that must stay mutually consistent (`.invert_num_Mpart()`, `.expand_Mparts_by()`, `.att_out_to_cens()`, the stacking in `.compute_vcov()`). `.expand_Mparts_by()` silently ignores the `Xtreat`/`A` it is handed in favor of ones it closed over, which is invisible from the outside and easy to trip over. A class with `stack()`/`invert()`/`expand()` methods and validity checks would make the contract explicit and the failure loud.
* **Replace the name-mangling method dispatch with a registry.** Methods are found by pasting a suffix onto a function name (`weightit2glm` + `.cens`/`.multi`/`.cont`), with `.weightit_methods` as a parallel table of capabilities that has to be updated by hand in lockstep. A proper registry (or S3 class per method) would let other packages contribute methods, and would remove the need to keep the table and the function names in sync.
* **Split `weightitMSM2cbps()`.** It now handles five treatment types inline in one function, with per-type behavior spread across a covariate-prep `switch()`, a `coef_ind` `switch()`, and five parallel closure lists. It is also less featureful than the point-treatment version (no `link` support, empirical rather than closed-form `Sigma()`). Factor the closure construction into per-type constructors and close the feature gap.
* **Scale the distance-based methods past ~10k units.** `method = "energy"` and `method = "cfd"` build dense n x n matrices and hand them to a QP. Nystrom approximation or random Fourier features would make them usable on large data, where they are currently infeasible.

### New weighting methods

* **Kernel balancing (`method = "kbal"`)** as a first-class method rather than the user-defined-function example in `?method_user`.
* **Transportability / generalizability weights**, targeting an external population's covariate distribution rather than the study sample. Falls out almost free from the targets refactor above.
* **Cross-fitted / debiased weight estimation**, where nuisance models are fit out-of-fold. Relevant to the machine-learning methods (`gbm`, `bart`, `super`), whose plug-in weights have no valid asymptotic variance at present.
* **Weight uncertainty from a posterior** rather than from M-estimation, e.g. propagating BART posterior draws through the weights so effect estimates reflect the full posterior. Would complement the existing M-estimation and bootstrap machinery.
* **Survival-outcome weighting beyond IPCW**, e.g. time-varying censoring weights on a long-format dataset rather than the current one-row-per-unit wide format.

### Diagnostics and testing

* ~~**Target-based balance assessment.**~~ Done in *cobalt* 5.0.0, which compares the weighted uncensored units against the at-risk sample and handles a censoring model sitting among longitudinal treatments. `DESCRIPTION` and `test-censoring.R` now name that version rather than a placeholder.
* **A fast test tier that runs on CRAN.** 277 of 286 `test_that()` blocks are behind `skip_on_cran()`, so CRAN effectively runs 9 of them. A small fast subset that always runs would catch platform-specific breakage that currently only surfaces locally.
* **Coverage in CI**, with the `NOT_CRAN=true` and non-parallel caveats baked in (see below) so the number is meaningful.
