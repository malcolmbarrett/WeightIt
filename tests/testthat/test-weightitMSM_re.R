# lme4-style random effects terms in the formulas supplied to weightitMSM(), e.g.
# `A_1 ~ X1_0 + (1 | clus)`. These are for clustering within a time point; the
# longitudinal structure is still handled by fitting one model per time point and
# multiplying the resulting weights, exactly as without random effects.
#
# The parsing helpers (.find_re_bars()/.no_re_bars()) and the single-time-point
# behavior of each treatment type are tested in test-method_glm_re.R; what is tested
# here is that weightitMSM() routes them per time point and that they compose with
# `by`, `s.weights`, censoring, and stabilization.

skip_on_cran()
skip_if_not_installed("lme4")

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

data("msmdata", package = "WeightIt")

# msmdata ships without a grouping variable, so one is built here (in-test only,
# as elsewhere in the MSM tests) with a real random-intercept signal: the
# treatments are redrawn so cluster membership shifts the treatment probability.
# Without that the variance component is estimated at the boundary and the
# multilevel fits are indistinguishable from the fixed-effects ones.
d <- local({
  set.seed(123)

  d <- msmdata
  nclus <- 40L
  d$clus <- factor(sample(seq_len(nclus), nrow(d), replace = TRUE))

  u <- rnorm(nclus, 0, 1.5)

  for (a in c("A_1", "A_2", "A_3")) {
    d[[a]] <- rbinom(nrow(d), 1L, plogis(-0.5 + 0.2 * d$X1_0 + u[as.integer(d$clus)]))
  }

  d$G <- factor(sample(c("a", "b"), nrow(d), replace = TRUE))
  d$SW <- runif(nrow(d), 0.5, 1.5)
  d$C2 <- rbinom(nrow(d), 1L, 0.1)

  d
})

msm_formulas <- list(
  A_1 ~ X1_0 + X2_0,
  A_2 ~ X1_1 + X2_1 + A_1,
  A_3 ~ X1_2 + X2_2 + A_2
)

msm_formulas_re <- list(
  A_1 ~ X1_0 + X2_0 + (1 | clus),
  A_2 ~ X1_1 + X2_1 + A_1 + (1 | clus),
  A_3 ~ X1_2 + X2_2 + A_2 + (1 | clus)
)

test_that("A multilevel model is fit at each time point", {
  expect_no_condition({
    W <- weightitMSM(msm_formulas_re, data = d, method = "glm",
                     include.obj = TRUE)
  })

  expect_s3_class(W, "weightitMSM")
  expect_length(W$obj, 3L)

  for (i in seq_along(W$obj)) {
    expect_s4_class(W$obj[[i]], "glmerMod")
    #The random intercept is actually estimated, not pinned at 0
    expect_gt(as.numeric(lme4::VarCorr(W$obj[[i]])[["clus"]]), 0.1)
  }

  expect_true(all(is.finite(W$weights) & W$weights > 0))

  #The grouping factor is not treated as a covariate to be balanced
  expect_false("clus" %in% unlist(lapply(W$covs.list, names)))

  #Different from the fixed-effects-only fits
  expect_not_equal(W$weights,
                   weightitMSM(msm_formulas, data = d, method = "glm")$weights)

  #Mixed models supply no M-estimation parts, so the object carries none
  expect_null(attr(W, "Mparts"))
  expect_null(attr(W, "Mparts.list"))
})

test_that("The product over time points equals the per-time-point weightit() fits", {
  W <- weightitMSM(msm_formulas_re, data = d, method = "glm")

  w_each <- lapply(msm_formulas_re, function(f) {
    weightit(f, data = d, method = "glm")$weights
  })

  expect_equal(W$weights, Reduce("*", w_each), tolerance = eps)
})

test_that("Random effects can be used at a subset of the time points", {
  expect_no_condition({
    W <- weightitMSM(list(A_1 ~ X1_0 + X2_0 + (1 | clus),
                          A_2 ~ X1_1 + X2_1 + A_1,
                          A_3 ~ X1_2 + X2_2 + A_2 + (1 | clus)),
                     data = d, method = "glm", include.obj = TRUE)
  })

  expect_s4_class(W$obj[[1L]], "glmerMod")
  expect_s3_class(W$obj[[2L]], "glm")
  expect_s4_class(W$obj[[3L]], "glmerMod")
})

test_that("A formula whose only terms are random effects is not an empty formula", {
  expect_no_condition({
    W <- weightitMSM(list(A_1 ~ (1 | clus), A_2 ~ A_1 + (1 | clus)),
                     data = d, method = "glm", include.obj = TRUE)
  })

  expect_s4_class(W$obj[[1L]], "glmerMod")
  expect_true(all(is.finite(W$weights) & W$weights > 0))

  #No fixed-effect covariates at the first time point
  expect_length(W$covs.list[[1L]], 0L)
})

test_that("Random effects compose with by, s.weights, and censoring", {
  #by: a separate multilevel model within each level of `G`
  expect_no_condition({
    W_by <- weightitMSM(msm_formulas_re, data = d, method = "glm", by = ~G,
                        include.obj = TRUE)
  })

  expect_true(all(is.finite(W_by$weights) & W_by$weights > 0))
  expect_named(W_by$by, "G")

  for (i in seq_along(W_by$obj)) {
    expect_length(W_by$obj[[i]], nlevels(d$G))
    expect_s4_class(W_by$obj[[i]][[1L]], "glmerMod")
  }

  #Matches fitting each `by` group separately
  W_a <- weightitMSM(msm_formulas_re, data = d[d$G == "a", ], method = "glm")
  expect_equal(W_by$weights[d$G == "a"], W_a$weights, tolerance = eps)

  #s.weights: passed through to glmer()'s `weights`
  expect_no_condition({
    W_sw <- weightitMSM(msm_formulas_re, data = d, method = "glm",
                        s.weights = "SW", include.obj = TRUE)
  })

  expect_s4_class(W_sw$obj[[1L]], "glmerMod")
  expect_not_equal(W_sw$weights,
                   weightitMSM(msm_formulas_re, data = d, method = "glm")$weights)

  #Censoring: a censoring model may itself be multilevel
  expect_no_condition({
    W_c <- weightitMSM(list(A_1 ~ X1_0 + X2_0 + (1 | clus),
                            .cens(C2) ~ X1_1 + A_1 + (1 | clus),
                            A_2 ~ X1_1 + X2_1 + A_1 + (1 | clus)),
                       data = d, method = "glm", include.obj = TRUE)
  })

  expect_named(W_c$obj, c("A_1", "C2", "A_2"))
  expect_s4_class(W_c$obj[["C2"]], "glmerMod")

  #Censored units get a weight of exactly 0, everyone else a positive weight
  expect_identical(unname(which(W_c$weights == 0)), which(d$C2 == 1L))
  expect_true(all(is.finite(W_c$weights)))
})

test_that("Random effects compose with stabilization", {
  W0 <- weightitMSM(msm_formulas_re, data = d, method = "glm")

  expect_no_condition({
    W <- weightitMSM(msm_formulas_re, data = d, method = "glm", stabilize = TRUE)
  })

  #The denominators keep their random effects; the numerators are the saturated
  #prior-treatment models, which have none
  expect_identical(vapply(W$stabilization, deparse1, character(1L)),
                   c("~1", "~A_1", "~A_1 + A_2 + A_1:A_2"))

  num_w <- Reduce("*", lapply(seq_along(msm_formulas), function(i) {
    treat_i <- d[[all.vars(rlang::f_lhs(msm_formulas[[i]]))]]
    p_i <- fitted(glm(update(msm_formulas[[i]], W$stabilization[[i]]),
                      data = d, family = binomial()))
    get_w_from_ps(p_i, treat_i, estimand = "ATE")
  }))

  expect_equal(W$weights, W0$weights / num_w, tolerance = eps)

  #Random effects on both sides of the ratio
  expect_no_condition({
    W_re <- weightitMSM(msm_formulas_re, data = d, method = "glm",
                        stabilize = ~ (1 | clus))
  })

  expect_true(all(is.finite(W_re$weights) & W_re$weights > 0))
  expect_not_equal(W_re$weights, W$weights)
})

test_that("A mixed-model numerator suppresses M-estimation", {
  #Fixed-effects denominators, which do supply M-estimation parts, and a
  #mixed-model numerator, which does not. Keeping the denominator parts alone would
  #leave an Mparts stack that recomputes the *unstabilized* weights.
  W <- weightitMSM(msm_formulas, data = d, method = "glm",
                   stabilize = ~ (1 | clus))

  expect_null(attr(W, "Mparts"))
  expect_null(attr(W, "Mparts.list"))

  #The fixed-effects numerator at the same specification does keep them
  W_fe <- weightitMSM(msm_formulas, data = d, method = "glm", stabilize = ~X2_0)

  expect_length(attr(W_fe, "Mparts.list"), 2L * length(msm_formulas))

  skip_if_not_installed("rootSolve")
  expect_M_parts_okay(W_fe, tolerance = eps)
})

test_that("Random effects are rejected by methods that do not support them", {
  #`re_ok = FALSE` for every method but "glm" and "bart"
  expect_error(weightitMSM(msm_formulas_re, data = d, method = "cbps"),
               "[Rr]andom effects terms in `formula`")

  expect_error(weightitMSM(msm_formulas_re, data = d, method = "ipt"),
               "[Rr]andom effects terms in `formula`")

  #The check fires for the time point that has them, not only the first
  expect_error(weightitMSM(list(A_1 ~ X1_0 + X2_0,
                                A_2 ~ X1_1 + X2_1 + A_1 + (1 | clus)),
                           data = d, method = "cbps"),
               "[Rr]andom effects terms in `formula`")
})

test_that("Errors in the random effects specification are caught", {
  #Grouping variable not in `data`
  expect_error(weightitMSM(list(A_1 ~ X1_0 + (1 | nope)), data = d,
                           method = "glm"),
               "nope")

  #Missing values in the grouping variable
  d_na <- d
  is.na(d_na$clus[1:5]) <- TRUE

  expect_error(weightitMSM(list(A_1 ~ X1_0 + (1 | clus)), data = d_na,
                           method = "glm"),
               "[Mm]issing values are not allowed")
})

test_that("Objects with random effects round-trip through update() and print()", {
  W <- weightitMSM(msm_formulas_re, data = d, method = "glm")

  #`formula.list` is stored as supplied, bar terms included
  expect_identical(vapply(W$formula.list, deparse1, character(1L)),
                   vapply(msm_formulas_re, deparse1, character(1L)))

  W2 <- update(W, data = d)
  expect_equal(W2$weights, W$weights, tolerance = eps)

  out <- capture.output(print(W))
  expect_true(any(grepl("+ time 1 (A_1): X1_0, X2_0", out, fixed = TRUE)))
  expect_false(any(grepl("clus", out, fixed = TRUE)))
})
