# Longitudinal treatment with censoring at more than one time point.
#
# Three treatment time points with a censoring model before the second and before
# the third, so the risk set shrinks twice:
#
#   A_1  ->  C_2  ->  A_2  ->  C_3  ->  A_3
#
# Once a unit is censored it is gone: its later censoring indicators, treatments,
# and covariates are all NA, and its final weight is exactly 0.

library("WeightIt")

# The `plot()` calls below print a ggplot. Printing one with no device open opens the
# default device, which outside an interactive session is `pdf()` -- leaving an
# `Rplots.pdf` behind in the working directory. Send those to the null device instead;
# interactively the usual device is used and the plots show up as they should.
if (!interactive()) {
  grDevices::pdf(NULL)
}

set.seed(31)

n <- 3000

d <- data.frame(X1_0 = rnorm(n),
                X2_0 = rbinom(n, 1L, .45))

# --- time 1: treatment, nobody censored yet ---------------------------------
d$A_1 <- rbinom(n, 1L, plogis(-0.2 + 0.4 * d$X1_0 + 0.3 * d$X2_0))

d$X1_1 <- rnorm(n, 0.5 * d$A_1 + 0.4 * d$X1_0)
d$X2_1 <- rbinom(n, 1L, plogis(-0.3 + 0.5 * d$A_1))

# --- censoring before time 2 -------------------------------------------------
# 1 = censored (drops out), 0 = still under observation
d$C_2 <- rbinom(n, 1L, plogis(-2.0 + 0.5 * d$X1_1 - 0.4 * d$A_1))

out_2 <- d$C_2 == 1L

# --- time 2: treatment, among those still under observation ------------------
d$A_2 <- rbinom(n, 1L, plogis(-0.2 + 0.4 * d$X1_1 + 0.5 * d$A_1))

d$X1_2 <- rnorm(n, 0.5 * d$A_2 + 0.4 * d$X1_1)
d$X2_2 <- rbinom(n, 1L, plogis(-0.3 + 0.5 * d$A_2))

# --- censoring before time 3, among those who survived the first -------------
d$C_3 <- rbinom(n, 1L, plogis(-2.0 + 0.5 * d$X1_2 - 0.4 * d$A_2))

out_3 <- out_2 | d$C_3 == 1L

# --- time 3: treatment -------------------------------------------------------
d$A_3 <- rbinom(n, 1L, plogis(-0.2 + 0.4 * d$X1_2 + 0.5 * d$A_2))

# Blank out everything a unit could not have contributed after it dropped out.
# `C_3` included: a unit censored at time 2 has no time-3 censoring indicator.
is.na(d[out_2, c("C_3", "A_2", "X1_2", "X2_2", "A_3")]) <- TRUE
is.na(d[out_3, "A_3"]) <- TRUE

cat(sprintf("n = %s | censored at C_2: %s | additionally at C_3: %s | still in at time 3: %s\n",
            n, sum(out_2), sum(out_3) - sum(out_2), sum(!out_3)))

# --- the models, in temporal order -------------------------------------------
# Censoring models are interleaved with the treatment models by wrapping the
# indicator in .cens() on the left side.
f <- list(A_1        ~ X1_0 + X2_0,
          .cens(C_2) ~ X1_1 + X2_1 + A_1,
          A_2        ~ X1_1 + X2_1 + A_1,
          .cens(C_3) ~ X1_2 + X2_2 + A_2 + A_1,
          A_3        ~ X1_2 + X2_2 + A_2 + A_1)

W <- weightitMSM(f, data = d, method = "glm")

W

# Censoring is a treatment type, not a separate kind of model: the indicators sit
# among the treatments in `treat.list`, and their covariates in `covs.list`, in the
# order the models were fit. There is no separate `cens.list`.
names(W$treat.list)
vapply(W$treat.list, WeightIt:::get_treat_type, character(1))

# Everyone censored at any point has a final weight of exactly 0
stopifnot(all(W$weights[out_3] == 0), all(W$weights[!out_3] > 0))

# `at.risk` has one column per model, in the same order
colSums(W$at.risk)


# ---------------------------------------------------------------------------
# summary(): one entry per model, treatments and censoring alike
# ---------------------------------------------------------------------------

s <- summary(W)

names(s)      # A_1, C_2, A_2, C_3, A_3

summary(W)

# Each entry summarizes the same weights -- the product across all five models --
# over a different sample: the units that model was fit on.
#
# WORTH A LOOK. With censoring at more than one point, a censoring entry is not
# free of zero weights. `C_2` covers the 2635 units still under observation at
# that point, but 339 of them are censored later at `C_3`, and their *final*
# weight is 0. Only the last censoring model has none.
#
# The alternative would be for each entry to report the weight contributed by its
# own model rather than the final product, but those are not the weights anyone
# uses, and the treatment entries have always reported the product.
# The effective sample size a censoring entry reports is measured against the
# units that model was fit on -- those under observation entering it -- so it
# shows what that censoring cost rather than comparing to a sample the units were
# never drawn from.
data.frame(at_risk_for_model = colSums(W$at.risk),
           ess_row_reports = sapply(s, function(x) sum(x$effective.sample.size["Unweighted", ])),
           n_zero_weights = sapply(s, function(x) sum(x$num.zeros)),
           row.names = names(s))


# ---------------------------------------------------------------------------
# which.time: restrict the summary, as in cobalt::bal.tab()
# ---------------------------------------------------------------------------

# Both censoring models, by name
summary(W, which.time = c("C_2", "C_3"))

# By position in `formula.list`. Note that positions are kept: A_3 is model 5
# however few models are asked for.
summary(W, which.time = 5)

# Omit the argument for all of them
length(summary(W))


# ---------------------------------------------------------------------------
# plot(): one model at a time
# ---------------------------------------------------------------------------

# A censoring model plots only the units still under observation, since the rest
# have a weight of 0 and are not part of the weighted sample
plot(s, which.time = "C_3")

# A treatment measured after a censoring model likewise leaves out the units that
# had already dropped out
plot(s, which.time = "A_3")

# `time` is the former name of this argument and still works
plot(s, time = "C_3")


# ---------------------------------------------------------------------------
# Balance
# ---------------------------------------------------------------------------

if (requireNamespace("cobalt", quietly = TRUE)) {
  # Recent cobalt assesses each model on the sample it was fit on, censoring
  # models included
  print(cobalt::bal.tab(W))
}
