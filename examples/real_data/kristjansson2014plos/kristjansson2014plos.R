##########################################################################
# Kristjansson, Johannesson & Thornton (2014): working with a fitted model.
# Reads the saved fits in fits/ (made by fit_kristjansson2014plos.R; no refitting),
# then looks at the data, convergence, the posterior, accuracy and
# posterior predictive checks.
# kristjansson2014plos.qmd renders this script; the "## ----" lines mark
# its chunks. Slow chunks are cached by the .qmd.
##########################################################################

## ---- load
source("../../../functions/set_up_FoMo.R")

# the fits this page reads: fit_kristjansson2014plos.R (in this folder) makes them
fit_files <- c("fits/kristjansson2014_fourier124.rds",
                "fits/kristjansson2014_fourier124_symmetric.rds")
if (!all(file.exists(fit_files))) stop("Run fit_kristjansson2014plos.R first: it makes the fits this page reads")

n_post_draws <- 20   # posterior draws for the posterior predictive check

## ---- data
d <- import_data("kristjansson2014plos")

d$stim |>
  summarise(people = n_distinct(person), trials = n_distinct(trial),
            items_per_trial = n() / n_distinct(trial), .by = condition) |>
  left_join(summarise(d$found, selections = n(), .by = condition), by = "condition") |>
  knitr::kable()

## ---- fig-trial
plot_trial(d$found, d$stim)

## ---- fig-runs
plot_runs(compute_runs(d$found))

## ---- trial-stats
trial_stats <- compute_trial_stats(d$found)
head(trial_stats) |> knitr::kable(digits = 2)

## ---- read-fit
# the fit, and the data list it was made from (fit_fomo() stores it with the fit)
fit <- readRDS("fits/kristjansson2014_fourier124.rds")
sl <- attr(fit, "fomo_sl")

## ---- fit-summary
fomo_summary(fit, cond_labels = levels(d$stim$condition))

## ---- fig-traces
plot_post_traces(fit)

## ---- extract
post <- extract_post(fit, d$found, d$stim)

## ---- fig-post
plot_post(post)

## ---- fig-post-weights
plot_post_spatial(post, show = "weights")

## ---- fig-person
plot_post_per_person_spatial(post)

## ---- fig-person-cat
plot_post_per_person_cat(post)

## ---- predict
# the fit doesn't store per-selection predictions (save_selections =
# FALSE), so predict them for a subset of draws
pred <- predict_fomo(fit, sl, n_draws = 200)
acc <- extract_post(pred, d$found, d$stim)$acc

## ---- fig-accuracy
plot_post_accuracy(acc, d$stim, d$found)

## ---- sim-posterior
# whole trials simulated from the fit; the model doesn't predict when a
# trial stops, so keep as many simulated selections as the data have
sim_post <- sim_fomo(sl, fit, iter = n_post_draws)
extract_pred(sim_post, d$stim) |>
  semi_join(d$found, by = c("trial", "n")) -> found_post

## ---- fig-sim-trial
plot_trial(found_post, d$stim)

## ---- fig-sim-stats
sim_stats <- compute_trial_stats(found_post)
plot_trial_stats(trial_stats, sim_stats)

## ---- compare-load
# the same data fitted with a symmetric transition matrix: one stay
# probability per condition, every other class equally likely
# (prep_stan_list(symmetric_v = TRUE))
fit_sym <- readRDS("fits/kristjansson2014_fourier124_symmetric.rds")
fits <- list(full = fit, symmetric = fit_sym)
# loo for each fit, from the per-trial log-likelihoods the fits store
loos <- map(fits, \(f) f$loo())

## ---- compare
cmp <- fomo_compare(fits, d, loos = loos)
format_compare(cmp, columns = c("model", "deterministic", "nn", "gain_nn", "elpd", "elpd_diff",
                                "se_diff", "pct_high_k", "robust"))

## ---- fig-compare
# where the difference comes from: the loo difference summed over each
# condition's trials
plot_post_elpd_by_condition(loos, d$stim)
