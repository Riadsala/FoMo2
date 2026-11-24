##########################################################################
# Haggouch, Doré-Mazars & Tagu (2026): working with a fitted model.
# Reads the saved fit in fits/ (made by fit_haggouch2026acta.R; no refitting),
# then looks at the data, convergence, the posterior, accuracy and
# posterior predictive checks.
# haggouch2026acta.qmd renders this script; the "## ----" lines mark
# its chunks. Slow chunks are cached by the .qmd.
#
# The largest dataset here (340,480 selections), so the predictions and
# simulations use fewer posterior draws than the other pages: enough for
# accuracy and trial statistics. The predictions are about as big as
# clarke2022qjep's (50 draws x 340,480 selections against 200 x 87,760).
##########################################################################

## ---- load
source("../../../functions/set_up_FoMo.R")

# the fit this page reads: fit_haggouch2026acta.R (in this folder) makes it
fit_files <- "fits/haggouch2026acta_fourier124.rds"
if (!all(file.exists(fit_files))) stop("Run fit_haggouch2026acta.R first: it makes the fit this page reads")

n_post_draws <- 2    # posterior draws for the posterior predictive check
n_pred_draws <- 50   # posterior draws for per-selection predictions (accuracy)

## ---- data
d <- import_data("haggouch2026acta")

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
fit <- readRDS("fits/haggouch2026acta_fourier124.rds")
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
# 304 people: one row each, too many to label
no_person_labels <- theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank())
plot_post_per_person_spatial(post) + no_person_labels

## ---- fig-person-cat
plot_post_per_person_cat(post) + no_person_labels

## ---- predict
# the fit doesn't store per-selection predictions (save_selections =
# FALSE), so predict them for a subset of draws
pred <- predict_fomo(fit, sl, n_draws = n_pred_draws)
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

