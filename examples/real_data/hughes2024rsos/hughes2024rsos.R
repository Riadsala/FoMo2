##########################################################################
# Hughes, Statham & Clarke (2024): working with a fitted model.
# Reads the saved fits in fits/ (made by fit_hughes2024rsos.R; no refitting),
# then looks at the data, convergence, the posterior, accuracy and
# posterior predictive checks.
# hughes2024rsos.qmd renders this script; the "## ----" lines mark
# its chunks. Slow chunks are cached by the .qmd.
##########################################################################

## ---- load
source("../../../functions/set_up_FoMo.R")

# the fits this page reads: fit_hughes2024rsos.R (in this folder) makes them
fit_files <- c("fits/hughes2024rsos_fourier124.rds",
                "fits/hughes2024rsos_scarcity_fourier124.rds")
if (!all(file.exists(fit_files))) stop("Run fit_hughes2024rsos.R first: it makes the fits this page reads")

n_post_draws <- 20   # posterior draws for the posterior predictive check

## ---- data
# feature and conjunction foraging, pooled over the scarcity levels
# (scarcity = FALSE), as fitted
d <- import_data("hughes2024rsos", scarcity = FALSE)

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
fit <- readRDS("fits/hughes2024rsos_fourier124.rds")
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
# the same trials with each scarcity level its own condition (K = 6: task x
# scarcity), against the main fit (K = 2: task)
d6 <- import_data("hughes2024rsos", scarcity = TRUE)
fit6 <- readRDS("fits/hughes2024rsos_scarcity_fourier124.rds")
fits <- list(`K = 2` = fit, `K = 6` = fit6)
# loo for each fit, from the per-trial log-likelihoods the fits store
loos <- map(fits, \(f) f$loo())

## ---- compare
# accuracy needs each fit's data list (sls) to predict the selections
cmp <- fomo_compare(fits, d, loos = loos, n_draws = 200,
                    sls = list(`K = 2` = sl, `K = 6` = attr(fit6, "fomo_sl")))
format_compare(cmp, columns = c("model", "deterministic", "nn", "gain_nn", "elpd", "elpd_diff",
                                "se_diff", "pct_high_k", "robust"))

## ---- fig-compare
# where the difference comes from: by task x scarcity
plot_post_elpd_by_condition(loos[c("K = 6", "K = 2")], d6$stim)

## ---- papers-scale
# The paper's parameters (see the single-level cat example): b_a, the bias
# towards class A (1), is the average of the transition matrix's two rows'
# log-odds of class A over class B; b_s, the stick/switch rate, is half
# their difference. Computed for every posterior draw of the population
# values in each task x scarcity condition
post6 <- extract_post(fit6, d6$found, d6$stim)

papers_scale <- post6$cat |>
  filter(dist == "posterior") |>
  select(.draw, condition, prev_class, class, value) |>
  pivot_wider(names_from = class, names_prefix = "p", values_from = value) |>
  mutate(logodds = log(p1 / p2)) |>
  select(.draw, condition, prev_class, logodds) |>
  pivot_wider(names_from = prev_class, names_prefix = "after_", values_from = logodds) |>
  mutate(b_a = (after_1 + after_2) / 2,
         b_s = (after_1 - after_2) / 2) |>
  separate(condition, c("task", "scarcity"), sep = "_")

## ---- fig-scarcity
# each scarce condition minus the equal one, within each task. A bias
# towards the scarce type would push b_a up when A is scarce and down when
# B is scarce
papers_scale |>
  select(.draw, task, scarcity, b_a, b_s) |>
  pivot_longer(c(b_a, b_s), names_to = "parameter") |>
  pivot_wider(names_from = scarcity) |>
  mutate(`A scarce` = scarceA - equal,
         `B scarce` = scarceB - equal) |>
  pivot_longer(c(`A scarce`, `B scarce`), names_to = "comparison", values_to = "difference") |>
  mutate(parameter = recode(parameter, b_a = "b_a: bias towards A", b_s = "b_s: stick/switch")) |>
  ggplot(aes(difference, comparison)) +
  theme_fomo() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggdist::stat_halfeye(.width = FOMO_HPDI_WIDTHS, point_interval = "median_hdci") +
  facet_grid(task ~ parameter, scales = "free_x") +
  labs(x = "difference from the equal condition (population)", y = NULL)
