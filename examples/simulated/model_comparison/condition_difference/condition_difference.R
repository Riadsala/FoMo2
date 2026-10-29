##########################################################################
# Do two conditions differ?
#
# We simulate three datasets in which ten people each do two conditions:
#   "average":    on average, people stick with the item class they just
#                 picked more in condition one than in condition two
#   "individual": no difference on average, but each person behaves a
#                 little differently in the two conditions
#   "none":       each person behaves the same in both conditions
# and fit each dataset with three models:
#   "condition":          separate population values for each condition
#   "person x condition": one set of population values, but each person's
#                         values differ between the two conditions
#   "pooled":             no condition effect at all
# Then we ask (1) how big the differences are, from the condition model,
# and (2) which model predicts each dataset best (loo).
#
# Fitting the nine models takes about an hour, so this script saves them in
# condition_difference_fits.rds, and the page (condition_difference.qmd)
# loads that file instead of fitting. To refit, delete the file and run:
#   Rscript condition_difference.R
##########################################################################

## ---- load
source("../../../../functions/set_up_FoMo.R")

# sampler settings: warmup and sampling iterations per chain (the same in
# every example)
n_iter <- 500

# where the fitted models are saved (see "fit" below)
fits_file <- "condition_difference_fits.rds"

## ---- stimuli
# 10 people, each doing 10 trials in each of two conditions, on an 8 x 8
# grid with two item classes
spec <- list(K = 2, nrow = 8, ncol = 8, prop_missing = 0.25, v_type = "cat", n_classes = 2)
stim <- generate_stim(spec, n_trials = 10, n_people = 10)

## ---- relabel
# The "pooled" and "person x condition" models are the same model as the
# "condition" model, fitted to relabelled data.

# merge_conditions(): call both conditions "both", so the model estimates
# one set of values that applies to every trial
merge_conditions <- function(d) {
  mutate(d, condition = factor("both"))
}

# person_by_condition(): as above, but also give each person a new number
# in each condition, so that person 1 in condition one and person 1 in
# condition two are treated as two different people
person_by_condition <- function(d) {
  d |>
    mutate(person = paste(person, condition),
           person = as.integer(factor(person, levels = unique(person)))) |>
    merge_conditions()
}

## ---- simulate
# Stickiness: the probability of picking the same item class again.
# "sticky" stays with probability 0.8, "no_bias" with 0.5 (no preference).
# (Each column of the matrix gives the next class after one class.)
sticky  <- matrix(c(0.8, 0.2,
                    0.2, 0.8), nrow = 2)
no_bias <- matrix(0.5, nrow = 2, ncol = 2)

# Everything else is the same everywhere: attraction to nearby items,
# turning, and a preference for moving right and along the grid. People
# vary around these population values (sd_lambda sets how much: a larger
# value means less spread). n_cond: the number of conditions
other_params <- function(n_cond) {
  list(rho_delta_mu = rep(1, n_cond),
       rho_psi_mu = rep(-0.5, n_cond),
       fourier_coef_mu = matrix(c(0.8, 0, 0, 0, 0.5, 0), nrow = n_cond, ncol = 6, byrow = TRUE),
       sd_lambda = list(b_v = 3,
                        rho_delta = rep(6, n_cond),
                        rho_psi = rep(6, n_cond),
                        fourier_coef = c(8, 8, 8)))
}

# Simulate one dataset, then put the real person and condition labels back
# (the "none" data are simulated with the conditions merged)
simulate <- function(stim_for_sim, params) {
  sl <- prep_stan_list(stim_for_sim, params)
  sim <- sim_fomo(sl, iter = 1, chains = 1)
  sim$draws()   # read the simulated values into memory, so they can be saved

  found <- extract_pred(sim, stim_for_sim) |>
    truncate_trials(lambda = 2, min_keep = 10) |>
    select(-person, -condition) |>
    left_join(distinct(stim, trial, person, condition), by = "trial")

  list(found = found, sim = sim)
}

sims <- list(
  # each condition has its own population values, and each person's values
  # are drawn separately in each condition
  average    = simulate(stim, c(other_params(2), list(b_v_pr = list(sticky, no_bias)))),
  individual = simulate(stim, c(other_params(2), list(b_v_pr = list(no_bias, no_bias)))),
  # simulated with the conditions merged: each person has one set of
  # values, used in both conditions
  none       = simulate(merge_conditions(stim), c(other_params(1), list(b_v_pr = list(no_bias)))))

# the selections in each dataset
found <- map(sims, "found")

## ---- fig-runs
# Runs of same-class selections: sticky behaviour gives fewer, longer runs.
# One point per person and condition (averaged over their trials)
runs <- bind_rows(map(found, compute_runs), .id = "data") |>
  summarise(n_runs = mean(n_runs),
            longest_run = mean(max_run_length),
            .by = c(data, person, condition))

ggplot(runs, aes(n_runs, longest_run, colour = condition)) +
  theme_fomo() +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~data) +
  labs(x = "number of runs (mean per trial)", y = "longest run (mean per trial)")

## ---- fit
# The three models: the same multi-level model, with the data relabelled
models <- list(condition            = identity,   # the data as they are
               `person x condition` = person_by_condition,
               pooled               = merge_conditions)

# Fit every model to every dataset (nine fits, about two hours) and save
# them. centred_params = "none": with only 10 trials per person and
# condition, the model samples better this way (see prep_stan_list()).
# pathfinder_init = TRUE: start each chain from a pathfinder draw, not near
# the priors; from the priors, one chain sometimes got stuck in early warmup
if (!file.exists(fits_file)) {
  fits <- list()
  for (data in names(found)) {
    for (model in names(models)) {
      relabel <- models[[model]]
      sl <- prep_stan_list(relabel(stim), found = relabel(found[[data]]),
                           centred_params = "none")
      fits[[data]][[model]] <- fit_fomo(sl, iter = n_iter, threads_per_chain = 2,
                                        pathfinder_init = TRUE, seed = FOMO_SEED)
    }
  }
  saveRDS(list(found = found, sims = sims, fits = fits), fits_file)
}

## ---- load-fits
# Load the saved fits, and check they were made from these simulated data
saved <- readRDS(fits_file)
stopifnot("The saved fits were made from different data: delete condition_difference_fits.rds and run this script again" =
            identical(saved$found, found))
fits <- saved$fits
sims <- saved$sims

## ---- stay
# Stickiness on the model's own scale: the stay probabilities, the chance
# of picking the same class again after class 1, and after class 2 (the
# diagonal of the class-transition matrix, b_v). extract_post() gives the
# transition probabilities, with columns prev_class (the class just
# picked), class (the next one) and value (the probability).
stay_probabilities <- function(transitions) {
  transitions |>
    filter(prev_class == class) |>
    mutate(class = paste("after class", class),
           condition = as.character(condition)) |>
    rename(stay = value)
}

# the condition model's posterior for each dataset...
post <- map2(fits, found, \(f, d) extract_post(f$condition, d, stim))
# ...and the values each simulation used (the "none" data were simulated
# with the conditions merged, so they are read with the conditions merged)
truth <- list(average    = extract_post(sims$average$sim, found$average),
              individual = extract_post(sims$individual$sim, found$individual),
              none       = extract_post(sims$none$sim, merge_conditions(found$none)))

# the population stay probabilities (every posterior draw)
population <- bind_rows(map(post, \(p) stay_probabilities(filter(p$cat, dist == "posterior"))),
                        .id = "data")

# each person's stay probabilities: estimated (every posterior draw), and
# the values the simulation used
person <- bind_rows(map(post, \(p) stay_probabilities(filter(p$u$cat, dist == "posterior"))),
                    .id = "data")
person_true <- bind_rows(map(truth, \(p) stay_probabilities(p$u$cat)), .id = "data")

# in the "none" data each person's values apply to both conditions
both <- filter(person_true, condition == "both")
person_true <- bind_rows(filter(person_true, condition != "both"),
                         mutate(both, condition = "one"),
                         mutate(both, condition = "two"))

## ---- fig-stay
# The population stay probabilities, with the simulated values (dashed):
# 0.8 in condition one of the "average" data, 0.5 everywhere else
true_stay <- expand_grid(data = names(found), condition = c("one", "two")) |>
  mutate(stay = if_else(data == "average" & condition == "one", 0.8, 0.5))

ggplot(population, aes(stay, fill = class)) +
  theme_fomo() +
  ggdist::stat_halfeye(.width = FOMO_HPDI_WIDTHS, point_interval = "median_hdci", alpha = 0.6) +
  geom_vline(aes(xintercept = stay), data = true_stay, linetype = "dashed") +
  facet_grid(condition ~ data, labeller = label_both) +
  labs(x = "stay probability (population)", y = NULL, fill = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

## ---- fig-person-stay
# Each person's stay probability (the average of the two classes): the
# estimate (median and 97% interval) against the value the simulation used
person_estimates <- person |>
  summarise(stay = mean(stay), .by = c(data, condition, person, .draw)) |>
  summarise(median = median(stay),
            lower = hdci_bounds(stay, 0.97)[["lower"]],
            upper = hdci_bounds(stay, 0.97)[["upper"]],
            .by = c(data, condition, person))

person_values <- person_true |>
  summarise(true = mean(stay), .by = c(data, condition, person))

left_join(person_estimates, person_values, by = c("data", "condition", "person")) |>
  ggplot(aes(true, median, colour = condition)) +
  theme_fomo() +
  geom_abline(linetype = "dotted") +
  geom_pointrange(aes(ymin = lower, ymax = upper), alpha = 0.7) +
  facet_wrap(~data) +
  labs(x = "true stay probability (simulation)", y = "estimated stay probability")

## ---- fig-difference
# (1) How big is each difference? Condition two minus condition one, for
# the two stay probabilities, and for distance attraction (rho_delta),
# which is the same in both conditions in every dataset
stay_difference <- population |>
  select(data, .draw, condition, class, stay) |>
  pivot_wider(names_from = condition, values_from = stay) |>
  mutate(parameter = class, difference = two - one)

rho_delta_difference <- bind_rows(map(post, "spatial"), .id = "data") |>
  filter(dist == "posterior", parameter == "rho_delta") |>
  select(data, .draw, condition, value) |>
  pivot_wider(names_from = condition, values_from = value) |>
  mutate(parameter = "rho_delta", difference = two - one)

differences <- bind_rows(stay_difference, rho_delta_difference)

# the simulated differences: 0.5 - 0.8 for the stay probabilities in the
# "average" data, 0 everywhere else
true_difference <- expand_grid(data = names(found), parameter = unique(differences$parameter)) |>
  mutate(difference = if_else(data == "average" & parameter != "rho_delta", 0.5 - 0.8, 0))

ggplot(differences, aes(difference)) +
  theme_fomo() +
  ggdist::stat_halfeye(.width = FOMO_HPDI_WIDTHS, point_interval = "median_hdci") +
  geom_vline(aes(xintercept = difference), data = true_difference, linetype = "dashed") +
  facet_grid(data ~ parameter, scales = "free_x") +
  labs(x = "condition two minus condition one (population)", y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

## ---- compare
# (2) Which model predicts each dataset best? fomo_compare() gives loo and
# accuracy for the three models fitted to one dataset
comparison <- bind_rows(map2(fits, found, \(f, d) fomo_compare(f, list(stim = stim, found = d),
                                                               n_draws = 200)),
                        .id = "data")

format_compare(comparison, group = "data",
               columns = c("data", "model", "deterministic", "elpd_diff", "se_diff", "robust"))
