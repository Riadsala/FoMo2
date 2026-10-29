##########################################################################
# Multi-level cat (categorical item features) model: simulate data for
# several people with known parameters, fit the model, check parameter
# recovery at the condition and per-person level, and run posterior
# predictive checks. cat.qmd renders this script; the "## ----" lines
# mark its chunks.
##########################################################################

## ---- load
# import packages and load our functions
source("../../../../functions/set_up_FoMo.R")

# sampler settings: warmup and sampling iterations per chain (the same in
# every example)
n_iter <- 500
n_post_draws <- 50   # posterior draws for posterior predictive simulation

## ---- stimuli
spec <- list(
  K = 2,
  nrow = 8, ncol = 8,
  prop_missing = 0.25,
  v_type = "cat",
  n_classes = 2)

stim <- generate_stim(spec, n_trials = 10, n_people = 6)

## ---- simulate
simparams_cat <- list(
  b_v_pr = list(matrix(c(0.85, 0.15,
                         0.45, 0.55),
                       nrow = spec$n_classes),
                matrix(c(0.50, 0.50,
                         0.20, 0.80),
                       nrow = spec$n_classes)))

simparams_spa <- list(
  rho_delta_mu = c(0.75, 1.50),
  rho_psi_mu = c(-1, 1),
  # absolute direction, Fourier-124 ({cos, sin} per harmonic): condition
  # one prefers the four cardinal directions (h4), two rightwards (h1)
  #                       h1: a, b   h2: a, b   h4: a, b
  fourier_coef_mu = rbind(c(0,   0,     0, 0,     0.8, 0),
                          c(1.5, 0,     0, 0,     0,   0)))

# random effects: each SD is drawn from exponential(sd_lambda)
simparams_sd <- list(
  sd_lambda = list(
    b_v = 2,
    rho_delta = c(3, 3),
    rho_psi = c(10, 10),
    fourier_coef = c(4, 4, 4)))

simparams <- c(simparams_sd, simparams_spa, simparams_cat)

sl_sim <- prep_stan_list(stim, simparams)

sim <- sim_fomo(sl_sim, iter = 1, chains = 1)

found <- extract_pred(sim, stim)

## ---- fig-trial
plot_trial(found, stim)

## ---- fig-runs
plot_runs(compute_runs(found))

## ---- trial-stats
trial_stats <- compute_trial_stats(found)

head(trial_stats)

## ---- groundtruth
# each person's simulated values
post <- extract_post(sim, found)
groundtruth <- groundtruth_from_sim(simparams, post)

## ---- fig-check-bv
ggplot(post$u$cat, aes(prev_class, value, fill = class)) +
  theme_fomo() +
  geom_boxplot(alpha = 0.5) +
  geom_point(data = post$cat, aes(colour = class), size = 5) +
  facet_wrap(~class)

## ---- fit_prep
sl <- prep_stan_list(stim, found = found,
                     centred_params = "none",
                     save_selections = TRUE)

## ---- fit
# a fixed seed, so the chains are the same on every run
fit <- fit_fomo(sl, iter = n_iter, chains = 4, seed = FOMO_SEED)

## ---- fit-summary
# estimates and convergence, and parameter recovery: are the values used
# to simulate the data inside the posterior's 97% HPDIs?
fomo_summary(fit, sim, random_effects = TRUE, cond_labels = levels(stim$condition))

## ---- fig-traces
plot_post_traces(fit)

## ---- loo
loo_result <- fit$loo()
loo_result

# the trials with a high Pareto k, worst first
check_high_pareto_k(loo_result, sl)

## ---- fig-post
post <- extract_post(fit, found, stim)
plot_post(post, groundtruth = groundtruth)

## ---- fig-person-spatial
plot_post_per_person_spatial(post, groundtruth)

## ---- fig-accuracy
plot_post_accuracy(post$acc, stim, found)

## ---- sim-posterior
sim_post <- sim_fomo(sl, fit, iter = n_post_draws)

# the model doesn't predict when a trial stops: keep as many simulated
# selections per trial as the data have
extract_pred(sim_post, stim) |>
  semi_join(found, by = c("trial", "n")) -> found_post

## ---- fig-sim-trial
plot_trial(found_post, stim)

## ---- fig-sim-feature
plot_features(found, feature = "v", sim = found_post)

## ---- fig-sim-stats
# compare per-trial foraging statistics between the training data and data
# simulated from the fitted posterior
sim_stats <- compute_trial_stats(found_post)
plot_trial_stats(trial_stats, sim_stats)
