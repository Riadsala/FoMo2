##########################################################################
# von Mises vs Fourier directional models (single-level spatial, one
# condition): simulate one dataset from each, fit both models to both
# datasets, compare their directional curves and parameter recovery, then
# compare the models with fomo_compare() (accuracy, loo, weights,
# robustness) and, once, with K-fold. fourier_v_vonmises.qmd renders this
# script; the "## ----" lines mark its chunks.
##########################################################################

## ---- load
source("../../../../functions/set_up_FoMo.R")

# sampler settings: warmup and sampling iterations per chain (the same in
# every example)
n_iter <- 500

## ---- stimuli
# one condition, 100 trials of 40 items (an 8 x 8 grid, 24 cells empty),
# as in Kristjansson et al. (2014)
spec <- list(K = 1, nrow = 8, ncol = 8, prop_missing = 0.375, v_type = "none")
stim <- generate_stim(spec, n_trials = 100)

## ---- simulate_vm
# von Mises: one log-weight (log_theta) each for right, up, left and down.
# Values: the von Mises fit to Kristjansson et al.'s feature condition
sl_sim_vm <- prep_stan_list(stim, 
                            list(rho_delta_mu = 2.41, 
                                 rho_psi_mu = -2.14,
                                 log_theta_mu = rbind(c(-3.07, -0.10, -5.07, -0.03))),
                            directional_model = "von_mises")
sim_vm <- sim_fomo(sl_sim_vm, iter = 1, chains = 1)
found_vm <- extract_pred(sim_vm, stim)

## ---- simulate_f
# Fourier-124: {cos, sin} coefficients for harmonics 1, 2 and 4. Values:
# the Fourier-124 fit to the same data
sl_sim_f <- prep_stan_list(stim, 
                           list(rho_delta_mu = 2.50, 
                                rho_psi_mu = -2.12,
                                fourier_coef_mu = rbind(c(0.27, 0.04, -0.45, 0.19, 0.16, -0.04))),
                           directional_model = "fourier")
sim_f <- sim_fomo(sl_sim_f, iter = 1, chains = 1)
found_f <- extract_pred(sim_f, stim)

## ---- fig-true-curves
truth <- list(vm = compute_directional_curve(extract_post(sim_vm, found_vm)),
              f = compute_directional_curve(extract_post(sim_f, found_f)))

bind_rows(`von Mises` = truth$vm, `Fourier-124` = truth$f, .id = "model") |>
  ggplot(aes(angle, log_weight, colour = model)) +
  theme_fomo() +
  geom_line() +
  scale_colour_model() +
  scale_x_direction() +
  labs(y = "log-weight (centred)", colour = "simulated from")

## ---- fit
# each dataset fitted with each directional model: fit_<data>_<model>
fit_with <- function(found, model) {
  fit_fomo(prep_stan_list(stim, 
                          found = found, 
                          directional_model = model), 
           iter = n_iter, seed = FOMO_SEED)
}

fit_vm_vm <- fit_with(found_vm, "von_mises")
fit_vm_f  <- fit_with(found_vm, "fourier")
fit_f_vm  <- fit_with(found_f, "von_mises")
fit_f_f   <- fit_with(found_f, "fourier")

## ---- recovery-vm
# each dataset's own model: are the simulation values inside the 97% HPDIs?
fomo_summary(fit_vm_vm, sim_vm)

## ---- recovery-f
fomo_summary(fit_f_f, sim_f)

## ---- fig-curves-vm
bind_rows(`von Mises`   = compute_directional_curve(extract_post(fit_vm_vm, found_vm, stim)),
          `Fourier-124` = compute_directional_curve(extract_post(fit_vm_f, found_vm, stim)),
          .id = "model") |>
  plot_post_directional_curves(truth$vm)

## ---- fig-curves-f
bind_rows(`von Mises`   = compute_directional_curve(extract_post(fit_f_vm, found_f, stim)),
          `Fourier-124` = compute_directional_curve(extract_post(fit_f_f, found_f, stim)),
          .id = "model") |>
  plot_post_directional_curves(truth$f)

## ---- compare
# fomo_compare() compares the models fitted to one dataset: one call each
compare_vm <- fomo_compare(list(`von Mises` = fit_vm_vm, `Fourier-124` = fit_vm_f),
                           data = list(stim = stim, found = found_vm))
compare_f  <- fomo_compare(list(`von Mises` = fit_f_vm, `Fourier-124` = fit_f_f),
                           data = list(stim = stim, found = found_f))

bind_rows(`von Mises data` = compare_vm, `Fourier-124 data` = compare_f, .id = "dataset") |>
  format_compare(columns = c("dataset", "model", "deterministic", "gain_nn", "elpd_diff",
                             "se_diff", "pct_high_k", "stacking", "robust"),
                 group = "dataset")

## ---- kfold
# the same comparison with 3-fold cross-validation (3 fits per model)
kfold_vm <- fomo_compare(list(`von Mises` = fit_vm_vm, `Fourier-124` = fit_vm_f),
                         data = list(stim = stim, found = found_vm), kfold_cv = 3)
format_compare(kfold_vm, columns = c("model", "elpd_diff", "se_diff",
                                     "elpd_diff_kfold", "se_diff_kfold"))
