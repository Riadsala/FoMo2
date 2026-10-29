##########################################################################
# kfold_fomo() on simulated data:
#   1. single-level spatial, simulated from Fourier-124 with a clear 4th
#      harmonic: k-fold gives one finite elpd per trial, every fold holds
#      every condition, the Fourier-124 k-fold elpd is close to its loo
#      from a fit to all the data, and loo_compare() doesn't prefer
#      Fourier-12 (which can't make the 4th harmonic) by more than 2 SE
#   2. multi-level cat (3 people): every fold holds every person and
#      condition, one finite elpd per trial; too large a K is refused
#   3. prep_stan_list_like() (each fold's data list from a fit's own):
#      on the full data it reproduces the original data list exactly
#      (spatial with xy_scaling 7, cat with priors and harmonics,
#      multi-level cat von Mises with dir_offset), no fits
# Fits 3 + 3 + 1 + 3 small models (a few minutes). Exits 1 on any failure.
#   Rscript tests/test_kfold_fomo.R
##########################################################################

root <- rprojroot::find_root(rprojroot::has_dir("models"))
setwd(file.path(root, "tests"))
suppressMessages(source("../functions/set_up_FoMo.R"))

fit_args <- list(iter = 200, chains = 2, threads_per_chain = 1, save_dir = tempdir(),
                 refresh = 0, show_messages = FALSE)
kfold <- function(d, K, prep_args) {
  suppressMessages(do.call(kfold_fomo, c(list(d, K = K, prep_args = prep_args), fit_args)))
}
failures <- character(0)
report <- function(label, expr) {
  res <- tryCatch(expr, error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(sprintf("%-28s %s\n", label, res))
  if (startsWith(res, "FAILED")) failures <<- c(failures, label)
}

################################################################
# 1. single-level spatial, Fourier-12 vs Fourier-124
################################################################

report("single-level spatial", {
  spec <- list(K = 1, nrow = 8, ncol = 8, prop_missing = 0.25, v_type = "none")
  stim <- generate_stim(spec, n_trials = 24)
  sp <- list(rho_delta_mu = 1.5, rho_psi_mu = -0.5,
             #                    h1: a, b   h2: a, b   h4: a, b
             fourier_coef_mu = rbind(c(0.8, 0,   0.3, 0,    0.8, 0)))
  sl <- suppressMessages(prep_stan_list(stim, sp, directional_model = "fourier"))
  d <- list(stim = stim, found = extract_pred(suppressMessages(sim_fomo(sl, iter = 1)), stim))

  k12  <- kfold(d, 3, list(directional_model = "fourier", harmonics = c(1, 2)))
  k124 <- kfold(d, 3, list(directional_model = "fourier", harmonics = c(1, 2, 4)))
  n_trials <- n_distinct(d$stim$trial)

  # loo from a fit to all the data, for comparison
  sl_all <- suppressMessages(prep_stan_list(d$stim, found = d$found, directional_model = "fourier",
                                            save_selections = FALSE))
  fit_all <- suppressMessages(do.call(fit_fomo, c(list(sl_all), fit_args)))
  elpd_loo <- fit_all$loo()$estimates["elpd_loo", "Estimate"]
  elpd_kf <- k124$elpd$estimates["elpd", "Estimate"]

  # Fourier-124 minus Fourier-12 (positive: Fourier-124 better); loo_compare()
  # gives the SE of the difference on the row of the worse model
  cmp <- loo::loo_compare(list(fourier12 = k12$elpd, fourier124 = k124$elpd))
  diff_124 <- k124$elpd$estimates["elpd", "Estimate"] - k12$elpd$estimates["elpd", "Estimate"]
  se_diff <- max(cmp[, "se_diff"])
  stopifnot("not one elpd per trial" = nrow(k124$elpd$pointwise) == n_trials,
            "non-finite elpd" = all(is.finite(k124$elpd$pointwise[, "elpd"])),
            "a fold misses a condition" = all(count(k124$folds, fold, condition)$n > 0),
            "k-fold elpd far from loo (over 10%)" = abs(elpd_kf - elpd_loo) < 0.1 * abs(elpd_loo),
            "k-fold prefers Fourier-12 by more than 2 SE" = diff_124 > -2 * se_diff)
  sprintf("k-fold %.1f vs loo %.1f; Fourier-124 - Fourier-12 = %.1f (SE %.1f)",
          elpd_kf, elpd_loo, diff_124, se_diff)
})

################################################################
# 2. multi-level cat: folds hold every person
################################################################

report("multi-level cat", {
  spec <- list(K = 2, nrow = 5, ncol = 5, prop_missing = 0, v_type = "cat", n_classes = 2)
  stim <- generate_stim(spec, n_trials = 3, n_people = 3)
  sl <- suppressMessages(prep_stan_list(stim))
  d <- list(stim = stim, found = extract_pred(suppressMessages(sim_fomo(sl, iter = 1)), stim))

  k <- kfold(d, 3, list(centred_params = "none"))
  stopifnot("not one elpd per trial" = nrow(k$elpd$pointwise) == n_distinct(d$stim$trial),
            "non-finite elpd" = all(is.finite(k$elpd$pointwise[, "elpd"])),
            "a fold misses a person x condition" =
              nrow(count(k$folds, fold, person, condition)) == 3 * 3 * 2,
            "K larger than the trials per cell not refused" =
              inherits(try(kfold_fomo(d, K = 4), silent = TRUE), "try-error"))
  sprintf("elpd %.1f (SE %.1f), %d folds; max R-hat %.2f",
          k$elpd$estimates["elpd", "Estimate"], k$elpd$estimates["elpd", "SE"],
          nrow(k$fits), max(k$fits$max_rhat))
})

################################################################
report("prep_stan_list_like()", {
  same <- function(a, b) identical(lapply(a[sort(names(a))], unname), lapply(b[sort(names(b))], unname))
  sim <- function(spec, n_trials, n_people = 1) {
    stim <- generate_stim(spec, n_trials = n_trials, n_people = n_people)
    sl <- suppressMessages(prep_stan_list(stim))
    list(stim = stim, found = extract_pred(suppressMessages(sim_fomo(sl, iter = 1)), stim))
  }
  cases <- list(
    spatial = list(d = sim(list(K = 2, nrow = 5, ncol = 5, prop_missing = 0.2, v_type = "none"), 4),
                   args = list(xy_scaling = 7)),
    cat = list(d = sim(list(K = 1, nrow = 5, ncol = 5, prop_missing = 0, v_type = "cat", n_classes = 2), 4),
               args = list(harmonics = c(1, 2), priors = list(rho_delta_mu = 2))),
    cat_multi = list(d = sim(list(K = 2, nrow = 5, ncol = 5, prop_missing = 0.2, v_type = "cat", n_classes = 3), 3, 3),
               args = list(directional_model = "von_mises", dir_offset = pi / 4, centred_params = "none")))
  ok <- map_lgl(cases, function(cs) {
    sl <- suppressMessages(do.call(prep_stan_list, c(list(cs$d$stim, found = cs$d$found), cs$args)))
    same(sl, prep_stan_list_like(sl, cs$d))
  })
  stopifnot("prep_stan_list_like() doesn't reproduce the data list" = all(ok))
  paste("reproduced:", paste(names(ok), collapse = ", "))
})

if (length(failures) == 0) {
  cat("\nALL kfold_fomo() TESTS PASSED\n")
} else {
  cat("\n", length(failures), "FAILURE(S):", paste(failures, collapse = ", "), "\n")
  quit(status = 1)
}
