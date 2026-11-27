##########################################################################
# Smoke test for every production Stan model. A correctness check, not a
# benchmark or parameter-recovery test.
#
# Stage 1: compile all 16 models in models/{single,multi}_level/
#          {fit,simulate}/{von_mises,fourier}/{spatial,cat}.stan
#          (fit models with stan_threads + O1, as load_stan_model() does)
# Stage 2: for each of the 8 {single,multi}_level x {von_mises,fourier} x
#          {spatial,cat} combinations, on a tiny generated stimulus set:
#          prep_stan_list() -> sim_fomo() -> extract_pred() ->
#          prep_stan_list(found = ...) -> fit_fomo() (a few iterations; the
#          default centring, FOMO_CENTRED_DEFAULT, for multi-level),
#          then check that log_lik came back finite and that extract_post() /
#          fomo_summary() work on the fit.
# Stage 3: the 4 multi-level combinations again with every per-person
#          effect centred (prep_stan_list(centred_params = "all")), plus a
#          posterior predictive simulation from that fit (sim_fomo(sl, fit))
#          and prep_stan_list()'s centred_params errors and default flags.
# Stage 4: multi-level sim_params defaults and prior overrides (the
#          Fourier fit starting from pathfinder, pathfinder_init = TRUE)
# Stage 5: Fourier with every sine term (cos_only = NULL; the default
#          leaves out harmonic 4's), each level and item type
# Stage 6: cat_sym.stan (prep_stan_list(symmetric_v = TRUE), multi-level
#          Fourier cat fits only): data simulated from cat.stan, fitted
#          with the default centring and with every effect centred; the
#          save_selections / predict_fomo() checks and extract_post()
# Stage 2 also checks save_selections: with it off, the same draws
# give the same per-trial log_lik and no P, W or loglik, and
# predict_fomo() gets them back
#
# Run after any change to models/, models/include/ or the model-facing
# functions (prep.R, simulate.R, fit.R, extract.R, summarise.R):
#   Rscript tests/smoke_test_models.R
# Exits with status 1 if anything fails, so it can be used in scripts.
# Takes a few minutes once models are compiled (longer on first compile).
##########################################################################

root <- rprojroot::find_root(rprojroot::has_dir("models"))
setwd(file.path(root, "tests"))
source("../functions/set_up_FoMo.R")

levels_ <- c("single_level", "multi_level")
dir_models <- c("von_mises", "fourier")
types <- c("spatial", "cat")

failures <- character(0)

################################################################
# Stage 1: compile everything
################################################################

cat("\n=== Stage 1: compile all models ===\n")

for (lv in levels_) for (task in c("fit", "simulate")) for (dm in dir_models) for (ty in types) {
  path <- file.path(root, "models", lv, task, dm, paste0(ty, ".stan"))
  label <- paste(lv, task, dm, ty, sep = "/")
  threaded <- task == "fit"
  ok <- tryCatch({
    cmdstan_model(path,
                  cpp_options = if (threaded) list(stan_threads = TRUE) else list(),
                  stanc_options = if (threaded) list("O1") else list(),
                  quiet = TRUE)
    TRUE
  }, error = function(e) { cat("  ", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("%-45s %s\n", label, if (ok) "OK" else "FAILED"))
  if (!ok) failures <- c(failures, paste("compile", label))
}

################################################################
# Stage 2: simulate -> fit on tiny data
################################################################

cat("\n=== Stage 2: simulate -> fit ===\n")

# save_selections: the fit model's generated quantities on the same draws
# with the flag on and off, and the per-selection loglik sums to each
# trial's log_lik
check_save_selections <- function(m_fit, sl_fit) {
  sl_fit$grainsize <- 1
  mod <- load_stan_model(sl_fit, model_task = "fit")
  # the fit's parameter draws as predict_fomo() passes them (sum-to-zero
  # vectors re-centred after CmdStan's rounding)
  dr <- parameter_draws(m_fit, mod, sl_fit)
  gq <- function(flag) {
    sl <- sl_fit; sl$save_selections <- flag
    mod$generate_quantities(fitted_params = dr, data = sl, threads_per_chain = 1, seed = 1)
  }
  g_on <- gq(1L); g_off <- gq(0L)
  ll_on <- unclass(g_on$draws("log_lik", format = "matrix"))
  stopifnot("save_selections changes log_lik" =
              max(abs(ll_on - unclass(g_off$draws("log_lik", format = "matrix")))) < 1e-8,
            "save_selections = 0 still outputs P/W/loglik" =
              # size-0 variables are left out of cmdstanr's metadata
              all(map_dbl(c("P", "W", "loglik"),
                          \(v) prod(g_off$metadata()$stan_variable_sizes[[v]] %||% 0)) == 0))
  sel <- unclass(g_on$draws("loglik", format = "matrix"))
  sums <- sapply(split(seq_len(ncol(sel)), sl_fit$trial), \(i) rowSums(sel[, i, drop = FALSE]))
  stopifnot("per-selection loglik doesn't sum to log_lik" = max(abs(sums - ll_on)) < 1e-6)
  sl_off <- sl_fit; sl_off$save_selections <- 0L
  g <- predict_fomo(m_fit, sl_off, n_draws = 3, threads_per_chain = 1, seed = 2)
  stopifnot("predict_fomo() didn't return W for every selection" =
              prod(g$metadata()$stan_variable_sizes$W) == sl_fit$Ns)
}

specs <- list(
  spatial = list(K = 1, nrow = 4, ncol = 4, prop_missing = 0, v_type = "none"),
  cat     = list(K = 1, nrow = 4, ncol = 4, prop_missing = 0, v_type = "cat", n_classes = 2))

run_case <- function(lv, dm, ty, centred_params = "default") {
  n_people <- if (lv == "multi_level") 3 else 1
  n_trials <- if (lv == "multi_level") 2 else 3
  stim <- generate_stim(specs[[ty]], n_trials = n_trials, n_people = n_people)

  sl <- prep_stan_list(stim, directional_model = dm)
  # iter = 1: more draws would duplicate the trial data in extract_pred()
  m_sim <- sim_fomo(sl, iter = 1, chains = 1)
  df <- extract_pred(m_sim, stim)

  sl_fit <- prep_stan_list(stim, found = df, directional_model = dm, centred_params = centred_params)
  m_fit <- fit_fomo(sl_fit, iter = 5, chains = 1, save_dir = tempdir(),
                    grainsize = 1, threads_per_chain = 1)
  if (identical(centred_params, "all")) {
    stopifnot("centred_params = \"all\" didn't set every centre_* flag" =
                all(unlist(sl_fit[str_starts(names(sl_fit), "centre_")]) == 1))
    # the simulate models (always non-centred) read the fit's draws by position
    extract_pred(sim_fomo(sl_fit, m_fit, iter = 1), stim)
  }

  ll <- m_fit$draws("log_lik", format = "matrix")
  stopifnot("log_lik missing or not finite" = length(ll) > 0 && all(is.finite(ll)))
  if (identical(centred_params, "default")) check_save_selections(m_fit, sl_fit)

  # post-fit helpers must work for both directional models
  post <- extract_post(m_fit, df)
  stopifnot("extract_post() returned no spatial parameters" = !is.null(post$spatial))
  stopifnot("post$vonmises should be present only for von Mises fits" =
              is.null(post$vonmises) == (dm == "fourier"))
  stopifnot("post$fourier should be present only for Fourier fits" =
              is.null(post$fourier) == (dm == "von_mises"))
  # the first selection's class weights (the stationary distribution of
  # b_v): population-level in single-level cat fits, per person in
  # multi-level ones
  stopifnot("post$first_v should be present only for single-level cat fits" =
              is.null(post$first_v) == (ty != "cat" || lv == "multi_level"))
  if (dm == "fourier") {
    stopifnot("Fourier harmonics not labelled from the fit" =
                all(str_detect(levels(post$fourier$harmonic), "^h\\d+$")),
              "default Fourier fit should leave out harmonic 4's sine (cos_only)" =
                !any(post$fourier$harmonic == "h4" & post$fourier$term == "sin"))
  }
  if (lv == "multi_level") {
    stopifnot("per-person Fourier coefficients missing" =
                is.null(post$u$fourier) == (dm == "von_mises"))
    stopifnot("per-person first_v missing" = is.null(post$u$first_v) == (ty != "cat"))
  }
  if (ty == "cat") {
    first_v <- if (lv == "multi_level") post$u$first_v else post$first_v
    sums <- first_v |> group_by(pick(any_of(c(".draw", "condition", "person")))) |>
      summarise(s = sum(value), .groups = "drop")
    stopifnot("first_v probabilities don't sum to 1" = all(abs(sums$s - 1) < 1e-6))
  }
  fomo_summary(m_fit)
  invisible(TRUE)
}

for (lv in levels_) for (dm in dir_models) for (ty in types) {
  label <- paste(lv, dm, ty, sep = "/")
  ok <- tryCatch({
    suppressMessages(capture.output(run_case(lv, dm, ty)))
    TRUE
  }, error = function(e) { cat("  ", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("%-45s %s\n", label, if (ok) "OK" else "FAILED"))
  if (!ok) failures <- c(failures, paste("sim/fit", label))
}

################################################################
# Stage 3: multi-level with every per-person effect centred
################################################################

cat("\n=== Stage 3: multi-level, centred_params = \"all\" ===\n")

for (dm in dir_models) for (ty in types) {
  label <- paste("multi_level", dm, ty, "centred", sep = "/")
  ok <- tryCatch({
    suppressMessages(capture.output(run_case("multi_level", dm, ty, centred_params = "all")))
    TRUE
  }, error = function(e) { cat("  ", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("%-45s %s\n", label, if (ok) "OK" else "FAILED"))
  if (!ok) failures <- c(failures, paste("centred", label))
}

# prep_stan_list() must refuse a group the model doesn't have, and
# centring for one-person data (no multi-level model)
expect_error <- function(expr) inherits(tryCatch(suppressMessages(expr), error = identity), "error")
ok <- expect_error(prep_stan_list(generate_stim(specs$spatial, n_trials = 2, n_people = 3), centred_params = "v")) &&
  expect_error(prep_stan_list(generate_stim(specs$cat, n_trials = 2, n_people = 1), centred_params = "all"))

# ... and the default centres FOMO_CENTRED_DEFAULT's groups (delta and v)
# where the model has them, and is quietly ignored for one-person data
flags <- function(ty, n_people = 3) {
  sl <- suppressMessages(prep_stan_list(generate_stim(specs[[ty]], n_trials = 2, n_people = n_people)))
  unlist(sl[str_starts(names(sl), "centre_")])
}
ok <- ok &&
  identical(names(which(flags("spatial") == 1)), "centre_delta") &&
  identical(names(which(flags("cat") == 1)), c("centre_delta", "centre_v")) &&
  is.null(flags("cat", n_people = 1))
cat(sprintf("%-45s %s\n", "centred_params errors", if (ok) "OK" else "FAILED"))
if (!ok) failures <- c(failures, "centred_params errors")

################################################################
# Stage 4: sim_params that leave out the directional parameters and their
# sd_lambda (prep_stan_list() must fill in defaults), then a fit with one
# prior overridden
################################################################

cat("\n=== Stage 4: multi-level, sim_params and prior overrides ===\n")

for (dm in dir_models) {
  label <- paste("multi_level", dm, "spatial, sim_params", sep = "/")
  ok <- tryCatch({
    suppressMessages(capture.output({
      stim <- generate_stim(specs$spatial, n_trials = 2, n_people = 3)
      simparams <- list(rho_delta_mu = 1, rho_psi_mu = 0,
                        sd_lambda = list(rho_delta = 3, rho_psi = 3))
      sl <- prep_stan_list(stim, simparams, directional_model = dm)
      m_sim <- sim_fomo(sl, iter = 1, chains = 1)
      df <- extract_pred(m_sim, stim)
      sl_fit <- prep_stan_list(stim, found = df, priors = list(rho_delta_sd = 0.5),
                               directional_model = dm)
      stopifnot(all(sl_fit$rho_delta_sd == 0.5))
      # assigned, so capture.output() doesn't print the fit (cmdstanr's
      # print of a 5-iteration fit warns "NAs introduced by coercion")
      # Fourier: from pathfinder's draws (fit_fomo(pathfinder_init = TRUE)),
      # whose sum-to-zero vectors must be re-centred to be valid starts
      m_fit <- fit_fomo(sl_fit, iter = 5, chains = 2, save_dir = tempdir(),
               grainsize = 1, threads_per_chain = 1, pathfinder_init = dm == "fourier")
    }))
    TRUE
  }, error = function(e) { cat("  ", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("%-45s %s\n", label, if (ok) "OK" else "FAILED"))
  if (!ok) failures <- c(failures, paste("sim_params", label))
}

################################################################
# Stage 5: Fourier with every sine term (prep_stan_list(cos_only = NULL);
# stages 2-4 use the default, harmonic 4 cosine only): simulate -> fit,
# the coefficients' labels, and predict_fomo()
################################################################

cat("\n=== Stage 5: Fourier, every sine term (cos_only = NULL) ===\n")

for (lv in levels_) for (ty in types) {
  label <- paste(lv, "fourier", ty, "all sines", sep = "/")
  ok <- tryCatch({
    suppressMessages(capture.output({
      n_people <- if (lv == "multi_level") 3 else 1
      stim <- generate_stim(specs[[ty]], n_trials = 2, n_people = n_people)
      sl <- prep_stan_list(stim, list(fourier_coef_mu = c(0.5, 0, 0, 0.3, 0.4, 0)), cos_only = NULL)
      df <- extract_pred(sim_fomo(sl, iter = 1, chains = 1), stim)
      sl_fit <- prep_stan_list(stim, found = df, cos_only = NULL)
      m_fit <- fit_fomo(sl_fit, iter = 5, chains = 1, save_dir = tempdir(),
                        grainsize = 1, threads_per_chain = 1)
      post <- extract_post(m_fit, df)
      stopifnot("all sines: wrong Fourier coefficients" =
                  identical(paste(post$fourier$harmonic, post$fourier$term) |> unique(),
                            c("h1 cos", "h1 sin", "h2 cos", "h2 sin", "h4 cos", "h4 sin")),
                "all sines: log_lik not finite" = all(is.finite(m_fit$draws("log_lik"))))
      compute_directional_curve(post)
      g <- predict_fomo(m_fit, n_draws = 2, threads_per_chain = 1, seed = 2)
    }))
    TRUE
  }, error = function(e) { cat("  ", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("%-45s %s\n", label, if (ok) "OK" else "FAILED"))
  if (!ok) failures <- c(failures, paste("all sines", label))
}

################################################################
# Stage 6: cat_sym.stan, the symmetric transition matrix (multi-level
# Fourier cat fits only; no simulate model, so the data come from cat.stan)
################################################################

cat("\n=== Stage 6: multi-level Fourier cat_sym (symmetric_v = TRUE) ===\n")

for (centred in c("default", "all")) {
  label <- paste("multi_level/fourier/cat_sym", if (centred == "all") "centred" else "", sep = "/") |>
    str_remove("/$")
  ok <- tryCatch({
    suppressMessages(capture.output({
      stim <- generate_stim(specs$cat, n_trials = 2, n_people = 3)
      df <- extract_pred(sim_fomo(prep_stan_list(stim), iter = 1, chains = 1), stim)
      sl_fit <- prep_stan_list(stim, found = df, symmetric_v = TRUE, centred_params = centred)
      stopifnot("symmetric_v didn't set the flag" = sl_fit$symmetric_v == 1)
      m_fit <- fit_fomo(sl_fit, iter = 5, chains = 1, save_dir = tempdir(),
                        grainsize = 1, threads_per_chain = 1)
      stopifnot("cat_sym: no stay log-odds in the fit" = "b_stay" %in% m_fit$metadata()$stan_variables,
                "cat_sym: log_lik missing or not finite" = all(is.finite(m_fit$draws("log_lik"))))
      check_save_selections(m_fit, sl_fit)
      post <- extract_post(m_fit, df)
      stopifnot("cat_sym: no Fourier coefficients" = !is.null(post$fourier))
    }))
    TRUE
  }, error = function(e) { cat("  ", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("%-45s %s\n", label, if (ok) "OK" else "FAILED"))
  if (!ok) failures <- c(failures, paste("cat_sym", label))
}

################################################################
# Summary
################################################################

if (length(failures) == 0) {
  cat("\nALL SMOKE TESTS PASSED (16 compiles, 8 simulate -> fit runs, 4 centred runs, 2 sim_params runs, 4 all-sines runs, 2 cat_sym runs)\n")
} else {
  cat("\n", length(failures), "FAILURE(S):\n", paste(" -", failures, collapse = "\n"), "\n")
  quit(status = 1)
}
