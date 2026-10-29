##########################################################################
# predict_fomo() and the accuracy rule, on simulated data: single- and
# multi-level spatial and cat models (Fourier-124). For each:
#   1. on a fit that stores P, W and loglik (save_selections = TRUE),
#      predict_fomo() with every draw gives the same W and loglik per
#      selection (averaged over draws, so draw order doesn't matter) and
#      the same stochastic accuracy (P) up to simulation noise
#   2. on a fit without them (save_selections = FALSE) it returns them for
#      every selection, and extract_post(pred, found, stim) has accuracy
#   3. the accuracy rule: the default ("first_and_forced") leaves out one
#      first pick per trial plus every forced last pick, "none" leaves out
#      nothing, and without stim it says it can't find forced picks
#   4. deterministic accuracy beats chance
# and plot_post_accuracy() draws with its baselines on the kept selections.
# A few minutes (8 small fits). Exits 1 on any failure.
#   Rscript tests/test_predict_fomo.R
##########################################################################

root <- rprojroot::find_root(rprojroot::has_dir("models"))
setwd(file.path(root, "tests"))
suppressMessages(source("../functions/set_up_FoMo.R"))

specs <- list(
  spatial = list(K = 1, nrow = 5, ncol = 5, prop_missing = 0, v_type = "none"),
  cat     = list(K = 1, nrow = 5, ncol = 5, prop_missing = 0, v_type = "cat", n_classes = 2))

per_selection <- function(gq) {
  # mean over draws of each selection's accuracy and loglik
  gq |>
    summarise(det = mean(W == id), sto = mean(P == id), ll = mean(loglik),
              .by = c(trial, n))
}

run_case <- function(lv, ty) {
  n_people <- if (lv == "multi_level") 3 else 1
  n_trials <- if (lv == "multi_level") 3 else 6
  stim <- generate_stim(specs[[ty]], n_trials = n_trials, n_people = n_people)
  sl <- suppressMessages(prep_stan_list(stim))
  df <- extract_pred(suppressMessages(sim_fomo(sl, iter = 1, chains = 1)), stim)

  quiet_fit <- function(sl, iter) {
    suppressMessages(fit_fomo(sl, iter = iter, chains = 2, save_dir = tempdir(),
                              threads_per_chain = 1, refresh = 0, show_messages = FALSE))
  }

  # 1. stored vs predicted, every draw
  sl_on <- suppressMessages(prep_stan_list(stim, found = df, centred_params = "none",
                                           save_selections = TRUE))
  fit_on <- quiet_fit(sl_on, 200)
  n_draws <- posterior::ndraws(fit_on$draws("lp__"))
  pred <- predict_fomo(fit_on, sl_on, n_draws = n_draws, threads_per_chain = 1, seed = 3)
  vars_fit <- get_variables(fit_on); vars_pred <- get_variables(pred)
  stored <- per_selection(extract_gq(fit_on, df, vars_fit, stim))
  redone <- per_selection(extract_gq(pred, df, vars_pred, stim))
  stopifnot("predict_fomo() doesn't use every draw" = n_distinct(extract_gq(pred, df, vars_pred)$.draw) == n_draws,
            "predicted W differs from the stored W" = max(abs(stored$det - redone$det)) < 1e-9,
            "predicted loglik differs from the stored loglik" = max(abs(stored$ll - redone$ll)) < 1e-4,
            "stochastic accuracy differs by more than noise" = abs(mean(stored$sto) - mean(redone$sto)) < 0.03)

  # 2. a fit without per-selection output
  sl_off <- suppressMessages(prep_stan_list(stim, found = df, centred_params = "none",
                                            save_selections = FALSE))
  fit_off <- quiet_fit(sl_off, 20)
  stopifnot("save_selections = FALSE fit still has W" =
              is.null(fit_off$metadata()$stan_variable_sizes$W))
  pred_off <- predict_fomo(fit_off, sl_off, n_draws = 10, threads_per_chain = 1)
  stopifnot("predict_fomo() didn't return W for every selection" =
              prod(pred_off$metadata()$stan_variable_sizes$W) == nrow(df),
            "extract_post(pred, found, stim) has no accuracy" =
              !is.null(suppressMessages(extract_post(pred_off, df, stim))$acc))
  # the accuracy plot, its baselines on the same selections as accuracy
  post_off <- suppressMessages(extract_post(pred_off, df, stim))
  stopifnot("first picks in post$acc" = !any(post_off$acc$n == 1))
  invisible(ggplot2::ggplot_build(plot_post_accuracy(post_off$acc, stim, df)))

  # 3. the accuracy rule
  gq <- extract_gq(pred, df, vars_pred, stim)
  per_draw <- function(x) nrow(x) / n_distinct(x$.draw)
  n_trials_total <- n_distinct(df$trial)
  n_forced <- df |> left_join(count(stim, trial, name = "n_items"), by = "trial") |>
    filter(n == n_items) |> nrow()
  stopifnot("accuracy_exclude \"none\" left something out" = per_draw(keep_for_accuracy(gq, "none")) == nrow(df),
            "default rule: wrong number of selections left out" =
              per_draw(keep_for_accuracy(gq)) == nrow(df) - n_trials_total - n_forced,
            "default rule kept a first or forced pick" =
              !any(keep_for_accuracy(gq)$first | keep_for_accuracy(gq)$forced))
  msg <- character(0)
  withCallingHandlers(keep_for_accuracy(extract_gq(pred, df, vars_pred)),
                      message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") })
  stopifnot("no message without stim" = any(grepl("forced last picks can't be identified", msg)))
  stopifnot("unknown rule not caught" = inherits(try(keep_for_accuracy(gq, "nonsense"), silent = TRUE), "try-error"))

  # 4. better than chance
  kept <- keep_for_accuracy(gq)
  acc <- c(deterministic = mean(kept$W == kept$id), chance = mean(kept$chance))
  stopifnot("deterministic accuracy not above chance" = acc[["deterministic"]] > acc[["chance"]])
  sprintf("n = %d selections, %d forced; deterministic %.2f vs chance %.2f",
          nrow(df), n_forced, acc[["deterministic"]], acc[["chance"]])
}

failures <- character(0)
for (lv in c("single_level", "multi_level")) for (ty in names(specs)) {
  label <- paste(lv, ty, sep = "/")
  res <- tryCatch(run_case(lv, ty), error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(sprintf("%-25s %s\n", label, res))
  if (startsWith(res, "FAILED")) failures <- c(failures, label)
}

if (length(failures) == 0) {
  cat("\nALL predict_fomo() TESTS PASSED\n")
} else {
  cat("\n", length(failures), "FAILURE(S):", paste(failures, collapse = ", "), "\n")
  quit(status = 1)
}
