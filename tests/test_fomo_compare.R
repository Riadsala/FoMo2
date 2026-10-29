##########################################################################
# fomo_compare() on simulated data, two datasets (one call each):
#   spatial: simulated from Fourier-124 with a clear 4th harmonic, fitted
#     with Fourier-124 (save_selections = FALSE, so its accuracy needs the
#     data list fit_fomo() stores with the fit) and Fourier-12
#   cat: von Mises and Fourier-124 (both storing P/W)
# Checks: one row per model; accuracy for every row; Fourier-124 best on
# the spatial data; the best model's row has elpd_diff 0 and is robust;
# the k-fold warning fires when its thresholds say so, but not for loo's
# "N < 100" (a message instead); a tibble, which format_compare() makes a
# kable; models fitted to different trials are refused; several datasets
# in one call refused; kfold_cv = K adds the k-fold columns (and no
# warnings); bad K, group and strata give clear errors. A few minutes (5
# small fits, and 4 for k-fold). Exits 1 on any failure.
#   Rscript tests/test_fomo_compare.R
##########################################################################

root <- rprojroot::find_root(rprojroot::has_dir("models"))
setwd(file.path(root, "tests"))
suppressMessages(source("../functions/set_up_FoMo.R"))

quiet_fit <- function(sl) {
  suppressMessages(fit_fomo(sl, iter = 200, chains = 2, threads_per_chain = 1, save_dir = tempdir(),
                            refresh = 0, show_messages = FALSE))
}
simulate <- function(spec, n_trials, sp = NULL, ...) {
  stim <- generate_stim(spec, n_trials = n_trials)
  sl <- suppressMessages(prep_stan_list(stim, sp, ...))
  list(stim = stim, found = extract_pred(suppressMessages(sim_fomo(sl, iter = 1)), stim))
}

# spatial, Fourier-124 vs Fourier-12
d_sp <- simulate(list(K = 1, nrow = 8, ncol = 8, prop_missing = 0.25, v_type = "none"), 20,
                 list(rho_delta_mu = 1.5, rho_psi_mu = -0.5,
                      fourier_coef_mu = rbind(c(0.8, 0, 0.3, 0, 0.8, 0))),
                 directional_model = "fourier")
f124 <- quiet_fit(suppressMessages(prep_stan_list(d_sp$stim, found = d_sp$found, save_selections = FALSE)))
f12  <- quiet_fit(suppressMessages(prep_stan_list(d_sp$stim, found = d_sp$found, harmonics = c(1, 2),
                                                  save_selections = TRUE)))

# cat, von Mises vs Fourier-124
d_cat <- simulate(list(K = 1, nrow = 6, ncol = 6, prop_missing = 0, v_type = "cat", n_classes = 2), 12)
cat_vm <- quiet_fit(suppressMessages(prep_stan_list(d_cat$stim, found = d_cat$found, directional_model = "von_mises",
                                                    save_selections = TRUE)))
cat_f  <- quiet_fit(suppressMessages(prep_stan_list(d_cat$stim, found = d_cat$found, save_selections = TRUE)))

fits <- list(spatial = list(fourier124 = f124, fourier12 = f12),
             cat = list(vonmises = cat_vm, fourier124 = cat_f))
data <- list(spatial = d_sp, cat = d_cat)

failures <- character(0)
check <- function(label, ok) {
  cat(sprintf("%-55s %s\n", label, if (isTRUE(ok)) "OK" else "FAILED"))
  if (!isTRUE(ok)) failures <<- c(failures, label)
}

# one call per dataset, bound with a dataset column
compare_all <- function(...) {
  imap_dfr(fits, \(f, ds) suppressWarnings(suppressMessages(fomo_compare(f, data[[ds]], ...))) |>
             mutate(dataset = ds, .before = 1))
}
tab <- compare_all()
print(tab |> select(dataset, model, deterministic, nn, elpd_diff, se_diff, robust), width = 200)
check("no k-fold columns without kfold_cv", !any(c("check_kfold", "kfold_reason", "elpd_diff_kfold") %in% names(tab)))
check("one row per model", nrow(tab) == 4 && all(c("spatial", "cat") %in% tab$dataset))
check("accuracy for every row (incl. the fit without P/W)", all(!is.na(tab$deterministic)))
check("model beats the nearest-neighbour rule on spatial data",
      all(tab$gain_nn[tab$dataset == "spatial"] > 0))
check("Fourier-124 best on Fourier-124 data",
      tab$model[tab$dataset == "spatial" & tab$elpd_diff == 0] == "fourier124")
check("best model: elpd_diff 0 and robust",
      all(tab$robust[tab$elpd_diff == 0]) && sum(tab$elpd_diff == 0) == 2)
check("model weights sum to 1 per dataset",
      all(abs(summarise(tab, s = sum(stacking), .by = dataset)$s - 1) < 1e-6))

# the k-fold advice, forced on with extreme thresholds
msg <- character(0); warn <- character(0)
forced <- withCallingHandlers(fomo_compare(fits$spatial, data$spatial, kfold_se = Inf, kfold_high_k = -1),
                    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") },
                    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
check("k-fold advice fires on its thresholds (the non-best model)",
      nrow(attr(forced, "kfold")) == 1 && grepl("close", attr(forced, "kfold")$reason))
check("a warning for the flagged comparison", sum(grepl("fomo_compare\\(kfold_cv = K\\)", warn)) == 1)
check("a tibble, format_compare() a kable",
      inherits(tab, "tbl_df") && inherits(format_compare(tab, columns = c("dataset", "model", "elpd_diff"),
                                                         group = "dataset"), "knitr_kable"))

# the spatial dataset has under 100 trials: loo's "N < 100" is a caution
# (a message, not a column), never a reason for k-fold
check("no small_n column", !("small_n" %in% names(tab)))
check("N < 100 is not a k-fold reason", !any(grepl("N < 100", attr(forced, "kfold")$reason)))
check("a small-N message", any(grepl("fewer than 100 trials", msg)))

half <- d_sp$stim |> distinct(person, condition, trial) |> slice_head(n = 10)
f_half <- quiet_fit(suppressMessages(prep_stan_list(subset_trials(d_sp, half)$stim,
                                                    found = subset_trials(d_sp, half)$found)))
check("models fitted to different trials refused",
      inherits(try(suppressMessages(fomo_compare(list(a = f124, b = f_half), d_sp)), silent = TRUE), "try-error"))
check("several datasets in one call refused",
      inherits(try(suppressMessages(fomo_compare(fits, data)), silent = TRUE), "try-error"))
check("models = subset",
      nrow(suppressMessages(fomo_compare(fits$spatial, d_sp, models = "fourier124"))) == 1)

# k-fold on request: every model on the same folds, the fit's own settings
kf_warn <- character(0)
kf <- withCallingHandlers(suppressMessages(fomo_compare(fits$spatial, d_sp, kfold_cv = 2,
                                                        kfold_args = list(save_dir = tempdir()))),
                          warning = function(w) { kf_warn <<- c(kf_warn, conditionMessage(w)); invokeRestart("muffleWarning") })
check("kfold_cv adds elpd_diff_kfold and se_diff_kfold",
      all(c("elpd_diff_kfold", "se_diff_kfold") %in% names(kf)) && sum(kf$elpd_diff_kfold == 0) == 1)
check("kfold_cv: each fold fit's diagnostics", nrow(attr(kf, "kfold_fits")) == 4)
check("kfold_cv: no 'verify with k-fold' warnings", !any(grepl("kfold_cv = K", kf_warn)))

# clear errors
err <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
check("K = 1 refused", grepl("at least 2", err(kfold_fomo(d_sp, K = 1, sl = attr(f124, "fomo_sl")))))
check("group refused (not yet supported)", grepl("isn't supported", err(kfold_fomo(d_sp, K = 2, group = "person"))))
check("strata without person refused", grepl("strata must", err(kfold_fomo(d_sp, K = 2, strata = "condition"))))

if (length(failures) == 0) {
  cat("\nALL fomo_compare() TESTS PASSED\n")
} else {
  cat("\n", length(failures), "FAILURE(S):", paste(failures, collapse = ", "), "\n")
  quit(status = 1)
}
