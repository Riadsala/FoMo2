##########################################################################
# kfold_fomo()'s folds, with stand-in fit_fomo() and predict_fomo() (no
# Stan, a few seconds). The stand-in prediction gives each held-out trial
# log_lik = 100 * person + condition, so the returned elpd shows whether
# every trial was held out once and predicted as the right person and
# condition. On kristjansson2014plos (16 people x 2 conditions x 20 trials):
#   1. full data, K = 10
#   2. uneven cells (one down to 2 trials), K = 10: allowed, with a message
#   3. people 1 and 16 (first and last) down to 5 and 3 trials per
#      condition, K = 20: held-out folds that miss a person keep the other
#      people's numbers (subset_trials(renumber_people = FALSE),
#      prep_stan_list_like(all_people = FALSE))
#   4. refused: a 1-trial cell, K larger than any cell, people not 1..L,
#      a held-out fold with only 1 person
# The real fits and predictions are tested in test_kfold_fomo.R.
# Exits 1 on any failure.
#   Rscript tests/test_kfold_folds.R
##########################################################################

root <- rprojroot::find_root(rprojroot::has_dir("models"))
setwd(file.path(root, "tests"))
suppressMessages(source("../functions/set_up_FoMo.R"))

failures <- character(0)
report <- function(label, expr) {
  res <- tryCatch(expr, error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(sprintf("%-34s %s\n", label, res))
  if (startsWith(res, "FAILED")) failures <<- c(failures, label)
}

# stand-ins, found by kfold_fomo() in the global environment: the fit
# carries its data list; the prediction checks it matches the fit's people
# and conditions, and returns 4 identical draws per held-out trial
fit_fomo <- function(sl, ...) {
  list(sl = sl, diagnostic_summary = \(...) list(num_divergent = 0),
       summary = \(...) tibble(rhat = 1))
}
fomo_variables <- function(fit, ...) "stand-in"
predict_fomo <- function(fit, sl, ...) {
  stopifnot("held-out data list has other people or conditions than the fit" =
              sl$L == fit$sl$L && sl$K == fit$sl$K)
  m <- matrix(rep(100 * sl$z + sl$X, each = 4), nrow = 4)
  list(draws = \(...) m)
}

trials_of <- function(d) d$stim |> distinct(person, condition, trial)
# d's stim and found (d may hold other entries, e.g. dataset)
on_tables <- function(d, f) {
  d[c("stim", "found")] <- map(d[c("stim", "found")], f)
  d
}
keep_trials <- function(d, trials) {
  on_tables(d, \(x) semi_join(x, trials, by = c("person", "condition", "trial")))
}
# every trial held out exactly once, with its own person and condition;
# kfold_fomo()'s messages kept as attr(, "messages")
held_out_right <- function(d, K) {
  msgs <- character(0)
  out <- withCallingHandlers(kfold_fomo(d, K = K),
                             message = \(m) {
                               msgs <<- c(msgs, conditionMessage(m))
                               invokeRestart("muffleMessage")
                             })
  attr(out, "messages") <- msgs
  truth <- trials_of(d) |> arrange(trial) |> mutate(code = 100 * person + as.numeric(condition))
  stopifnot("not one elpd per trial" = nrow(out$elpd$pointwise) == nrow(truth),
            "a held-out trial predicted as the wrong person or condition" =
              all(out$elpd$pointwise[, "elpd"] == truth$code),
            "a trial in no fold or in two" = nrow(out$folds) == nrow(truth))
  out
}
refused <- function(expr, pattern) {
  msg <- tryCatch({ expr; "" }, error = \(e) conditionMessage(e))
  grepl(pattern, msg, fixed = TRUE)
}

d <- suppressWarnings(suppressMessages(import_data("kristjansson2014plos")))
first_n <- function(d, keep) {
  # keep(person, condition, i): whether the i-th trial of a cell is kept
  trials_of(d) |> arrange(person, condition, trial) |>
    filter(keep(person, as.numeric(condition), row_number()), .by = c(person, condition))
}

report("1. full data, K = 10", {
  out <- held_out_right(d, 10)
  sprintf("%d trials, %d per fold", nrow(out$folds), unique(count(out$folds, fold)$n))
})

report("2. uneven cells, K = 10", {
  # person 5's condition 1 down to 2 trials, 2 and 3 trials off two others
  d2 <- keep_trials(d, first_n(d, \(p, k, i) !(p == 5 & k == 1 & i > 2) &
                                             !(p == 1 & k == 1 & i <= 2) &
                                             !(p == 3 & k == 2 & i <= 3)))
  msgs <- attr(held_out_right(d2, 10), "messages")
  stopifnot("no message about the short cell" =
              any(grepl("1 of 32 person x condition cells have fewer than K = 10", msgs, fixed = TRUE)))
  sprintf("%d trials, short cell reported", nrow(trials_of(d2)))
})

report("3. people missing from folds", {
  d3 <- keep_trials(d, first_n(d, \(p, k, i) !(p == 1 & i > 5) & !(p == 16 & i > 3)))
  out <- held_out_right(d3, 20)
  missing_from <- \(p) sum(!(1:20 %in% out$folds$fold[out$folds$person == p]))
  stopifnot("person 1 should be missing from 15 folds" = missing_from(1) == 15,
            "person 16 should be missing from 17 folds" = missing_from(16) == 17)
  sprintf("%d trials; person 1 missing from %d folds, person 16 from %d",
          nrow(trials_of(d3)), missing_from(1), missing_from(16))
})

report("4. refusals", {
  one <- keep_trials(d, first_n(d, \(p, k, i) !(p == 5 & k == 1 & i > 1)))
  gappy <- on_tables(d, \(x) mutate(x, person = person * 2))
  # people 2-16 down to 2 trials per cell: folds 3-10 hold out only person 1
  lonely <- keep_trials(d, first_n(d, \(p, k, i) p == 1 | i <= 2))
  stopifnot("a 1-trial cell not refused" = refused(kfold_fomo(one, K = 3, quiet = TRUE), "only 1 trial"),
            "K larger than any cell not refused" = refused(kfold_fomo(d, K = 21, quiet = TRUE), "would be empty"),
            "people not 1..L not refused" = refused(kfold_fomo(gappy, K = 3, quiet = TRUE), "numbered 1..16"),
            "a 1-person held-out fold not refused" =
              refused(kfold_fomo(lonely, K = 10, quiet = TRUE), "trials of only 1 person"))
  "1-trial cell, K too large, people not 1..L, 1-person fold"
})

if (length(failures) == 0) {
  cat("\nALL kfold_fomo() FOLD TESTS PASSED\n")
} else {
  cat("\n", length(failures), "FAILURE(S):", paste(failures, collapse = ", "), "\n")
  quit(status = 1)
}
