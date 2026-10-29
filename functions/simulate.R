# =============================================================================
# SIMULATE
#
# Simulated data: stimuli (generate_stim()), selections from a simulate
# model (sim_fomo(), then extract_pred()) and early stopping
# (truncate_trials()). sim_fomo() with a fit instead simulates from its
# posterior (posterior predictive checks).
# =============================================================================

#' Generate experimental stimuli
#'
#' Simulates trial-level stimulus grids for one or more conditions and
#' participants. Grid dimensions, missingness, and feature type are controlled
#' via \code{spec}. Supports categorical (\code{"cat"}) and spatial-only
#' feature types.
#'
#' @param spec Named list describing the stimulus structure. Required fields:
#' \describe{
#'   \item{nrow}{Integer or vector. Number of grid rows per condition.}
#'   \item{ncol}{Integer or vector. Number of grid columns per condition.}
#'   \item{prop_missing}{Numeric or vector in \[0, 1\]. Proportion of grid
#'     cells to remove per condition.}
#'   \item{v_type}{Character. Feature type: \code{"cat"}, or anything else
#'     for spatial-only.}
#'   \item{K}{Integer. Number of conditions (optional; inferred from
#'     \code{ncol} if absent).}
#'   \item{n_classes}{Integer or vector. Number of item classes per condition
#'     (required for \code{v_type = "cat"}). If scalar, the same value is
#'     used for all conditions.}
#' }
#' @param n_trials Integer or vector. Number of trials per condition. If
#'   scalar, the same value is used for all conditions. Defaults to 1.
#' @param n_people Integer. Number of participants to simulate. Defaults to 1.
#' @param cond_labels Optional character vector of condition labels. If
#'   \code{NULL}, defaults to \code{c("one", "two", "three", "four")[1:K]}
#'   (or \code{"cond1"}, \code{"cond2"}, ... if K > 4).
#'   Must have length >= K.
#'
#' @return Tibble with one row per item per trial, containing columns \code{x},
#'   \code{y}, \code{id}, \code{trial}, \code{condition}, \code{person}, and
#'   \code{v} (feature values; absent for spatial-only).
#'
#' @seealso \code{\link{fit_fomo}}, \code{\link{sim_fomo}}
generate_stim <- function(spec, n_trials = 1, n_people = 1, cond_labels = NULL) {

  stopifnot("stimuli specification (spec) should be a list" = is.list(spec))
  v_type <- spec$v_type
  stopifnot("spec$v_type: one feature type" = length(v_type) == 1)

  K <- spec$K %||% length(spec$ncol)
  # one value per condition: a single value is used for every condition
  per_condition <- function(x) if (length(x) == 1) rep(x, K) else x
  ncol <- per_condition(spec$ncol)
  nrow <- per_condition(spec$nrow)
  prop_missing <- per_condition(spec$prop_missing)
  n_trials <- per_condition(n_trials)

  make_trial <- function(person, condition, trial) {

    k <- condition
    grid <- expand.grid(x = seq(0, 1, length.out = ncol[k]),
                        y = seq(0, 1, length.out = nrow[k])) |>
      as_tibble() |>
      sample_frac(1 - prop_missing[k]) |>
      mutate(id = row_number(), trial = trial, condition = k)

    if (v_type == "cat") {
      n_classes <- per_condition(spec$n_classes)[k]
      classes <- rep(seq_len(n_classes), each = ceiling(nrow(grid) / n_classes))
      # every trial's factor has the same levels, so combining trials can't
      # reorder (relabel) the classes when a trial misses one
      grid$v <- factor(sample(classes, nrow(grid)), levels = seq_len(max(spec$n_classes)))
    }

    mutate(grid, person = person)
  }

  # trials numbered 1, 2, ... over people, then conditions, then trials
  expand_grid(person = seq_len(n_people), condition = seq_len(K)) |>
    mutate(trial = map(condition, \(k) seq_len(n_trials[k]))) |>
    unnest(trial) |>
    mutate(trial = row_number()) |>
    pmap(make_trial) |>
    list_rbind() -> stim

  if (is.null(cond_labels)) {
    cond_labels <- if (K <= 4) c("one", "two", "three", "four")[1:K] else paste0("cond", 1:K)
  }
  stim$condition <- factor(stim$condition, labels = cond_labels)

  stim
}

sim_fomo <- function(sl, fit = NA, iter = 500, chains = 1) {

  mod <- load_stan_model(sl, model_task = "simulate")

  if (!inherits(fit, "CmdStanMCMC")) {
    # no fitted model provided: sample from priors. Quietly: no progress
    # messages, and no sampler diagnostics (E-BFMI etc. mean nothing here,
    # and warn when iter < 3)
    return(mod$sample(data = sl,
                      iter_warmup = 1000,
                      iter_sampling = iter,
                      chains = chains,
                      refresh = 0,
                      show_messages = FALSE,
                      diagnostics = NULL))
  }

  # fitted model provided: run generated quantities on a subsample of draws.
  # The fit must use the same directional model as sl, or the draws won't
  # line up with the simulate model's parameters
  fit_model <- if (any(str_detect(fit$metadata()$variables, "^fourier_coef"))) "fourier" else "von_mises"
  stopifnot("fitted model and sl use different directional models - build sl with the same directional_model as the fit" =
            fit_model == get_directional_model(sl))
  draws_matrix <- parameter_draws(fit, mod, sl)
  idx <- sample(nrow(draws_matrix), iter)

  run_quietly(\() mod$generate_quantities(fitted_params = draws_matrix[idx, ], data = sl))
}

#' Extract predicted sequences from a simulation fit
#'
#' Converts posterior draws of predicted item sequences (\code{Q}) into a tidy
#' tibble aligned with the observed data. Observations with \code{id == 0}
#' (padding) are removed. If only one draw is present, the \code{.draw} column
#' is dropped.
#'
#' @param m CmdStanMCMC or CmdStanGQ object containing \code{Q} draws,
#'   typically from \code{\link{sim_fomo}}.
#' @param ds Data frame of observed trial data. Must contain columns
#'   \code{trial}, \code{id}, \code{condition}, and \code{person}.
#'
#' @return Tibble with one row per predicted item per draw, containing columns
#'   from \code{ds} plus \code{n} (position within trial) and \code{.draw}
#'   (omitted if only one draw).
#'
#' @seealso \code{\link{sim_fomo}}
extract_pred <- function(sim, stim) {

  Q <- get_draws(sim, "Q")

  # number of items per trial (used to reconstruct trial indices)
  s <- count(stim, trial)$n

  Q |>
    pivot_longer(-.draw, names_to = "n", values_to = "id") |>
    mutate(n     = parse_number(n),
           trial = rep(rep(seq_along(s), s), max(Q$.draw))) |>
    left_join(stim, by = join_by(id, trial)) |>
    group_by(.draw, trial, condition) |>
    mutate(n = row_number()) |>
    filter(id != 0) |>
    ungroup() -> dsim

  if (n_distinct(dsim$.draw) == 1) dsim <- select(dsim, -.draw)

  dsim
}

#' Remove the last few selections of each trial
#'
#' Simulated sequences (extract_pred()) are exhaustive: every item is
#' selected. Real foragers often stop early, so this drops a random number
#' of selections from the end of each trial: Poisson(lambda) of them, but
#' always keeping at least min_keep (or the whole trial, if it is shorter).
#' Pass the result to prep_stan_list() as found with the same stim: the
#' removed items stay on the display, unselected.
#'
#' @param found Selections, e.g. from extract_pred(): columns trial and n
#'   (position in the trial, 1, 2, ...); with a .draw column, each draw's
#'   trials are truncated separately
#' @param lambda Mean number of selections removed per trial
#' @param min_keep Fewest selections a trial keeps
#'
#' @return found without the removed selections
truncate_trials <- function(found, lambda = 2, min_keep = 10) {

  found |>
    filter(n <= max(n() - rpois(1, lambda), min_keep),
           .by = c(any_of(".draw"), trial))
}
