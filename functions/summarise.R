# =============================================================================
# SUMMARISE
#
# Summaries of a fit: fomo_summary() (estimates, convergence, recovery),
# HPDIs, the directional and weight curves, loo diagnostics
# (check_high_pareto_k()), and the simulation values (ground truth) in the
# shapes the summaries and plots need.
# =============================================================================

# Population-level parameters Stan samples, for both directional models
# (excludes generated quantities and per-person random-effect draws). Used by
# fomo_variables() and analysis scripts that summarise fits.
FOMO_MAIN_PARAMS_REGEX <- paste0(
  "^(rho_delta|rho_psi|log_theta|fourier_coef|b_v|b_v_logodds|",
  "sigma_u_delta|sigma_u_psi|sigma_u_log_theta|sigma_u_fourier_coef|sigma_u_v)(\\[|$)")

# Every variable name in a fit. A cmdstanr fit's names come from its
# metadata: tidybayes::get_variables() loads every draw first, which ran
# out of memory on fits that store the per-selection P/W
fit_variables <- function(fit) {
  if (inherits(fit, "CmdStanFit")) fit$metadata()$variables else get_variables(fit)
}

# The parameters worth looking at in a fit: the fixed effects (the
# population-level parameters, FOMO_MAIN_PARAMS_REGEX, without the
# random-effect SDs) and, with random_effects = TRUE, the random effects
# too: their SDs (sigma_u_*) and each person's values (u_*). Used by
# fomo_summary() and plot_post_traces()
fomo_variables <- function(fit, random_effects = FALSE) {
  vars <- fit_variables(fit)
  fixed <- str_detect(vars, FOMO_MAIN_PARAMS_REGEX) & !str_starts(vars, "sigma_u_")
  random <- str_detect(vars, "^(sigma_u|u)_\\w+\\[")
  vars[fixed | (random_effects & random)]
}

#' Summary of a fit: estimates, convergence and, for simulated data, recovery
#'
#' Prints one table of the fixed effects (fomo_variables()): median, HPDI,
#' R-hat and bulk/tail ESS, and, with sim, the value the data were
#' simulated from (true) and whether it is inside the HPDI. Ends with the
#' recovery count and a convergence line: the worst R-hat and ESS over
#' every parameter, random effects included (so a problem there is still
#' named), and the number of divergent transitions.
#'
#' With random_effects = TRUE (multi-level fits), the table also has the
#' random-effect SDs (sigma_u_*), and each per-person parameter (u_*) gets
#' one row: how many values, how many were recovered, the worst R-hat and
#' lowest ESS. For every value of one parameter: variables = "u_psi".
#'
#' @param fit A fit (fit_fomo())
#' @param sim The simulation the data came from (sim_fomo(..., iter = 1)).
#'   It holds every true value, including the group-level SDs and each
#'   person's parameters (a simparams list doesn't)
#' @param random_effects Include the random effects in the tables?
#' @param variables Summarise just these parameters instead (e.g. "u_v"),
#'   one row per value
#' @param cond_labels Condition labels (e.g. levels(stim$condition));
#'   otherwise conditions are numbered
#' @param width HPDI width (default: the wide FOMO_HPDI_WIDTHS, 97%)
#' @return Invisibly, the summary: one row per parameter value
fomo_summary <- function(fit, sim = NULL, random_effects = FALSE, variables = NULL,
                         cond_labels = NULL, width = max(FOMO_HPDI_WIDTHS)) {

  summary <- summarise_fit(fit, sim, variables, cond_labels, width)
  labelled <- label_components(summary, fit)
  shown <- if (random_effects || !is.null(variables)) labelled else filter(labelled, !random_effect)
  has_truth <- !is.null(sim)
  interval <- paste0(100 * width, "% HPDI")
  estimates <- c("lower", "median", "upper", if (has_truth) c("true", "inside"),
                 "rhat", "ess_bulk", "ess_tail")

  print_table <- function(title, table) {
    cat("\n**", title, "**\n\n", sep = "")
    if ("component" %in% names(table)) {
      table <- if (all(is.na(table$component))) select(table, -component) else mutate(table, component = replace_na(component, ""))
    }
    if (knitr::is_html_output()) cat(styled_summary_table(table), "\n") else print(knitr::kable(table, digits = 2))
  }

  if (!is.null(variables)) {
    print_table(paste0(paste(variables, collapse = ", "), " (median and ", interval, ")"),
                select(labelled, parameter, condition, person, component, all_of(estimates)))
  } else {
    print_table(paste0(if (random_effects) "Population-level parameters" else "Fixed effects",
                       " (median and ", interval, ")"),
                select(filter(shown, !per_person), parameter, condition, component, all_of(estimates)))

    if (any(shown$per_person)) {
      shown |>
        filter(per_person) |>
        summarise(values = n(),
                  inside = if (has_truth) sum(inside, na.rm = TRUE) else NULL,
                  max_rhat = max(rhat, na.rm = TRUE),
                  min_ess_bulk = min(ess_bulk, na.rm = TRUE),
                  min_ess_tail = min(ess_tail, na.rm = TRUE),
                  .by = parameter) -> per_person
      print_table("Per-person parameters (one row per parameter)", per_person)
    }
  }

  if (has_truth) {
    n_inside <- sum(shown$inside, na.rm = TRUE)
    n_true <- sum(!is.na(shown$inside))
    cat(sprintf("\nRecovery: %d of %d (%.0f%%) simulation values inside the %s; about %g%% expected by chance.\n",
                n_inside, n_true, 100 * n_inside / n_true, interval, 100 * width))
  }

  # the parameter value with the worst R-hat or ESS (which.max/min skip
  # NAs), over every parameter, random effects included
  worst <- function(column, pick) summary$variable[pick(summary[[column]])]
  cat(sprintf("\nConvergence (all parameters): largest R-hat %.3f (%s); smallest bulk ESS %.0f (%s), tail ESS %.0f (%s); %d divergent transitions.\n",
              max(summary$rhat, na.rm = TRUE), worst("rhat", which.max),
              min(summary$ess_bulk, na.rm = TRUE), worst("ess_bulk", which.min),
              min(summary$ess_tail, na.rm = TRUE), worst("ess_tail", which.min),
              sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)))

  invisible(summary)
}

styled_summary_table <- function(table) {

  # fomo_summary()'s tables for rendered (HTML) pages: a vertical line after
  # each group of columns (parameter | condition, person, component |
  # lower, median, upper | values, true, inside | rhat | ess_bulk, ess_tail),
  # a horizontal line only between parameters, the inside cells shaded
  # pale green/orange and the R-hat and ESS cells by FOMO_RHAT_OK and
  # FOMO_ESS_OK
  column_groups <- list("parameter", c("condition", "person", "component"),
                        c("lower", "median", "upper"), c("values", "true", "inside"),
                        c("rhat", "max_rhat"))
  last_in_group <- map_int(column_groups, \(g) max(c(0L, match(g, names(table))), na.rm = TRUE))
  border_after <- setdiff(last_in_group, 0L)

  # first row of each parameter after the first
  new_parameter <- which(table$parameter != lag(table$parameter))

  # parameter and numbers right-aligned, other text centred
  align <- if_else(map_lgl(table, is.numeric) | names(table) == "parameter", "r", "c")
  styled <- kableExtra::kbl(table, digits = 2, align = align) |>
    kableExtra::kable_styling(bootstrap_options = "condensed", full_width = FALSE, position = "center",
                              font_size = 12) |>
    kableExtra::column_spec(border_after, border_right = "1px solid #999", include_thead = TRUE) |>
    kableExtra::row_spec(seq_len(nrow(table)), extra_css = "border-top: none; border-bottom: none;") |>
    kableExtra::row_spec(new_parameter, extra_css = "border-top: 1px solid #999;")

  # the parameter names and medians shaded light grey
  grey_cols <- match(intersect(c("parameter", "median"), names(table)), names(table))
  styled <- kableExtra::column_spec(styled, grey_cols, background = "#f2f2f2")

  if ("inside" %in% names(table) && is.logical(table$inside)) {
    styled <- styled |>
      kableExtra::column_spec(match("inside", names(table)),
                              background = case_when(table$inside ~ "#dff0d8",
                                                     !table$inside ~ "#fde2c4",
                                                     .default = "white"))
  }

  # good / acceptable / bad
  shade <- function(styled, col, good, bad) {
    if (!(col %in% names(table))) return(styled)
    kableExtra::column_spec(styled, match(col, names(table)),
                            background = case_when(good ~ "#dff0d8", bad ~ "#f8c9c4",
                                                   !is.na(good) ~ "#fcf0c4", .default = "white"))
  }
  for (col in c("rhat", "max_rhat")) {
    r <- table[[col]]
    styled <- shade(styled, col, r < FOMO_RHAT_OK[1], r >= FOMO_RHAT_OK[2])
  }
  for (col in c("ess_bulk", "ess_tail", "min_ess_bulk", "min_ess_tail")) {
    e <- table[[col]]
    styled <- shade(styled, col, e >= FOMO_ESS_OK[2], e < FOMO_ESS_OK[1])
  }
  styled
}

label_components <- function(summary, fit) {

  # summarise_fit()'s index as a component for the printed tables: Fourier
  # coefficients as h1 cos, h1 sin, ... (from the fit's harmonics); von
  # Mises components by preferred direction (before the condition's
  # offset, as direction_labels()); anything else as the index itself
  is_fourier <- str_detect(summary$parameter, "fourier_coef")
  is_theta <- str_detect(summary$parameter, "theta")
  j <- suppressWarnings(as.integer(summary$index))

  harmonics <- if (any(is_fourier)) as.integer(posterior::as_draws_matrix(fit$draws("fourier_harmonics"))[1, ])
  pair <- paste0("h", harmonics[(j + 1) %/% 2])
  n_dircomps <- if (any(is_theta)) max(j[is_theta], na.rm = TRUE) else 1

  summary |>
    mutate(component = case_when(
      is_fourier & !is.na(j) ~ paste(pair, if_else(j %% 2 == 1, "cos", "sin")),
      is_theta & !is.na(j) ~ direction_labels(n_dircomps)[pmax(j, 1)],
      .default = index))
}

summarise_fit <- function(fit, sim, variables, cond_labels, width) {

  # one row per parameter value of fit: median, HPDI, R-hat, ESS, and the
  # labels parsed from its name (parameter, condition, person, index); with
  # sim, the true value and whether it is inside the HPDI
  vars <- if (is.null(variables)) {
    fomo_variables(fit, random_effects = TRUE)
  } else {
    all_vars <- fit_variables(fit)
    all_vars[str_detect(all_vars, paste0("^(", paste(variables, collapse = "|"), ")(\\[|$)"))]
  }

  posterior::summarise_draws(fit$draws(variables = vars),
                             median = median,
                             \(x) hdci_bounds(x, width),
                             rhat = posterior::rhat,
                             ess_bulk = posterior::ess_bulk,
                             ess_tail = posterior::ess_tail) |>
    mutate(across(c(ess_bulk, ess_tail), round),
           parameter = str_remove(variable, "\\[.*"),
           per_person = str_starts(parameter, "u_"),
           random_effect = per_person | str_starts(parameter, "sigma_u_"),
           indices = str_split(str_extract(variable, "(?<=\\[)[0-9,]+"), ",")) -> summary

  # the first index is always the condition; per-person parameters also
  # carry a person index - last, except for the categorical u_v[k, l, c1, c]
  # and u_first_v[k, l, c], where it is second (and u_v's logits,
  # u_logit_v) - and any remaining indices
  # (direction/coefficient, dimension, class) go in `index`
  person_position <- map2_int(summary$parameter, summary$indices, function(p, i) {
    if (!str_starts(p, "u_")) return(NA_integer_)
    if (p == "u_first_v" || (p %in% c("u_v", "u_logit_v") && length(i) == 4)) return(2L)
    length(i)
  })

  summary |>
    mutate(condition = map_int(indices, \(i) as.integer(i[1])),
           person = map2_int(indices, person_position, \(i, pos) if (is.na(pos)) NA_integer_ else as.integer(i[pos])),
           index = map2_chr(indices, person_position, function(i, pos) {
             rest <- i[-c(1, if (!is.na(pos)) pos)]
             if (length(rest) == 0) NA_character_ else paste(rest, collapse = ",")
           })) |>
    select(-indices) -> summary

  if (!is.null(cond_labels)) summary$condition <- label_condition(summary$condition, cond_labels)

  if (!is.null(sim)) {
    stopifnot("sim should have a single draw - use sim_fomo(..., iter = 1)" =
                posterior::ndraws(sim$draws()) == 1)
    true_vars <- intersect(summary$variable, fit_variables(sim))
    truth_draw <- posterior::as_draws_matrix(sim$draws(variables = true_vars))
    truth <- setNames(as.numeric(truth_draw[1, ]), colnames(truth_draw))
    summary <- mutate(summary,
                      true = as.numeric(truth[variable]),
                      inside = true >= lower & true <= upper)
  }

  summary
}

# HPDI widths for every plot and summary (narrow, wide); the wide one is
# also fomo_summary()'s default
FOMO_HPDI_WIDTHS <- c(0.53, 0.97)

# Pareto k above this flags a trial whose loo estimate is unreliable
# (check_high_pareto_k(), fomo_compare())
FOMO_PARETO_K <- 0.7

# R-hat below the first is good, below the second acceptable, above it bad
# (fomo_summary()'s colours)
FOMO_RHAT_OK <- c(1.01, 1.05)

# ESS (bulk and tail) above the second is good, above the first acceptable,
# below it bad (Vehtari et al., 2021: 100 per chain, so 400 for 4 chains)
FOMO_ESS_OK <- c(100, 400)

#' Highest-density continuous interval of a vector of draws
#'
#' @param x Numeric vector of draws
#' @param width Interval width
#'
#' @return Named vector c(lower, upper)
hdci_bounds <- function(x, width) {
  h <- ggdist::hdci(as.numeric(x), .width = width)
  c(lower = h[1, 1], upper = h[1, 2])
}

compute_prior_hpdi <- function(draws, widths = FOMO_HPDI_WIDTHS) {

  # HPDIs of the prior draws (dist == "prior") in an extract_post() table,
  # per parameter and whatever else the plots facet by (each Fourier
  # harmonic and term, each von Mises direction, each class transition),
  # since their priors can differ (e.g. Fourier priors are per harmonic)
  draws |>
    filter(dist == "prior") |>
    group_by(across(any_of(c("parameter", "harmonic", "term", "k", "prev_class", "class")))) |>
    mean_hdci(value, .width = widths)
}

#' Population-level absolute-direction log-weight as a function of angle
#'
#' The term the directional model adds to each item's log-weight, for an
#' item at angle `angle` (radians, atan2(dy, dx), -pi to pi) from the last
#' selection, for each posterior draw and condition. Works for both
#' directional models, so their fits can be compared directly:
#'   von Mises: log(1 + sum_d theta_d * vonMises(angle; mu_d, kappa)),
#'     mu_d = 2 * pi * (d - 1) / n_dircomps + dir_offset, kappa = FOMO_KAPPA
#'   Fourier:   sum_h c_h * cos(h * a) + s_h * sin(h * a), a = angle - dir_offset
#' Also works on a simulation (extract_post() of sim_fomo(..., iter = 1)),
#' to get the true curve.
#'
#' @param post Output of extract_post() (or groundtruth_as_post())
#' @param n_angles Number of angles from -pi to pi
#' @param draws Which posterior draws to use (.draw values); NULL for all
#' @param centre Subtract each curve's mean over angle? A constant added
#'   to every item's log-weight doesn't change the selection
#'   probabilities, so only the shape of the curve matters; centring puts
#'   von Mises and Fourier curves on the same footing
#' @return Tibble: .draw, condition, angle, log_weight
compute_directional_curve <- function(post, n_angles = 181, draws = NULL, centre = TRUE) {

  angles <- seq(-pi, pi, length.out = n_angles)
  von_mises <- !is.null(post$vonmises)
  directional <- post$vonmises %||% post$fourier
  if (is.null(directional)) stop("post has no directional parameters (vonmises or fourier)")

  directional <- filter(directional, dist == "posterior")
  if (!is.null(draws)) directional <- filter(directional, .draw %in% draws)

  # one condition's curves: its (draws x components) weights times a
  # (components x angles) basis
  condition_curves <- function(d, key) {
    if (von_mises) {
      n_dir <- nlevels(d$k)
      d <- mutate(d, component = as.integer(k))
      preferred_direction <- 2 * pi * (seq_len(n_dir) - 1) / n_dir + d$offset[1]
      basis <- t(outer(angles, preferred_direction, vonmises_density, kappa = FOMO_KAPPA))
    } else {
      harmonics <- parse_number(levels(d$harmonic))
      # components in fourier_coef's order: {cos, sin} for each harmonic
      d <- mutate(d, component = 2 * as.integer(harmonic) - (term == "cos"))
      a <- angles - d$offset[1]
      basis <- do.call(rbind, lapply(harmonics, \(h) rbind(cos(h * a), sin(h * a))))
    }
    weights <- d |>
      select(.draw, component, value) |>
      pivot_wider(names_from = component, values_from = value) |>
      arrange(.draw)
    w <- as.matrix(weights[, as.character(seq_len(nrow(basis)))]) %*% basis
    if (von_mises) w <- log(1 + w)
    tibble(.draw = rep(weights$.draw, times = n_angles),
           angle = rep(angles, each = nrow(w)),
           log_weight = as.vector(w))
  }

  directional |>
    group_by(condition) |>
    group_modify(condition_curves) |>
    ungroup() -> curves

  if (centre) curves <- mutate(curves, log_weight = log_weight - mean(log_weight), .by = c(.draw, condition))
  select(curves, .draw, condition, angle, log_weight)
}

#' Weight curves for every term of the spatial model
#'
#' How much each term multiplies an item's selection weight, for each
#' posterior draw and condition (relative weights: the selection
#' probabilities are these, multiplied together and normalised):
#'   distance:  exp(-rho_delta * distance), 1 at distance 0. x is in
#'              display widths, up to half the display (the models see
#'              it times xy_scaling)
#'   turn:      exp(-rho_psi * psi), 1 straight ahead. x is the turn from
#'              the previous step's direction, in radians (psi = x / pi)
#'   direction: exp() of compute_directional_curve() (centred, so the
#'              geometric mean over directions is 1). x is in radians
#'
#' @param post Output of extract_post() (or groundtruth_as_post())
#' @param xy_scaling As in prep_stan_list() (default 10)
#' @param n_points Points per curve
#' @param n_draws Posterior draws to use (a random subset; enough for HPDIs)
#' @return Tibble: .draw, condition, term, x, weight
compute_weight_curves <- function(post, xy_scaling = 10, n_points = 101, n_draws = 500) {

  draws <- sample_draw_ids(post$spatial$.draw, n_draws)

  post$spatial |>
    filter(dist == "posterior", .draw %in% draws) |>
    select(.draw, condition, parameter, value) |>
    pivot_wider(names_from = parameter, values_from = value) -> rho

  distance <- rho |>
    cross_join(tibble(x = seq(0, 0.5, length.out = n_points))) |>
    mutate(term = "distance", weight = exp(-rho_delta * x * xy_scaling))

  turn <- rho |>
    cross_join(tibble(x = seq(0, pi, length.out = n_points))) |>
    mutate(term = "turn", weight = exp(-rho_psi * x / pi))

  direction <- if (!is.null(post$vonmises) || !is.null(post$fourier)) {
    compute_directional_curve(post, n_angles = n_points, draws = draws) |>
      mutate(term = "direction", x = angle, weight = exp(log_weight))
  }

  bind_rows(distance, turn, direction) |>
    select(.draw, condition, term, x, weight)
}

#' Trials with high Pareto k in a loo result
#'
#' log_lik is one value per trial, so loo leaves out a whole trial at a
#' time; a high Pareto k marks a trial the fit leans on heavily, whose loo
#' estimate is unreliable.
#'
#' @param loo_result Output of loo() on a fit's log_lik
#' @param stanlist The data list the fit used (prep_stan_list()), for each
#'   trial's person (z)
#' @param threshold Pareto k above which a trial is listed (default 0.67)
#' @return Tibble of those trials, worst first: trial_idx, z (person) and k
check_high_pareto_k <- function(loo_result, stanlist, threshold = FOMO_PARETO_K) {

  tibble(trial_idx = seq_along(loo_result$diagnostics$pareto_k),
         z = stanlist$z,
         k = loo_result$diagnostics$pareto_k) |>
    filter(trial_idx %in% pareto_k_ids(loo_result, threshold = threshold)) |>
    arrange(desc(k))
}

#' Parse ground truth parameters into plotting format
#'
#' Converts a named list of ground truth parameters (e.g. the \code{simparams}
#' used to simulate data in the simulated examples) into tidy tibbles for
#' \code{add_groundtruth()} / \code{plot_post()}. Handles the
#' spatial parameters (\code{rho_delta}, \code{rho_psi}) and, if present,
#' categorical transition
#' probabilities (\code{b_v_pr}) and von Mises directional weights
#' (\code{log_theta}). Fourier coefficients are not handled yet.
#'
#' @param gt Named list of ground truth parameters. Fields used:
#' \describe{
#'   \item{rho_delta_mu, rho_psi_mu}{Numeric vectors. One value per condition.}
#'   \item{b_v_pr}{Optional list of K n_classes x n_classes transition
#'     probability matrices.}
#'   \item{log_theta_mu}{Optional matrix of log-scale directional weights,
#'     one column per direction and one row per condition (or a single row,
#'     used for every condition).}
#' }
#' @param cond_labels Character vector of condition labels, used to convert
#'   integer condition indices to a factor.
#'
#' @return Named list with elements:
#' \describe{
#'   \item{spatial}{Tibble of spatial parameters: \code{parameter},
#'     \code{value}, and \code{condition} if there is more than one.}
#'   \item{vonmises}{Tibble of directional weights (\code{k},
#'     \code{condition}, \code{value}), or \code{NULL}.}
#'   \item{cat}{Tibble of transition probabilities (\code{condition},
#'     \code{class}, \code{prev_class}, \code{value}), or \code{NULL}.}
#' }
#'
#' @seealso \code{\link{add_groundtruth}}, \code{\link{sim_fomo}}
groundtruth_tables <- function(groundtruth, cond_labels) {

  gt_theta <- NULL
  gt_cat <- NULL

  as_tibble(groundtruth[c("rho_delta_mu", "rho_psi_mu")]) |>
    mutate(condition = row_number()) |>
    pivot_longer(-condition, names_to = "parameter") |>
    mutate(parameter = str_remove(parameter, "_mu")) -> gt_spatial

  # a condition column only if there is more than one condition
  if (any(gt_spatial$condition > 1)) {
    gt_spatial <- mutate(gt_spatial, condition = factor(condition, labels = cond_labels))
  } else {
    gt_spatial <- select(gt_spatial, -condition)
  }

  if ("b_v_pr" %in% names(groundtruth)) {
    # one matrix (used for every condition) or a list of one per condition
    if (is.matrix(groundtruth$b_v_pr)) groundtruth$b_v_pr <- rep(list(groundtruth$b_v_pr), length(cond_labels))
    n_classes <- ncol(groundtruth$b_v_pr[[1]])
    expand_grid(condition  = seq_along(groundtruth$b_v_pr),
                class      = seq_len(n_classes),
                prev_class = seq_len(n_classes)) |>
      mutate(value = pmap_dbl(list(condition, class, prev_class),
                              \(k, c, c1) groundtruth$b_v_pr[[k]][c, c1]),
             condition = factor(condition, labels = cond_labels),
             class = paste("class", class),
             prev_class = paste("prev. class", prev_class)) -> gt_cat
  }

  if ("log_theta_mu" %in% names(groundtruth)) {

    # a single row is used for every condition
    log_theta <- groundtruth$log_theta_mu
    if (nrow(log_theta) == 1) log_theta <- log_theta[rep(1, length(cond_labels)), , drop = FALSE]

    # the same facet labels as extract_vonmises() (direction_labels()), so
    # groundtruth vlines land on the right facet
    log_theta |>
      as_tibble(.name_repair = \(nms) paste0("k", seq_along(nms))) |>
      mutate(condition = row_number()) |>
      pivot_longer(-condition, names_to = "k") |>
      mutate(k = factor(parse_number(k), labels = direction_labels(ncol(log_theta))),
             condition = factor(condition, labels = cond_labels)) -> gt_theta
  }

  list(spatial = gt_spatial, vonmises = gt_theta, cat = gt_cat)
}

groundtruth_from_sim <- function(simparams, post) {

  # the values a multi-level simulation used: simparams plus each person's
  # parameters, read from the simulation's (single) draw
  groundtruth <- simparams
  groundtruth$u_spatial <- select(post$u$spatial, -.draw)
  if (!is.null(post$u$vonmises)) groundtruth$u_vonmises <- select(post$u$vonmises, -.draw)
  if (!is.null(post$u$cat)) groundtruth$u_cat <- select(post$u$cat, -.draw)

  groundtruth
}

#' Simulation values in the format of extract_post()
#'
#' The values data were simulated with (the simulated examples' simparams), as a single
#' "draw" of the spatial and directional parameters laid out like
#' extract_post()'s output, so functions written for posteriors (e.g.
#' compute_weight_curves()) give the true curves too. Directional
#' parameters come from log_theta_mu (von Mises) or fourier_coef_mu
#' (Fourier), whichever post has; a single row is used for every
#' condition. Condition labels, direction labels, harmonics and offsets
#' are taken from post.
#'
#' @param groundtruth Named list of simulation values (rho_delta_mu,
#'   rho_psi_mu, and log_theta_mu or fourier_coef_mu)
#' @param post Output of extract_post() for a fit to the simulated data
#' @return List with spatial and vonmises or fourier, as in extract_post()
groundtruth_as_post <- function(groundtruth, post) {

  cond_labels <- levels(post$spatial$condition)
  n_cond <- length(cond_labels)
  as_one_draw <- function(d) mutate(d, .draw = 1L, dist = "posterior")
  offsets <- distinct(bind_rows(post$vonmises, post$fourier), condition, offset)

  # values are expanded to one per condition as set_sim_params() does (a
  # single value or row is used for every condition)
  tibble(condition = label_condition(seq_len(n_cond), cond_labels),
         rho_delta = expand_vec(groundtruth$rho_delta_mu, n_cond, "rho_delta_mu"),
         rho_psi = expand_vec(groundtruth$rho_psi_mu, n_cond, "rho_psi_mu")) |>
    pivot_longer(c(rho_delta, rho_psi), names_to = "parameter") |>
    as_one_draw() -> spatial

  truth <- list(spatial = spatial)

  if (!is.null(post$vonmises) && !is.null(groundtruth$log_theta_mu)) {
    directions <- levels(post$vonmises$k)
    theta <- exp(expand_mat(groundtruth$log_theta_mu, n_cond, length(directions), "log_theta_mu"))
    truth$vonmises <- expand_grid(condition = seq_len(n_cond), dir = seq_along(directions)) |>
      mutate(value = theta[cbind(condition, dir)],
             k = factor(directions[dir], levels = directions),
             condition = label_condition(condition, cond_labels)) |>
      left_join(offsets, by = "condition") |>
      select(condition, k, value, offset) |>
      as_one_draw()
  }

  if (!is.null(post$fourier) && !is.null(groundtruth$fourier_coef_mu)) {
    harmonics <- parse_number(levels(post$fourier$harmonic))
    coef <- expand_mat(groundtruth$fourier_coef_mu, n_cond, 2 * length(harmonics), "fourier_coef_mu")
    truth$fourier <- expand_grid(condition = seq_len(n_cond), j = seq_len(ncol(coef))) |>
      mutate(value = coef[cbind(condition, j)]) |>
      label_fourier_coef(harmonics) |>
      mutate(condition = label_condition(condition, cond_labels)) |>
      left_join(offsets, by = "condition") |>
      select(condition, harmonic, term, value, offset) |>
      as_one_draw()
  }

  truth
}

# loo::loo_compare() as a tibble with model, elpd_diff, se_diff and
# diag_diff ("" where loo has no warning), whichever loo version: >= 2.10
# returns a data frame with a model column, earlier ones a matrix with the
# models as row names
fomo_loo_compare <- function(loos) {
  cmp <- as.data.frame(loo::loo_compare(loos))
  if (!"model" %in% names(cmp)) cmp <- rownames_to_column(cmp, "model")
  if (!"diag_diff" %in% names(cmp)) cmp$diag_diff <- ""
  as_tibble(cmp) |> mutate(diag_diff = replace_na(as.character(diag_diff), ""))
}

#' The same models compared by several methods (e.g. loo and k-fold), as
#' one table: method, model, elpd_diff, se_diff, best model first within
#' each method. For a styled table, pipe into format_compare(group =
#' "method").
#'
#' @param ... Named by method, each a named list (by model) of loo or
#'   k-fold elpd objects (kfold_fomo()$elpd) for the same trials
compute_elpd_comparison <- function(...) {
  methods <- list(...)
  stopifnot("name each method, e.g. loo = list(...)" = !is.null(names(methods)) && all(names(methods) != ""))
  map(methods, fomo_loo_compare) |>
    bind_rows(.id = "method") |>
    select(method, model, elpd_diff, se_diff)
}

#' Compare models fitted to one dataset
#'
#' A tibble with one row per model (format_compare() makes a table of it):
#' accuracy (deterministic, stochastic, chance, and the nearest-neighbour
#' rule, `nn`, with the model's gain over it), loo (elpd, the difference
#' from the best model and its SE, the share of trials with Pareto k above
#' FOMO_PARETO_K), model weights (stacking and pseudo-BMA+) and the
#' robustness check (the difference recomputed without the trials with
#' high k in either model). The models must be fitted to the same trials.
#' For several datasets, call it once per dataset (e.g. with imap()) and
#' bind the results.
#'
#' K-fold cross-validation (kfold_fomo()) costs K fits per model, so it
#' only runs when asked for (kfold_cv = K): then every model is refitted on
#' the same K folds, with the settings of its own fit (its data list and
#' sampler settings), and elpd_diff_kfold / se_diff_kfold are added
#' (attr(result, "kfold_fits"): each fold fit's diagnostics). Otherwise,
#' comparisons that should be checked with k-fold get a warning each, and
#' are listed in attr(result, "kfold") (model, reason): those that aren't
#' robust (without the high-k trials the best model no longer leads by more
#' than 2 SE), are close (under kfold_se SE) while more than kfold_high_k of
#' the trials have high k in either model, or that loo_compare() itself
#' warns about (diag_diff, apart from its "|elpd_diff| < 4", which is a
#' tie, and "N < 100", which k-fold, on the same trials, can't fix). A model
#' within 2 SE of the best is a tie (a message says so). With fewer than
#' 100 trials a message says to treat the SEs as approximate.
#'
#' @param fits Named list of fits by model, e.g. list(fourier124 = m1,
#'   vonmises = m2)
#' @param data The data they were fitted to: list(stim, found)
#' @param kfold_cv NULL (default): loo only; K: also K-fold
#'   cross-validation for every model
#' @param models Optional: compare only these models
#' @param sls Optional data lists by model, for fits made before fit_fomo()
#'   kept them (attr(fit, "fomo_sl")): accuracy from fits without stored
#'   P/W, and k-fold
#' @param loos Optional loo objects by model (by default fit$loo())
#' @param accuracy_exclude Which selections accuracy leaves out
#'   (keep_for_accuracy())
#' @param n_draws Draws for predict_fomo()
#' @param kfold_se, kfold_high_k Thresholds of the "close" k-fold rule
#' @param kfold_args Other kfold_fomo() arguments (e.g. strata, save_dir)
#' @param seed For the pseudo-BMA+ Bayesian bootstrap and predict_fomo()
fomo_compare <- function(fits, data, kfold_cv = NULL, models = NULL, sls = NULL, loos = NULL,
                         accuracy_exclude = "first_and_forced", n_draws = 500,
                         kfold_se = 4, kfold_high_k = 0.05, kfold_args = list(), seed = FOMO_SEED) {

  stopifnot("fits: a named list of fits to one dataset (for several datasets, call fomo_compare() for each)" =
              is.list(fits) && !is.null(names(fits)) && all(names(fits) != "") &&
              all(map_lgl(fits, \(x) inherits(x, "CmdStanFit"))),
            "data: list(stim, found)" = all(c("stim", "found") %in% names(data)))
  if (!is.null(models)) fits <- fits[intersect(names(fits), models)]
  sl_of <- function(m) sls[[m]] %||% attr(fits[[m]], "fomo_sl")

  # accuracy, and the nearest-neighbour rule on the same selections
  nn <- mean(compute_nn_baseline(data$stim, data$found, accuracy_exclude)$correct)
  acc <- imap_dfr(fits, function(fit, m) {
    stored <- prod(fit$metadata()$stan_variable_sizes$P %||% 0) > 0
    if (!stored && is.null(sl_of(m))) {
      message(m, ": no stored P/W and no data list, so no accuracy")
      return(tibble(model = m, deterministic = NA_real_, stochastic = NA_real_, chance = NA_real_))
    }
    compute_selection_accuracy(fit, data, sl_of(m), accuracy_exclude, n_draws, seed) |>
      summarise(deterministic = mean(deterministic), stochastic = mean(stochastic),
                chance = mean(chance)) |>
      mutate(model = m, .before = 1)
  })

  # loo
  lo <- imap(fits, \(fit, m) loos[[m]] %||% fit$loo())
  stopifnot("models fitted to different trials: compare only fits to the same data" =
              n_distinct(map_int(lo, \(x) nrow(x$pointwise))) == 1)
  out <- acc |>
    mutate(nn = nn, gain_nn = deterministic - nn) |>
    left_join(compare_loos(lo, kfold_se, kfold_high_k, seed), by = "model")

  # one model: no comparison columns
  comparison <- list(elpd_diff = NA_real_, se_diff = NA_real_, stacking = NA_real_,
                     pseudo_bma = NA_real_, elpd_diff_ok = NA_real_, se_diff_ok = NA_real_,
                     tie = NA, robust = NA, check_kfold = NA, kfold_reason = NA_character_)
  out[setdiff(names(comparison), names(out))] <- comparison[setdiff(names(comparison), names(out))]

  if (any(out$small_n)) message("fewer than 100 trials, so treat the SEs as approximate (for k-fold too)")
  ties <- filter(out, tie %in% TRUE)
  pwalk(list(ties$model, ties$elpd_diff, ties$se_diff), \(m, d, se)
        message(m, ": a tie with the best model on loo (", round(d, 1), " +- ", round(se, 1),
                ", within 2 SE): prefer the simpler model"))

  flagged <- filter(out, check_kfold %in% TRUE) |> select(model, reason = kfold_reason)
  kfold_fits <- NULL
  if (is.null(kfold_cv)) {
    # recommend k-fold where loo is doubtful, never run it unasked
    pwalk(flagged, \(model, reason)
          warning(model, ": verify this comparison with fomo_compare(kfold_cv = K) (", reason, ")",
                  call. = FALSE))
  } else {
    # every model on the same folds, with its own fit's settings
    missing_sl <- names(fits)[map_lgl(names(fits), \(m) is.null(sl_of(m)))]
    if (length(missing_sl) > 0) {
      stop("k-fold needs each fit's data list (fits from fit_fomo() keep it; ",
           "otherwise pass sls): missing for ", paste(missing_sl, collapse = ", "), call. = FALSE)
    }
    message("k-fold: ", length(fits) * kfold_cv, " fits (", length(fits), " models x K = ", kfold_cv, ")")
    kf <- imap(fits, function(fit, m) {
      md <- fit$metadata()
      do.call(kfold_fomo, c(list(data, K = kfold_cv, sl = sl_of(m), quiet = TRUE,
                                 filename = paste0("kfold_", make.names(m)),
                                 iter = md$iter_sampling, chains = fit$num_chains(),
                                 threads_per_chain = md$threads_per_chain %||% 2, refresh = 0),
                            kfold_args))
    })
    out <- left_join(out,
                     fomo_loo_compare(map(kf, "elpd")) |>
                       select(model, elpd_diff_kfold = elpd_diff, se_diff_kfold = se_diff),
                     by = "model")
    kfold_fits <- map_dfr(kf, "fits", .id = "model")
  }

  out <- out |>
    select(-check_kfold, -kfold_reason, -small_n) |>
    relocate(model, any_of(c("deterministic", "stochastic", "chance")), nn, gain_nn,
             elpd, elpd_diff, se_diff, any_of(c("elpd_diff_kfold", "se_diff_kfold")),
             pct_high_k, stacking, pseudo_bma, elpd_diff_ok, se_diff_ok, tie, robust)
  attr(out, "kfold") <- flagged
  attr(out, "kfold_fits") <- kfold_fits
  out
}

# fomo_compare()'s loo columns for one dataset's loo objects: elpd and the
# share of high-k trials for every model; with two or more models, the
# difference from the best, the weights, the robustness check and the
# k-fold flag
compare_loos <- function(lo, kfold_se = 4, kfold_high_k = 0.05, seed = FOMO_SEED) {

  high_k <- map(lo, \(x) x$diagnostics$pareto_k > FOMO_PARETO_K)
  out <- tibble(model = names(lo),
                elpd = map_dbl(lo, \(x) x$estimates["elpd_loo", "Estimate"]),
                pct_high_k = 100 * map_dbl(high_k, mean),
                # loo's "N < 100": too few trials for the SE's normal
                # approximation (a caution; k-fold can't fix it)
                small_n = map_int(lo, \(x) nrow(x$pointwise)) < 100)
  if (length(lo) < 2) return(out)

  cmp <- fomo_loo_compare(lo)
  best <- cmp$model[1]
  stacking <- loo::loo_model_weights(lo, method = "stacking")
  # seeded locally (see predict_fomo())
  pseudo_bma <- withr::with_seed(seed, loo::loo_model_weights(lo, method = "pseudobma", BB = TRUE))

  # the robustness check: the difference without the high-k trials
  ok <- map_dfr(names(lo), function(m) {
    bad <- high_k[[m]] | high_k[[best]]
    dd <- (lo[[m]]$pointwise[, "elpd_loo"] - lo[[best]]$pointwise[, "elpd_loo"])[!bad]
    tibble(model = m, elpd_diff_ok = sum(dd),
           se_diff_ok = if (m == best) 0 else sqrt(length(dd)) * sd(dd),
           pct_high_k_either = 100 * mean(bad))
  })

  out |>
    left_join(select(cmp, model, elpd_diff, se_diff, diag_diff), by = "model") |>
    left_join(ok, by = "model") |>
    mutate(stacking = as.numeric(stacking[model]),
           pseudo_bma = as.numeric(pseudo_bma[model]),
           # a tie: within 2 SE of the best on loo. Nothing to verify: k-fold
           # wouldn't separate them either, so no robustness check or flag
           tie = model != best & abs(elpd_diff) < 2 * se_diff,
           robust = if_else(tie, NA, model == best | elpd_diff_ok < -2 * se_diff_ok),
           close = !tie & model != best & abs(elpd_diff) < kfold_se * se_diff &
             pct_high_k_either > 100 * kfold_high_k,
           # loo's own warnings, except "|elpd_diff| < 4", which only says the
           # difference is small (a tie, above), and "N < 100", which says
           # the SE is approximate: k-fold, on the same trials, wouldn't help
           diag_diff = str_remove_all(diag_diff, "\\|elpd_diff\\| < 4|N < 100") |>
             str_replace_all("[;,] *[;,]", ";") |> str_remove_all("^[;, ]+|[;, ]+$"),
           loo_warns = model != best & diag_diff != "",
           check_kfold = !tie & (!(robust %in% TRUE) | close | loo_warns),
           kfold_reason = if_else(tie, "tie: within 2 SE",
                                  str_c(if_else(robust %in% FALSE, "not robust without high-k trials; ", ""),
                                        if_else(close, "close, with many high-k trials; ", ""),
                                        if_else(loo_warns, str_c("loo warns: ", diag_diff), "")) |>
                                    str_remove("; $"))) |>
    select(-close, -loo_warns, -diag_diff, -pct_high_k_either) |>
    arrange(desc(elpd_diff))
}

# the loo difference between two models fitted to the same trials, summed
# over each condition's trials: elpd_diff (first model - second) and its SE.
# loos: a named list of two loo objects; stim: the data they were fitted to
# (the data list orders trials as stim's trial numbers)
compute_elpd_by_condition <- function(loos, stim) {
  stopifnot("loos: a named list of two loo objects" = length(loos) == 2 && !is.null(names(loos)))
  trial_condition <- stim |> distinct(trial, condition) |> arrange(trial)
  stopifnot("loos and stim have different numbers of trials" =
              all(map_int(loos, \(x) nrow(x$pointwise)) == nrow(trial_condition)))
  trial_condition |>
    mutate(diff = loos[[1]]$pointwise[, "elpd_loo"] - loos[[2]]$pointwise[, "elpd_loo"]) |>
    summarise(elpd_diff = sum(diff), se = sqrt(n()) * sd(diff), .by = condition)
}

# fomo_compare() (or compute_elpd_comparison()) output as a rounded table:
# on rendered (HTML) pages unshaded, with a line under the header and,
# with group, between groups (rows with a new value of that column: e.g.
# "dataset" after binding several fomo_compare() results, "method" for
# compute_elpd_comparison()); elsewhere a plain knitr::kable(). columns:
# which to show (default all)
format_compare <- function(out, columns = NULL, group = NULL) {
  if (!is.null(columns)) out <- select(out, all_of(columns))
  out <- out |>
    mutate(across(any_of(c("deterministic", "stochastic", "chance", "nn", "gain_nn",
                           "stacking", "pseudo_bma")), ~ round(.x, 3)),
           across(any_of(c("elpd", "elpd_diff", "se_diff", "elpd_diff_kfold", "se_diff_kfold",
                           "elpd_diff_ok", "se_diff_ok", "pct_high_k")), ~ round(.x, 1)))
  if (!knitr::is_html_output()) return(knitr::kable(out))

  # first row of each group after the first
  new_group <- if (!is.null(group) && group %in% names(out)) which(out[[group]] != lag(out[[group]])) else integer(0)
  kableExtra::kbl(out) |>
    kableExtra::kable_styling(bootstrap_options = "condensed", full_width = FALSE, position = "center",
                              font_size = 12) |>
    kableExtra::row_spec(seq_len(nrow(out)), extra_css = "border-top: none; border-bottom: none;") |>
    kableExtra::row_spec(new_group, extra_css = "border-top: 1px solid #999;") |>
    kableExtra::row_spec(0, extra_css = "border-bottom: 1px solid #999;")
}
