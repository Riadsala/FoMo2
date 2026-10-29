# =============================================================================
# EXTRACT
#
# Draws out of a fit or simulation as tidy tables. extract_post() is the one
# to call: it returns every part the model has (spatial, cat, first_v,
# vonmises, fourier, per-person values u, random-effect SDs, accuracy);
# the extract_*() functions below build its parts.
# =============================================================================

extract_post <- function(fit, found, stim = NULL, accuracy_exclude = "first_and_forced") {

  # fit: a fit, a simulation, or predict_fomo() output. stim (optional) gives
  # each trial's number of items, so accuracy can leave out forced last
  # picks; accuracy_exclude: which selections accuracy leaves out (see
  # keep_for_accuracy())
  cond_labels <- levels(found$condition)
  vars <- fit_variables(fit)

  random_effects <- extract_random_effects(fit, vars, cond_labels)

  list(spatial  = extract_spatial(fit, vars, cond_labels),
       cat      = extract_cat(fit, vars, cond_labels),
       first_v  = extract_first_v(fit, vars, cond_labels),
       vonmises = extract_vonmises(fit, vars, cond_labels),
       fourier  = extract_fourier(fit, vars, cond_labels),
       acc      = accuracy_by_draw(fit, found, vars, stim, accuracy_exclude),
       u        = purrr::compact(random_effects),
       sigmas   = random_effects$sigmas) |>
    # drop the parts this model doesn't have
    purrr::discard(\(x) is.null(x) || (is.list(x) && length(x) == 0))
}

#' Extract posterior draws as tibble
#'
#' Converts CmdStanR draws into a tidy tibble and removes chain and iteration
#' indices (retaining `.draw` as the primary index).
#'
#' @param f CmdStanMCMC object
#' @param v Character vector of variable names
#'
#' @return Tibble of posterior draws
get_draws <- function(fit, variables) {

  fit$draws(variables = variables, format = "df") |>
    as_tibble() |>
    select(-.chain, -.iteration)
}

sample_draw_ids <- function(draw_ids, n) {
  # up to n of the draw numbers (.draw) at random, or all if there are fewer
  draw_ids <- unique(draw_ids)
  draw_ids[sample.int(length(draw_ids), min(n, length(draw_ids)))]
}

tidy_fomo_draws <- function(fit, vars, index_names) {
  # draws of vars in long format, one row per draw and variable: dist
  # ("prior" for the pr_ variables), parameter (the name without pr_ or
  # indices) and one integer column per Stan index, named index_names
  index_regex <- paste0("\\[", paste(rep("(\\d+)", length(index_names)), collapse = ","), "\\]")

  get_draws(fit, vars) |>
    pivot_longer(-.draw) |>
    mutate(dist = if_else(str_starts(name, "pr_"), "prior", "posterior"),
           parameter = name |> str_remove("^pr_") |> str_remove("\\[.*")) |>
    tidyr::extract(name, index_names, index_regex, convert = TRUE)
}

label_condition <- function(condition, cond_labels) {
  # condition index (1..K) -> factor labelled in condition order
  factor(condition, levels = seq_along(cond_labels), labels = cond_labels)
}

extract_spatial <- function(fit, vars, cond_labels) {

  vars_c <- vars[str_detect(vars, "rho_delta|rho_psi")]

  tidy_fomo_draws(fit, vars_c, "condition") |>
    mutate(condition = label_condition(condition, cond_labels),
           dist = factor(dist)) |>
    select(.draw, value, condition, parameter, dist)
}

extract_cat <- function(fit, vars, cond_labels) {

  # class-transition probabilities: b_v[k, c1] (and the prior draws,
  # pr_b_v) is a simplex over the next class, at both levels
  vars_c <- vars[str_detect(vars, "^(pr_)?b_v\\[")]
  if (length(vars_c) == 0) return(NULL)

  tidy_fomo_draws(fit, vars_c, c("condition", "prev_class", "class")) |>
    mutate(condition = label_condition(condition, cond_labels),
           prev_class = factor(prev_class),
           class = factor(class)) |>
    select(.draw, parameter, condition, prev_class, dist, class, value) |>
    arrange(.draw, dist, condition, prev_class, class)
}

extract_first_v <- function(fit, vars, cond_labels) {

  # the first selection's class weights (see first_item_weight()): the
  # stationary distribution of b_v, derived, on the probability scale.
  # Single-level only; multi-level models have them per person (u_first_v)
  vars_i <- vars[str_detect(vars, "^first_v\\[")]
  if (length(vars_i) == 0) return(NULL)

  tidy_fomo_draws(fit, vars_i, c("condition", "class")) |>
    mutate(condition = label_condition(condition, cond_labels),
           class = factor(class)) |>
    select(.draw, parameter, condition, class, dist, value) |>
    arrange(.draw, condition, class)
}

get_dir_offset <- function(fit, vars) {
  # directional offset for each condition (the same in every draw; models
  # without a directional_offset output have none)
  vars_o <- vars[str_detect(vars, "^directional_offset")]
  if (length(vars_o) == 0) return(0)
  as.numeric(get_draws(fit, vars_o)[1, vars_o])
}

extract_vonmises <- function(fit, vars, cond_labels) {

  vars_t <- vars[str_detect(vars, "^(pr_)*theta")]
  # Fourier fits have no theta
  if (length(vars_t) == 0) return(NULL)

  # one offset per condition (or a single shared value)
  offset <- get_dir_offset(fit, vars)

  tidy_fomo_draws(fit, vars_t, c("condition", "k")) |>
    mutate(offset = offset[pmin(condition, length(offset))],
           condition = label_condition(condition, cond_labels),
           k = factor(k, labels = direction_labels(max(k)))) |>
    select(.draw, condition, k, value, dist, offset)
}

get_fourier_harmonics <- function(fit, vars) {
  # the harmonics a Fourier model used (the same in every draw)
  vars_h <- vars[str_detect(vars, "^fourier_harmonics")]
  as.integer(get_draws(fit, vars_h)[1, vars_h])
}

label_fourier_coef <- function(draws, harmonics) {
  # fourier_coef's j-th column is {cos, sin} of harmonic (j + 1) %/% 2, as
  # in compute_directional_curve()
  pair_labs <- paste0("h", harmonics)
  draws |>
    mutate(harmonic = factor(pair_labs[(j + 1) %/% 2], levels = pair_labs),
           term = factor(if_else(j %% 2 == 1, "cos", "sin"), levels = c("cos", "sin"))) |>
    select(-j)
}

extract_fourier <- function(fit, vars, cond_labels) {

  # Fourier coefficients (and their prior draws): one {cos, sin} pair per
  # harmonic and condition, used directly as an additive log-weight on
  # angle - offset (see compute_directional_curve() for the curve itself)
  vars_f <- vars[str_detect(vars, "^(pr_)?fourier_coef\\[")]
  # von Mises fits have no fourier_coef
  if (length(vars_f) == 0) return(NULL)

  harmonics <- get_fourier_harmonics(fit, vars)
  offset <- get_dir_offset(fit, vars)

  tidy_fomo_draws(fit, vars_f, c("condition", "j")) |>
    label_fourier_coef(harmonics) |>
    mutate(offset = offset[pmin(condition, length(offset))],
           condition = label_condition(condition, cond_labels)) |>
    select(.draw, condition, harmonic, term, dist, value, offset)
}

direction_labels <- function(n) {
  # labels for n equally spaced directions, starting at 0
  k <- seq_len(n) - 1
  case_when(k == 0     ~ "0",
            k == n / 2 ~ "pi",
            k == 1     ~ paste0("pi/", n / 2),
            .default   = paste0(k, "*pi/", n / 2))
}

extract_random_effects <- function(fit, vars, cond_labels) {

  list(spatial  = extract_u_spatial(fit, vars, cond_labels),
       sigmas   = extract_sigmas(fit, vars, cond_labels),
       vonmises = extract_u_vonmises(fit, vars, cond_labels),
       fourier  = extract_u_fourier(fit, vars, cond_labels),
       cat      = extract_u_cat(fit, vars, cond_labels),
       first_v  = extract_u_first_v(fit, vars, cond_labels))
}

extract_sigmas <- function(fit, vars, cond_labels) {

  # random-effect SDs (sigma_u_*) and, where the model draws them, their
  # prior samples (pr_u_*): parameter is delta, psi, v, log_theta or
  # fourier_coef; the first index is the condition, and any further
  # indices (direction/coefficient, dimension, class) go in `index`
  vars_s <- vars[str_detect(vars, "^(sigma_u|pr_u)_[a-z_]+\\[")]
  if (length(vars_s) == 0) return(NULL)

  get_draws(fit, vars_s) |>
    pivot_longer(-.draw) |>
    mutate(dist = if_else(str_detect(name, "^pr_"), "prior", "posterior"),
           parameter = str_extract(name, "(?<=_u_)[a-z_]+"),
           idx = str_extract(name, "(?<=\\[)[0-9,]+"),
           condition = label_condition(as.integer(str_extract(idx, "^\\d+")), cond_labels),
           index = na_if(str_remove(idx, "^\\d+,?"), "")) |>
    select(.draw, parameter, condition, index, dist, value)
}

extract_u_spatial <- function(fit, vars, cond_labels) {

  vars_r <- vars[str_detect(vars, "^u_(psi|delta)\\[")]
  if (length(vars_r) == 0) return(NULL)

  # one column per parameter: rho_delta, rho_psi
  tidy_fomo_draws(fit, vars_r, c("condition", "person")) |>
    mutate(parameter = str_replace(parameter, "^u_", "rho_"),
           condition = label_condition(condition, cond_labels),
           person = factor(person)) |>
    select(-dist) |>
    pivot_wider(names_from = parameter, values_from = value)
}

extract_u_vonmises <- function(fit, vars, cond_labels) {

  vars_theta <- vars[str_detect(vars, "^u_theta")]
  if (length(vars_theta) == 0) return(NULL)

  offset <- get_dir_offset(fit, vars)

  # k indexes the von Mises components; phi is each component's preferred
  # direction in radians, as in the von Mises models' preferred directions:
  # 2*pi*(k - 1)/n_dircomps + dir_offset[condition]
  tidy_fomo_draws(fit, vars_theta, c("condition", "k", "person")) |>
    mutate(phi = 2 * pi * (k - 1) / max(k) + offset[pmin(condition, length(offset))],
           person = factor(person),
           condition = label_condition(condition, cond_labels)) |>
    select(.draw, condition, k, person, theta = value, phi)
}

extract_u_fourier <- function(fit, vars, cond_labels) {

  # per-person Fourier coefficients, labelled as in extract_fourier()
  vars_f <- vars[str_detect(vars, "^u_fourier_coef\\[")]
  if (length(vars_f) == 0) return(NULL)

  tidy_fomo_draws(fit, vars_f, c("condition", "j", "person")) |>
    label_fourier_coef(get_fourier_harmonics(fit, vars)) |>
    mutate(person = factor(person),
           condition = label_condition(condition, cond_labels)) |>
    select(.draw, person, condition, harmonic, term, value)
}

extract_u_first_v <- function(fit, vars, cond_labels) {

  # per-person first-selection class weights: u_first_v, the stationary
  # distribution of each person's transition matrix (probabilities)
  vars_i <- vars[str_detect(vars, "^u_first_v\\[")]
  if (length(vars_i) == 0) return(NULL)

  tidy_fomo_draws(fit, vars_i, c("condition", "person", "class")) |>
    mutate(person = factor(person),
           class = factor(class),
           condition = label_condition(condition, cond_labels)) |>
    select(.draw, person, condition, class, value)
}

extract_u_cat <- function(fit, vars, cond_labels) {

  vars_c <- vars[str_detect(vars, "^u_v\\[")]
  if (length(vars_c) == 0) return(NULL)

  tidy_fomo_draws(fit, vars_c, c("condition", "person", "prev_class", "class")) |>
    mutate(person = factor(person),
           class = factor(class),
           prev_class = factor(prev_class),
           condition = label_condition(condition, cond_labels)) |>
    select(.draw, person, condition, prev_class, class, dist, value)
}

# post$acc: accuracy_by_draw(), the accuracy per
# condition, selection number and draw, computed from the draws matrices
# instead of a long table (one row per selection and draw: about 100 M
# rows for kristjansson2014's stored P/W/loglik, past R's memory limit).
# Each of P, W and loglik is read as a draws x selections matrix, cut to
# the selections accuracy_exclude keeps, and averaged within each
# condition x n group by a matrix product. NULL (with a message) if the
# fit has no per-selection output
accuracy_by_draw <- function(fit, found, vars, stim = NULL, exclude = "first_and_forced") {

  if (!any(str_detect(vars, "^W\\["))) {
    if (any(str_detect(vars, "^log_lik\\["))) {
      message("no accuracy: the fit doesn't store per-selection predictions ",
              "(save_selections = FALSE, the default). For accuracy, use ",
              "extract_post(predict_fomo(fit), found, stim)")
    }
    return(NULL)
  }

  kept <- flag_selections(found, stim) |>
    mutate(row = row_number()) |>
    keep_for_accuracy(exclude)
  groups <- distinct(kept, condition, n) |> arrange(condition, n) |> mutate(g = row_number())
  kept <- left_join(kept, groups, by = c("condition", "n"))
  # selections x groups: 1 / (selections in the group) where it belongs, so
  # draws x selections %*% it gives each group's mean per draw
  avg <- matrix(0, nrow(kept), nrow(groups))
  avg[cbind(seq_len(nrow(kept)), kept$g)] <- 1 / tabulate(kept$g, nrow(groups))[kept$g]

  # one variable's draws x the kept selections (each read once: they can
  # be large)
  read_kept <- function(var) {
    m <- unclass(posterior::as_draws_matrix(fit$draws(var)))
    stopifnot("fit and data have different selections" = ncol(m) == nrow(found))
    m[, kept$row, drop = FALSE]
  }
  hit <- \(m) m == rep(kept$id, each = nrow(m))
  stoch <- hit(read_kept("P")) %*% avg
  determ <- hit(read_kept("W")) %*% avg
  loglik <- read_kept("loglik")
  mean_prob <- exp(loglik) %*% avg
  mean_loglik <- loglik %*% avg
  rm(loglik); gc(verbose = FALSE)

  tibble(condition = rep(groups$condition, each = nrow(stoch)),
         n = rep(groups$n, each = nrow(stoch)),
         .draw = rep(seq_len(nrow(stoch)), times = nrow(groups)),
         stochastic = as.vector(stoch), deterministic = as.vector(determ),
         mean_prob = as.vector(mean_prob), mean_loglik = as.vector(mean_loglik)) |>
    arrange(condition, n, .draw)
}

extract_gq <- function(fit, found, vars, stim = NULL) {

  gq_vars <- vars[str_detect(vars, "^(P|W|loglik)\\[")]
  if (length(gq_vars) == 0) {
    if (any(str_detect(vars, "^log_lik\\["))) {
      message("no accuracy: the fit doesn't store per-selection predictions ",
              "(save_selections = FALSE, the default). For accuracy, use ",
              "extract_post(predict_fomo(fit), found, stim)")
    }
    return(NULL)
  }

  get_draws(fit, gq_vars) |>
    pivot_longer(-.draw) |>
    mutate(row  = parse_number(name),
           name = str_extract(name, "^(P|W|loglik)")) |>
    pivot_wider(names_from = name, values_from = value) |>
    right_join(flag_selections(found, stim) |> mutate(row = row_number()) |> select(-x, -y),
               by = "row") |>
    select(-row)
}

# each selection's place in its trial: first (the trial's first pick, made
# without a previous item), items_left (including the one picked), chance
# (1 / items_left) and forced (only one item left). items_left, chance and
# forced need stim (each trial's number of items); without it they are NA
flag_selections <- function(found, stim = NULL) {

  found <- mutate(found, first = n == 1)
  if (is.null(stim)) return(mutate(found, items_left = NA_real_, chance = NA_real_, forced = NA))

  found |>
    left_join(count(stim, trial, name = "n_items"), by = "trial") |>
    mutate(items_left = n_items - n + 1,
           chance = 1 / items_left,
           forced = items_left == 1) |>
    select(-n_items)
}

# which selections accuracy leaves out (accuracy_exclude; flag_selections()
# columns):
#   "first_and_forced" (default): leave out each trial's first pick and any
#     forced last pick (one item left), which say nothing about the model
#   "none": every selection
# Add new rules here
keep_for_accuracy <- function(gq, exclude = "first_and_forced") {

  switch(exclude,
         none = gq,
         first_and_forced = {
           if (all(is.na(gq$forced))) {
             message("accuracy: no stim, so forced last picks can't be identified; ",
                     "leaving out first picks only (pass stim to extract_post())")
           }
           filter(gq, !first, !(forced %in% TRUE))
         },
         stop("unknown accuracy_exclude \"", exclude, "\": use \"first_and_forced\" or \"none\""))
}

# each kept selection's accuracy, averaged over draws: the flag_selections()
# columns plus stochastic (mean(P == id)) and deterministic (mean(W == id)),
# leaving out what accuracy_exclude says (keep_for_accuracy()). P and W
# come from the fit if it stored them (every draw), otherwise from
# predict_fomo() on n_draws draws, which needs the fit's data list (sl; by
# default the one fit_fomo() stored with it). The one place accuracy is
# computed from a fit (fomo_compare() and the example pages use it)
compute_selection_accuracy <- function(fit, d, sl = attr(fit, "fomo_sl"), exclude = "first_and_forced",
                               n_draws = 500, seed = NULL) {

  stored <- prod(fit$metadata()$stan_variable_sizes$P %||% 0) > 0
  stopifnot("no stored P/W and no data list for predict_fomo()" = stored || !is.null(sl))
  src <- if (stored) fit else predict_fomo(fit, sl, n_draws = n_draws, seed = seed)

  # one matrix at a time, compared with the selected items column by column
  hit_rate <- function(var) {
    m <- unclass(posterior::as_draws_matrix(src$draws(var)))
    stopifnot("fit and data have different selections" = ncol(m) == nrow(d$found))
    colMeans(m == rep(d$found$id, each = nrow(m)))
  }
  flag_selections(d$found, d$stim) |>
    mutate(stochastic = hit_rate("P"), deterministic = hit_rate("W")) |>
    keep_for_accuracy(exclude) |>
    structure(n_draws = posterior::ndraws(src$draws("W")),
              source = if (stored) "stored" else "predict_fomo")
}

