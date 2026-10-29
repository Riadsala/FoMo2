################################################################################
# FOMO PLOTTING
#
# Naming: plot_<what>() plots selection data - observed, simulated or
# predicted (extract_pred()) - with the columns of import_data()'s found
# (plus stim where needed); where a second set can be overlaid, it is `sim`.
# plot_post_<what>() plots a fit: extract_post() output, or the fit itself.
# Helpers that build one piece of a plot are panel_*(); helpers that add
# layers to a plot are add_*(). Every plot_*() returns a ggplot or patchwork.
#
# Contents
# ------------------------------------------------------------------------------
# 1. Selection data (observed, simulated or predicted)
#    - plot_trial()        : Trials (one per condition by default): items and the selection path
#    - plot_features()     : A feature (distance, item class, ...) by selection
#    - plot_steps()        : Step lengths and directions by selection
#    - plot_runs()         : Run statistics (number of runs, longest run)
#    - plot_trial_stats()  : Observed vs simulated trial statistics
#
# 2. Fits
#    - plot_post()                     : Every population-level parameter
#    - plot_post_spatial()             : Spatial model: densities or weight curves
#    - plot_post_cat()                 : Class-transition probabilities
#    - plot_post_per_person_spatial()  : Per-person rho_delta and rho_psi
#    - plot_post_per_person_vonmises() : Per-person directional curves
#    - plot_post_directional_curves()  : Directional curves (compute_directional_curve())
#    - plot_post_accuracy()            : Accuracy and the chosen item's probability by selection
#    - plot_post_traces()              : Rank (or trace) plots of the worst-mixing parameters
#
# 3. Helpers: panel_densities(), panel_spatial_densities(),
#    panel_direction_densities(), panel_weight_curve(), panel_rose(),
#    panel_waves(), add_groundtruth(), add_hpdi_ribbons()
################################################################################

# conditions' colours (the default discrete palette in FoMo's plots)
FOMO_PALETTE <- as.character(paletteer_d("fishualize::Callanthias_australis"))

# FoMo's theme, with its discrete palette, added to every plot here (not set
# globally, so sourcing these functions leaves other ggplots unchanged)
theme_fomo <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          palette.colour.discrete = FOMO_PALETTE,
          palette.fill.discrete = FOMO_PALETTE)
}

fomo_ggplot <- function(...) ggplot(...) + theme_fomo()

# the colour that palette gives each condition (as the discrete scales do:
# in the order of the factor's levels), for when a plot can't map
# condition to colour
condition_colours <- function(conditions) {
  lv <- levels(factor(conditions))
  set_names(FOMO_PALETTE[seq_along(lv)], lv)
}

# models (plots comparing them) get their own colours, so they can't be
# mistaken for conditions: Okabe-Ito blue, bluish green, reddish purple,
# sky blue, vermillion, yellow (colour-blind safe)
MODEL_COLOURS <- c("#0072B2", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#E69F00")
scale_colour_model <- function(...) scale_colour_manual(values = MODEL_COLOURS, ...)
scale_fill_model <- function(...) scale_fill_manual(values = MODEL_COLOURS, ...)

# simulation values (groundtruth) in every plot: dashed, so they can't be
# mistaken for an estimate
TRUTH_LINE <- list(linetype = 2, linewidth = 0.7, alpha = 0.9)

# x-axis for absolute directions in radians, -pi to pi (atan2: 0 = right,
# pi/2 = up)
scale_x_direction <- function(name = "absolute direction\n(0 = right, pi/2 = up)") {
  scale_x_continuous(name, breaks = seq(-pi, pi, by = pi / 2),
                     labels = c("-pi", "-pi/2", "0", "pi/2", "pi"), expand = c(0, 0))
}

################################################################################
# 1. Selection data (observed, simulated or predicted)
################################################################################

#' Trials: every item, coloured (and sized) by its features, with the
#' selection path numbered in order and drawn in its condition's colour
#'
#' By default one random trial per condition, side by side (patchwork).
#' trial = one or more trial numbers plots those instead; one trial gives
#' a single ggplot, so plots still combine with +.
#'
#' @param found Selections (observed, simulated or predicted). With several
#'   draws (.draw), one is picked at random
#' @param stim Every item (default: found, i.e. every item was selected)
#' @param trial NA (default): one random trial per condition; or trial
#'   numbers
plot_trial <- function(found, stim = NULL, trial = NA) {

  if (is.null(stim)) stim <- found

  if (".draw" %in% names(found)) {
    found <- filter(found, .draw == sample(unique(found$.draw), 1))
  }

  trials <- if (all(is.na(trial))) {
    stim |>
      filter(trial %in% found$trial) |>
      distinct(condition, trial) |>
      slice_sample(n = 1, by = condition) |>
      arrange(condition) |>
      pull(trial)
  } else trial

  colours <- condition_colours(stim$condition)
  plots <- map(trials, \(t) panel_trial(filter(found, trial == t), filter(stim, trial == t), colours))

  if (length(plots) == 1) plots[[1]] else wrap_plots(plots)
}

panel_trial <- function(selections, items, colours) {

  # one trial for plot_trial()
  item_features <- if ("v" %in% names(selections)) "cat" else "none"

  plt <- switch(item_features,
    cat = fomo_ggplot(selections, aes(x, y, colour = v, size = 5)) +
      scale_colour_paletteer_d("PrettyCols::Bold"),
    none = fomo_ggplot(selections, aes(x, y)) +
      scale_size(range = c(1, 10)))

  plt <- plt +
    geom_point(data = items) +
    theme_minimal() +
    ggtitle(paste0("person: ", unique(items$person), ", condition: ",
                   unique(items$condition), ", trial: ", unique(items$trial))) +
    theme(plot.title = element_text(size = 10))

  if ("n" %in% names(selections)) {
    # inherit.aes = FALSE: the points may map size, which lines warn about
    plt <- plt +
      geom_path(data = filter(selections, is.finite(n)), aes(x, y),
                inherit.aes = FALSE, colour = colours[[as.character(unique(items$condition))]],
                linewidth = 0.8) +
      geom_text(aes(label = n), size = 2.5, colour = "grey20",
                position = position_nudge(x = -0.02, y = -0.02))
  }

  plt +
    coord_equal() +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          plot.background = element_rect(fill = "darkgrey", colour = "black"),
          panel.grid = element_blank())
}

#' A feature by selection number (1st, 2nd, ... selection in a trial)
#'
#' found is drawn as points with a loess line per person; sim as a median
#' line with its 97% HPDI over draws. Item classes are shown as the
#' proportion of each class, one facet per class; numeric features as the
#' mean.
#'
#' @param found Selections (with the feature as a column)
#' @param feature Column to plot: "delta" (step length; see
#'   compute_steps()), "v" (item features), ...
#' @param sim Simulated or predicted selections, with .draw
plot_features <- function(found = NULL, feature = "delta", sim = NULL) {

  y_label <- if (feature == "delta") expression(delta) else feature

  plt <- fomo_ggplot(mapping = aes(n, value, colour = condition, fill = condition,
                              group = interaction(person, condition))) +
    single_condition_guides((found %||% sim)$condition)

  if (!is.null(found)) {
    found <- summarise_over_trials(found, feature)
    plt <- plt +
      geom_point(data = found, alpha = 0.5) +
      geom_smooth(data = found, method = "loess", formula = "y ~ x",
                  se = FALSE, span = 1.5)
  }

  if (!is.null(sim)) {
    sim <- summarise_over_trials(sim, feature)
    plt <- plt +
      stat_lineribbon(data = sim, alpha = 0.5, .width = max(FOMO_HPDI_WIDTHS),
                      point_interval = "median_hdci", colour = NA)
  }

  # one facet per class (item classes)
  if ("class" %in% names(found %||% sim)) plt <- plt + facet_wrap(~class)

  plt + labs(y = y_label, x = expression(italic(n) * "th selected item"))
}

#' Step lengths by selection number, and the distribution of step
#' directions
#'
#' @param found Selections; the steps (delta, phi) are computed with
#'   compute_steps() unless already there
#' @param sim Simulated or predicted selections, shown alongside as "model"
#' @param direction_plot "waves" (smoothed densities) or "rose" (polar
#'   histogram)
plot_steps <- function(found, sim = NULL, direction_plot = c("waves", "rose")) {

  direction_plot <- match.arg(direction_plot)

  steps <- if ("delta" %in% names(found)) found else compute_steps(found)
  if (!is.null(sim)) {
    steps <- bind_rows(mutate(steps, dataset = "data"),
                       mutate(compute_steps(sim), dataset = "model"))
  }

  step_lengths <- plot_features(steps, feature = "delta")
  step_directions <- switch(direction_plot,
                            rose  = panel_rose(steps, fill = "condition"),
                            waves = panel_waves(steps))

  step_lengths + step_directions + guides(colour = "none") + plot_layout(guides = "collect")
}

#' Number of runs against the longest run, one point per trial
#'
#' @param runs Output of compute_runs()
plot_runs <- function(runs) {

  fomo_ggplot(runs, aes(n_runs, max_run_length, colour = condition)) +
    geom_jitter(height = 0.1, width = 0.1, alpha = 0.5) +
    single_condition_guides(runs$condition) +
    labs(x = "number of runs", y = "longest run")
}

#' Observed vs simulated trial statistics
#'
#' Posterior predictive check: each statistic from compute_trial_stats()
#' (best-r, PAO, ...) on the observed data against the same
#' statistic on simulated data, one facet per statistic, with y = x
#' (dashed) and a linear fit (solid). With several draws in sim, each point
#' is the median over draws with its HPDI as an error bar.
#'
#' @param obs compute_trial_stats() of the observed data
#' @param sim compute_trial_stats() of simulated data (optionally with .draw)
#' @param metrics Statistics to plot (default: every one in both)
#' @param ci_width Width of the error bars (default: the wide
#'   FOMO_HPDI_WIDTHS, 97%)
plot_trial_stats <- function(obs, sim, metrics = NULL, ci_width = max(FOMO_HPDI_WIDTHS)) {

  join_keys <- intersect(names(obs), c("trial", "person", "condition"))
  stopifnot(length(join_keys) > 0)

  if (is.null(metrics)) {
    all_metrics <- c("best_r", "mean_itd", "pao", "intersection_rate",
                     "n_runs", "max_run_length")
    metrics <- intersect(all_metrics, intersect(names(obs), names(sim)))
  }

  sim |>
    select(all_of(join_keys), any_of(".draw"), all_of(metrics)) |>
    pivot_longer(all_of(metrics), names_to = "metric", values_to = "simulated") |>
    group_by(across(all_of(c(join_keys, "metric")))) |>
    median_hdci(simulated, .width = ci_width) |>
    ungroup() -> sim_summary

  obs |>
    select(all_of(join_keys), all_of(metrics)) |>
    pivot_longer(all_of(metrics), names_to = "metric", values_to = "observed") |>
    inner_join(sim_summary, by = c(join_keys, "metric")) |>
    fomo_ggplot(aes(observed, simulated,
               colour = if ("condition" %in% join_keys) condition else NULL)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = .lower, ymax = .upper), alpha = 0.25, width = 0) +
    geom_point(alpha = 0.5) +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x,
                se = FALSE, colour = "black") +
    facet_wrap(~metric, scales = "free", nrow = 1) +
    labs(x = "observed", y = "simulated (posterior predictive)", colour = "condition") +
    single_condition_guides(obs$condition)
}

################################################################################
# 2. Fits
################################################################################

#' Every population-level parameter: the item-class parameters (cat
#' models) above the spatial model's (plot_post_spatial())
#'
#' @param post Output of extract_post()
#' @param groundtruth Optional simulation values (e.g. the simulated examples' simparams)
#' @param spatial "densities" or "weights": see plot_post_spatial()
plot_post <- function(post, groundtruth = NULL, spatial = c("densities", "weights")) {

  spatial <- match.arg(spatial)

  # one panel per part, stacked, each as tall as its rows of facets (so
  # every facet row has the same height)
  panels <- if (spatial == "densities") {
    list(panel_spatial_densities(post, groundtruth), panel_direction_densities(post, groundtruth))
  } else {
    list(plot_post_spatial(post, groundtruth, show = "weights"))
  }
  rows <- if (spatial == "densities") c(1, if (is.null(post$fourier)) 1 else 2) else 1

  if (!is.null(post$cat)) {
    panels <- c(list(plot_post_cat(post, groundtruth)), panels)
    rows <- c(n_distinct(post$cat$prev_class), rows)
  }

  # one legend: the panels' fills differ only in alpha (none with one
  # condition)
  conditions <- post$spatial$condition
  wrap_plots(panels, ncol = 1, heights = rows, guides = "collect") &
    (single_condition_guides(conditions) %||%
       guides(fill = guide_legend(override.aes = list(alpha = 0.4))))
}

#' The spatial model's parameters: distance (rho_delta), turn (rho_psi)
#' and direction (von Mises theta or Fourier coefficients)
#'
#' show = "densities" ("d"): posterior densities over the prior's HPDIs
#' (grey), with the simulation values as dashed lines.
#' show = "weights" ("w"): what each parameter does to an item's selection
#' weight (compute_weight_curves()): the negative exponential in distance,
#' the same in turn angle, and the directional curve, as 53% and 97% HPDI
#' bands; simulation values as dashed lines.
#'
#' @param post Output of extract_post()
#' @param groundtruth Optional simulation values (e.g. the simulated examples' simparams)
#' @param show "densities" or "weights" (or "d" / "w")
#' @param xy_scaling As in prep_stan_list() (default 10); only for weights
plot_post_spatial <- function(post, groundtruth = NULL, show = c("densities", "weights"),
                              xy_scaling = 10) {

  show <- match.arg(show)

  if (show == "densities") {
    return(panel_spatial_densities(post, groundtruth) /
             panel_direction_densities(post, groundtruth))
  }

  curves <- compute_weight_curves(post, xy_scaling)
  truth <- if (!is.null(groundtruth)) {
    compute_weight_curves(groundtruth_as_post(groundtruth, post), xy_scaling)
  }

  # one panel per term, with its own x-axis
  x_scales <- list(
    distance  = scale_x_continuous("distance\n(display widths)"),
    turn      = scale_x_continuous("turn from previous step\n(0 = straight on, pi = U-turn)",
                                   breaks = seq(0, pi, by = pi / 4),
                                   labels = c("0", "pi/4", "pi/2", "3pi/4", "pi"), expand = c(0, 0)),
    direction = scale_x_direction())

  imap(x_scales, function(x_scale, which_term) {
    panel_weight_curve(filter(curves, term == which_term),
                       if (!is.null(truth)) filter(truth, term == which_term),
                       x_scale)
  }) |>
    wrap_plots() +
    plot_layout(guides = "collect")
}

#' Class-transition probabilities (b_v), over the prior's HPDIs
#'
#' @param post Output of extract_post() for a cat model
#' @param groundtruth Optional simulation values (b_v_pr)
plot_post_cat <- function(post, groundtruth = NULL) {

  if (!any(post$cat$dist == "prior")) warning("no prior draws for the cat model")

  post$cat |>
    mutate(class = paste("class", class),
           prev_class = paste("prev. class", prev_class)) |>
    panel_densities(facet_grid(prev_class ~ class, scales = "free")) |>
    add_groundtruth(groundtruth, "cat", levels(post$cat$condition)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.33)
}

#' Per-person rho_delta and rho_psi (intervals), over the population-level
#' HPDIs (grey bands)
#'
#' @param post Output of extract_post() for a multi-level fit
#' @param groundtruth Optional: groundtruth_from_sim(), whose per-person
#'   values are drawn as squares
plot_post_per_person_spatial <- function(post, groundtruth = NULL) {

  post$u$spatial |>
    pivot_longer(c(rho_psi, rho_delta), names_to = "parameter", values_to = "value") -> per_person

  post$spatial |>
    filter(dist == "posterior") |>
    group_by(parameter, condition) |>
    median_hdci(value, .width = FOMO_HPDI_WIDTHS) -> population

  fomo_ggplot(per_person, aes(value, person, colour = condition)) +
    geom_rect(data = population, aes(xmin = .lower, xmax = .upper, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, alpha = 0.5) +
    stat_interval(.width = FOMO_HPDI_WIDTHS, point_interval = "median_hdci", alpha = 0.5) +
    facet_grid(rows = condition_rows(per_person$condition), cols = vars(parameter), scales = "free") +
    single_condition_guides(per_person$condition) -> plt

  if (is.null(groundtruth)) return(plt)

  groundtruth$u_spatial |>
    pivot_longer(c(rho_delta, rho_psi), names_to = "parameter", values_to = "value") -> true_per_person

  plt <- plt + geom_point(data = true_per_person, colour = "black", shape = 0, size = 3)
  add_groundtruth(plt, groundtruth, "spatial", levels(post$spatial$condition))
}

#' Per-person class-transition probabilities (u_v, intervals), over the
#' population-level HPDIs (grey bands): one column per transition (previous
#' class -> next class), one row per condition
#'
#' @param post Output of extract_post() for a multi-level cat fit
#' @param groundtruth Optional: groundtruth_from_sim(), whose per-person
#'   values are drawn as squares
plot_post_per_person_cat <- function(post, groundtruth = NULL) {

  stopifnot("needs extract_post() of a multi-level cat fit" = !is.null(post$u$cat))

  label_transition <- \(d) mutate(d, transition = factor(paste(prev_class, "->", class)))

  per_person <- label_transition(post$u$cat)

  post$cat |>
    filter(dist == "posterior") |>
    label_transition() |>
    group_by(transition, condition) |>
    median_hdci(value, .width = FOMO_HPDI_WIDTHS) -> population

  fomo_ggplot(per_person, aes(value, person, colour = condition)) +
    geom_rect(data = population, aes(xmin = .lower, xmax = .upper, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, alpha = 0.5) +
    stat_interval(.width = FOMO_HPDI_WIDTHS, point_interval = "median_hdci", alpha = 0.5) +
    facet_grid(rows = condition_rows(per_person$condition), cols = vars(transition)) +
    single_condition_guides(per_person$condition) +
    scale_x_continuous("probability of the next class", breaks = c(0, 0.5, 1)) +
    coord_cartesian(xlim = c(0, 1)) -> plt

  if (is.null(groundtruth$u_cat)) return(plt)

  plt + geom_point(data = label_transition(groundtruth$u_cat), colour = "black", shape = 0, size = 3)
}

#' Each person's directional preference curve (von Mises fits): the sum
#' over directions of theta * von Mises(angle; preferred direction,
#' FOMO_KAPPA), for 25 posterior draws, one row per person
#'
#' @param post Output of extract_post() for a multi-level von Mises fit
#' @param groundtruth Optional: groundtruth_from_sim(); true curves in black
#' @param yscale Vertical space per person
plot_post_per_person_vonmises <- function(post, groundtruth = NULL, yscale = 50) {

  person_curves <- function(u) {
    u |>
      crossing(angle = seq(-pi, pi, length.out = 50)) |>
      mutate(y = theta * vonmises_density(angle, phi, FOMO_KAPPA)) |>
      summarise(y = sum(y), .by = c(any_of(".draw"), condition, person, angle)) |>
      mutate(offset = as.numeric(person) * yscale)
  }

  some_draws <- sample_draw_ids(post$u$vonmises$.draw, 25)
  curves <- person_curves(filter(post$u$vonmises, .draw %in% some_draws))

  true_curves <- if (!is.null(groundtruth)) {
    geom_path(data = person_curves(groundtruth$u_vonmises),
              aes(group = interaction(condition, person)),
              linewidth = 1, colour = "black")
  }

  n_people <- n_distinct(curves$person)

  fomo_ggplot(curves, aes(angle, y + offset, colour = condition)) +
    geom_hline(aes(yintercept = offset), colour = "black") +
    geom_vline(xintercept = seq(-pi, pi, by = pi / 2), linetype = 2, alpha = 0.5) +
    geom_path(alpha = 0.25, aes(group = interaction(.draw, condition, person))) +
    scale_y_continuous("person",
                       breaks = seq(yscale, yscale * n_people, by = yscale),
                       labels = levels(curves$person)) +
    scale_x_direction() +
    coord_cartesian(ylim = c(yscale, yscale * (n_people + 1))) +
    true_curves +
    (if (n_distinct(curves$condition) > 1) facet_wrap(~condition, ncol = 2)) +
    single_condition_guides(curves$condition)
}

#' Directional log-weight curves: 53% and 97% HPDI bands, one facet per
#' condition (no facets with one condition)
#'
#' @param curves compute_directional_curve() output; with a `model`
#'   column (several fits bound together), one colour per model
#' @param truth Optional true curve (compute_directional_curve() of the
#'   simulation), drawn in black
plot_post_directional_curves <- function(curves, truth = NULL) {

  if (!("model" %in% names(curves))) curves$model <- "fit"

  curves |>
    group_by(model, condition, angle) |>
    median_hdci(log_weight, .width = FOMO_HPDI_WIDTHS) -> summary

  plt <- fomo_ggplot(summary, aes(angle, log_weight, colour = model, fill = model)) |>
    add_hpdi_ribbons(summary) +
    scale_colour_model() +
    scale_fill_model() +
    scale_x_direction() +
    scale_y_continuous("log-weight (centred)")
  if (n_distinct(summary$condition) > 1) plt <- plt + facet_wrap(~condition)

  if (!is.null(truth)) {
    plt <- plt + geom_line(data = truth, aes(angle, log_weight), inherit.aes = FALSE,
                           colour = "black", linetype = TRUTH_LINE$linetype,
                           linewidth = TRUTH_LINE$linewidth, alpha = TRUTH_LINE$alpha)
  }
  if (n_distinct(summary$model) == 1) plt <- plt + theme(legend.position = "none")

  plt
}

#' Accuracy by selection number
#'
#' By default deterministic accuracy (the model's most likely item is the
#' one chosen) and the mean probability the model gave the chosen item
#' (mean(exp(loglik))), as 53% and 97% HPDIs, with chance (dotted; the mean
#' of 1 / items left, so on the same scale as the probability) and, if
#' stim and found are given, the nearest-neighbour baseline (dashed;
#' compute_nn_baseline()). With stim and found, both baselines use the
#' observed selections that accuracy_exclude keeps, as post$acc does (pass
#' the value given to extract_post()).
#'
#' @param acc post$acc (extract_post(), via accuracy_by_draw())
#' @param stim Every item, for the baselines
#' @param found The observed selections, for the baselines
#' @param nn_baseline Draw the nearest-neighbour baseline?
#' @param accuracy_exclude Which selections to leave out (keep_for_accuracy())
#' @param metrics Any of "deterministic", "probability" (the default two),
#'   "stochastic" (a sampled item is the one chosen: on average the same as
#'   probability, with sampling noise) and "likelihood" (the geometric mean
#'   probability, exp(mean(loglik)))
plot_post_accuracy <- function(acc, stim = NULL, found = NULL, nn_baseline = TRUE,
                               accuracy_exclude = "first_and_forced",
                               metrics = c("deterministic", "probability")) {

  if (is.null(acc) || nrow(acc) == 0) {
    stop("no accuracy to plot: the fit doesn't store per-selection predictions ",
         "(save_selections = FALSE, the default). Get them with predict_fomo():\n",
         "  pred <- predict_fomo(fit)\n",
         "  post <- extract_post(pred, found, stim)\n",
         "  plot_post_accuracy(post$acc, stim, found)",
         call. = FALSE)
  }
  stopifnot("acc should be extract_post()$acc" = all(c("condition", "n") %in% names(acc)))
  metric_labels <- c(deterministic = "deterministic accuracy",
                     probability = "probability of the chosen item",
                     stochastic = "stochastic accuracy",
                     likelihood = "likelihood (geometric mean)")
  metrics <- match.arg(metrics, names(metric_labels), several.ok = TRUE)
  stopifnot("acc has no mean_prob: rerun extract_post()" =
              !("probability" %in% metrics) || "mean_prob" %in% names(acc))

  # chance, 1 / (items left): over the selections actually made if stim and
  # found are given; otherwise for exhaustive search, from the items in
  # each trial (stim) or assuming as many items as the most selections in a
  # trial
  chance <- if (!is.null(stim) && !is.null(found)) {
    flag_selections(found, stim) |>
      keep_for_accuracy(accuracy_exclude) |>
      summarise(value = mean(chance), .by = c(condition, n)) |>
      arrange(condition, n)
  } else if (!is.null(stim)) {
    compute_chance_baseline(stim)
  } else {
    compute_chance_baseline(n_items = summarise(acc, n_items = max(n), .by = condition))
  }

  if (!("mean_prob" %in% names(acc))) acc$mean_prob <- NA_real_  # older acc
  acc |>
    mutate(likelihood = exp(mean_loglik), probability = mean_prob) |>
    pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") |>
    mutate(metric = factor(metric_labels[metric], levels = metric_labels[metrics])) |>
    group_by(condition, n, metric) |>
    median_hdci(value, .width = FOMO_HPDI_WIDTHS) -> accuracy

  fomo_ggplot(accuracy, aes(n, value, fill = condition, group = interaction(condition, .width))) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5, linewidth = 0, colour = NA) +
    geom_path(data = chance, aes(group = condition), linetype = "dotted", colour = "black") +
    facet_wrap(~metric) +
    single_condition_guides(accuracy$condition) -> plt

  if (!nn_baseline) return(plt)
  if (is.null(stim) || is.null(found)) {
    message("no stim and found, so no nearest-neighbour baseline")
    return(plt)
  }

  # the nearest-neighbour baseline on the same selections as accuracy
  nn <- summarise_nn_baseline(compute_nn_baseline(stim, found, accuracy_exclude))

  plt +
    geom_path(data = nn,
              aes(colour = condition, group = condition),
              linetype = "dashed", alpha = 0.5)
}

#' MCMC diagnostics: rank plots (or traces) of the worst-mixing parameters
#'
#' The n parameters, from anywhere in the model (random effects included),
#' that mixed worst, in one row. "Worst" combines the two diagnostics, as
#' each catches what the other misses (chains that disagree: R-hat; chains
#' that each mix slowly: ESS): a parameter's score is the larger of
#' (R-hat - 1) / (FOMO_RHAT_OK[1] - 1) and FOMO_ESS_OK[2] / its smaller
#' (bulk or tail) ESS, so below 1 is fine on both. Each panel is labelled
#' with its R-hat and ESS, and coloured green, orange or red for good,
#' acceptable or bad (as in fomo_summary()'s table).
#'
#' Rank plots (type = "rank", bayesplot::mcmc_rank_overlay()): a histogram
#' of each chain's draws' ranks among all draws; well-mixed chains are all
#' flat (uniform), and a chain above or below the rest is stuck in part of
#' the posterior.
#'
#' @param fit A fit (fit_fomo())
#' @param n How many parameters
#' @param variables Parameter names: plot every element of these instead of
#'   the worst n (e.g. c("rho_delta", "fourier_coef"))
#' @param type "rank" (default) or "trace"
plot_post_traces <- function(fit, n = 5, variables = NULL, type = c("rank", "trace")) {

  type <- match.arg(type)

  vars <- if (is.null(variables)) {
    fomo_variables(fit, random_effects = TRUE)
  } else {
    all_vars <- fomo_variables(fit, random_effects = TRUE)
    all_vars[str_detect(all_vars, paste0("^(", paste(variables, collapse = "|"), ")(\\[|$)"))]
  }
  stopifnot("no such parameters in the fit" = length(vars) > 0)

  fit$summary(variables = vars, "rhat", "ess_bulk", "ess_tail") |>
    mutate(score = pmax((rhat - 1) / (FOMO_RHAT_OK[1] - 1),
                        FOMO_ESS_OK[2] / pmin(ess_bulk, ess_tail))) -> diag

  if (is.null(variables)) diag <- slice_max(diag, score, n = n, with_ties = FALSE)

  # good / acceptable / bad, as in fomo_summary()'s colours
  diag <- mutate(diag,
                 ess = pmin(ess_bulk, ess_tail),
                 status = case_when(rhat >= FOMO_RHAT_OK[2] | ess < FOMO_ESS_OK[1] ~ "red",
                                    rhat >= FOMO_RHAT_OK[1] | ess < FOMO_ESS_OK[2] ~ "orange",
                                    .default = "green"))

  # one plot per parameter, each in its own colour scheme (bayesplot's
  # schemes are global and fixed when a plot is made)
  previous_scheme <- bayesplot::color_scheme_get()
  on.exit(bayesplot::color_scheme_set(unname(unlist(previous_scheme))))

  plots <- pmap(diag, function(variable, rhat, ess_bulk, ess_tail, status, ...) {
    draws <- fit$draws(variables = variable)
    bayesplot::color_scheme_set(status)
    plt <- if (type == "rank") bayesplot::mcmc_rank_overlay(draws) else bayesplot::mcmc_trace(draws)
    # few x ticks: the panels are narrow
    plt <- suppressMessages(plt + scale_x_continuous(n.breaks = 3))
    # one parameter: no facet strip, so the label goes in the title
    plt +
      labs(title = sprintf("%s\nR-hat %.3f, ESS %.0f / %.0f", variable, rhat, ess_bulk, ess_tail)) +
      theme(legend.position = "none", plot.title = element_text(size = 7, hjust = 0.5))
  })

  wrap_plots(plots, nrow = if (is.null(variables)) 1 else NULL) +
    plot_annotation(title = if (is.null(variables)) paste("The", nrow(diag), "worst-mixing parameters (R-hat and ESS)"))
}

################################################################################
# 3. Helpers
################################################################################

# With one condition, a condition legend or a row of condition facets says
# nothing: condition_rows() gives facet_grid() its rows only when there are
# several conditions, and single_condition_guides() hides the condition
# legend when there is one (ggplot ignores the NULL otherwise)
condition_rows <- function(conditions) {
  if (n_distinct(conditions) > 1) vars(condition)
}

single_condition_guides <- function(conditions) {
  if (n_distinct(conditions) == 1) guides(colour = "none", fill = "none")
}

panel_densities <- function(d, facets, fill_alpha = 0.33) {

  # posterior densities by condition, over the prior's HPDIs (grey; assumed
  # the same in every condition)
  prior_hpdi <- compute_prior_hpdi(d)

  plt <- fomo_ggplot(filter(d, dist == "posterior"), aes(value))

  if (nrow(prior_hpdi) > 0) {
    plt <- plt +
      geom_rect(data = prior_hpdi, aes(ymin = 0, ymax = Inf, xmin = .lower, xmax = .upper),
                alpha = 0.15)
  }

  plt +
    geom_density(aes(fill = condition), alpha = fill_alpha) +
    facets +
    single_condition_guides(d$condition)
}

panel_spatial_densities <- function(post, groundtruth) {

  # rho_delta and rho_psi
  post$spatial |>
    panel_densities(facet_wrap(~parameter, scales = "free")) |>
    add_groundtruth(groundtruth, "spatial", levels(post$spatial$condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.33)
}

panel_direction_densities <- function(post, groundtruth) {

  # von Mises: log(theta) per direction. Fourier: each coefficient
  if (!is.null(post$vonmises)) {
    return(post$vonmises |>
             mutate(value = log(value)) |>
             panel_densities(facet_wrap(~k, scales = "free", nrow = 1), fill_alpha = 0.5) |>
             add_groundtruth(groundtruth, "vonmises", levels(post$vonmises$condition)))
  }

  plt <- panel_densities(post$fourier, facet_grid(term ~ harmonic, scales = "free"), fill_alpha = 0.5) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3))
  if (is.null(groundtruth$fourier_coef_mu)) return(plt)

  plt + geom_vline(data = groundtruth_as_post(groundtruth, post)$fourier,
                   aes(xintercept = value, colour = condition), linetype = TRUTH_LINE$linetype,
                   linewidth = TRUTH_LINE$linewidth, alpha = TRUTH_LINE$alpha)
}

panel_weight_curve <- function(curves, truth, x_scale) {

  # one term's weight curve (compute_weight_curves()): HPDI bands per
  # condition, with the true curve dashed
  curves |>
    group_by(condition, x) |>
    median_hdci(weight, .width = FOMO_HPDI_WIDTHS) -> summary

  plt <- fomo_ggplot(summary, aes(x, weight, colour = condition, fill = condition)) |>
    add_hpdi_ribbons(summary) +
    x_scale +
    scale_y_continuous("relative weight") +
    single_condition_guides(summary$condition)

  if (is.null(truth)) return(plt)

  plt + geom_line(data = truth, aes(x, weight, colour = condition), inherit.aes = FALSE,
                  linetype = TRUTH_LINE$linetype, linewidth = TRUTH_LINE$linewidth,
                  alpha = TRUTH_LINE$alpha)
}

panel_rose <- function(steps, fill = NULL, n_bins = 8) {

  # polar histogram of the step directions (phi, 0 to 2 pi)
  stopifnot(n_bins > 0)

  fomo_ggplot(steps, aes(x = phi, fill = if (!is.null(fill)) .data[[fill]] else NULL)) +
    geom_histogram(aes(y = after_stat(count / sum(count))),
                   breaks = (0:n_bins - 0.5) / n_bins * 2 * pi,
                   alpha = 0.25, position = "identity") +
    coord_polar(start = -pi / n_bins, direction = -1) +
    scale_x_continuous("step directions", breaks = seq(0, 3 * pi / 2, by = pi / 2)) +
    theme(axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    labs(fill = fill) -> plt

  if ("dataset" %in% names(steps)) plt <- plt + facet_wrap(~dataset)
  plt
}

panel_waves <- function(steps, bw = 25) {

  # smoothed (von Mises kernel) density of the step directions (phi), per
  # condition
  circular_density <- function(d) {
    angles <- circular::circular(d$phi, type = "angles", units = "radians",
                       modulo = "2pi", zero = 0, rotation = "counter")
    # density.circular() coerces its own output grid and always warns
    # about it; harmless
    density <- suppressWarnings(circular::density.circular(angles, kernel = "vonmises", bw = bw,
                                                 from = 0, to = 2 * pi, n = 512))
    tibble(phi = as.numeric(density$x), density = density$y)
  }

  steps |>
    group_by(condition) |>
    group_modify(\(d, key) circular_density(d)) |>
    ungroup() |>
    fomo_ggplot(aes(phi, density, colour = condition)) +
    geom_vline(xintercept = c(0, pi / 2, pi, 3 * pi / 2, 2 * pi), linetype = "dashed") +
    geom_line(linewidth = 2) +
    scale_x_continuous(expression(phi),
                       breaks = seq(0, 2 * pi, by = pi / 2),
                       labels = c("0", expression(pi / 2), expression(pi),
                                  expression(3 * pi / 2), expression(2 * pi)),
                       limits = c(0, 2 * pi), expand = c(0, 0)) +
    labs(y = "density", colour = "condition") +
    single_condition_guides(steps$condition)
}

add_groundtruth <- function(plt, groundtruth, params, cond_labels) {

  # simulation values (groundtruth_tables()) as vertical lines, coloured by
  # condition if the ground truth has more than one
  stopifnot(params %in% c("spatial", "vonmises", "cat"))

  if (is.null(groundtruth)) return(plt)

  truth <- groundtruth_tables(groundtruth, cond_labels)[[params]]
  line_aes <- if ("condition" %in% names(truth)) aes(xintercept = value, colour = condition) else aes(xintercept = value)

  plt + geom_vline(data = truth, line_aes, linetype = TRUTH_LINE$linetype,
                   linewidth = TRUTH_LINE$linewidth, alpha = TRUTH_LINE$alpha)
}

add_hpdi_ribbons <- function(plt, summary) {

  # 53% and 97% HPDI bands (FOMO_HPDI_WIDTHS), from median_hdci() output;
  # x, y and colours come from the plot's aes(). No median line: the
  # posterior is the bands, not one best-fitting curve
  wide <- filter(summary, .width == max(FOMO_HPDI_WIDTHS))
  narrow <- filter(summary, .width == min(FOMO_HPDI_WIDTHS))

  plt +
    geom_ribbon(data = wide, aes(ymin = .lower, ymax = .upper), alpha = 0.15, colour = NA) +
    geom_ribbon(data = narrow, aes(ymin = .lower, ymax = .upper), alpha = 0.3, colour = NA)
}

#' The loo difference between two models by condition
#'
#' compute_elpd_by_condition() as points with 2 SE bars: where in the
#' data one model beats the other.
#'
#' @param loos A named list of two loo objects (first - second)
#' @param stim The data both models were fitted to
plot_post_elpd_by_condition <- function(loos, stim) {
  compute_elpd_by_condition(loos, stim) |>
    fomo_ggplot(aes(elpd_diff, condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_pointrange(aes(xmin = elpd_diff - 2 * se, xmax = elpd_diff + 2 * se)) +
    labs(x = paste0("elpd difference (", names(loos)[1], " - ", names(loos)[2], ")"), y = NULL)
}
