# =============================================================================
# BEHAVIOUR
#
# Statistics of selection data, observed or simulated: per-trial foraging
# statistics (compute_trial_stats()), steps and runs, features by selection
# number, and the chance and nearest-neighbour accuracy baselines.
# =============================================================================

compute_trial_stats <- function(found, per_trial = FALSE, which_stats = c("runs", "spatial")) {

  # Per-trial foraging summary statistics, following the descriptive
  # measures used in Clarke & Hughes (2026) "FoMo: a unifying theory of
  # visual foraging" (best-r, mean ITD, PAO and intersection rate, via the
  # ForagingOrg package). Works on real data (d$found)
  # or simulated data with multiple posterior draws (extract_pred() output,
  # tagged with .draw) - grouping columns are detected automatically, same
  # convention as compute_steps()/compute_runs().
  # per_trial = TRUE returns one row per trial; otherwise data with more
  # than one person are averaged to one row per person x condition (below).
  # which_stats picks the statistics (n_found is always included): any of
  # "runs" (n_runs, max_run_length; needs item classes v), "spatial"
  # (best-r, mean ITD, PAO) and "intersections" (number and rate of path
  # self-intersections), or "all".
  # The default is "runs" and "spatial". "intersections" is slow:
  # ForagingOrg::number.intersections() takes about 0.06 s per trial, which
  # dominates on posterior predictive draws (minutes for 50 draws of 50
  # trials)

  which_stats <- match.arg(which_stats, c("all", "runs", "spatial", "intersections"),
                           several.ok = TRUE)
  if ("all" %in% which_stats) which_stats <- c("runs", "spatial", "intersections")

  grps  <- trial_grouping_cols(found)
  has_v <- "v" %in% names(found) && (is.factor(found$v) || is.character(found$v))
  if (identical(which_stats, "runs") && !has_v) {
    stop("which_stats = \"runs\" needs item classes: a factor or character column v")
  }

  found |>
    arrange(across(all_of(grps)), n) |>
    group_by(across(all_of(grps))) |>
    summarise({

      n_found <- length(x)
      out <- data.frame(n_found = n_found)

      if ("spatial" %in% which_stats) {
        out$best_r   <- tryCatch(ForagingOrg::best.r(x, y), error = \(e) NA_real_)
        out$mean_itd <- if (n_found >= 2) {
          tryCatch(ForagingOrg::meanITD(x, y), error = \(e) NA_real_)
        } else NA_real_
        out$pao      <- tryCatch(ForagingOrg::PAO(x, y), error = \(e) NA_real_)
      }

      if ("intersections" %in% which_stats) {
        ni <- tryCatch(ForagingOrg::number.intersections(x, y),
                       error = \(e) list(Number.of.intersections = NA_real_,
                                         Intersection.rate = NA_real_))
        out$n_intersections   <- ni$Number.of.intersections
        out$intersection_rate <- ni$Intersection.rate
      }

      if ("runs" %in% which_stats && has_v) {
        runs <- run_lengths(v)
        out$n_runs         <- length(runs)
        out$max_run_length <- max(runs)
      }

      out

    }, .groups = "drop") -> dout

  # per-trial values of best_r in particular are dominated by
  # noise from the single stochastic selection path in that one trial (even
  # two simulations from identical true parameters barely correlate
  # trial-by-trial on those two metrics).
  # With multiple people, that noise washes out once averaged over each
  # person's trials, so collapse to one row per person x condition (per
  # .draw, if present) rather than returning noisy single-trial values.
  if (per_trial || n_distinct(found$person) < 2) return(dout)

  dout |>
    select(-matches("trial")) |>
    group_by(across(all_of(intersect(names(dout), c("person", "condition", ".draw"))))) |>
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
              .groups = "drop")
}

trial_grouping_cols <- function(found) {
  # the columns that identify a trial (in real or simulated data)
  names(found)[str_detect(names(found), "trial|person|condition|\\.draw")]
}

run_lengths <- function(v) {
  # lengths of the runs of same-class selections in one trial
  rle(as.vector(v))$lengths
}

compute_steps <- function(found) {

  # inter-item distance (delta) and direction (phi, 0 to 2 pi) of each
  # selection from the one before, for a dataframe (d$found) of trials
  found |>
    group_by(across(all_of(trial_grouping_cols(found)))) |>
    mutate(x0 = lag(x),
           y0 = lag(y),
           delta = sqrt((x - x0)^2 + (y - y0)^2),
           phi = atan2(y - y0, x - x0) %% (2 * pi)) |>
    filter(n > 1)
}

compute_runs <- function(found) {

  # run statistics (max run length and number of runs) for a dataframe
  # (d$found) of trials
  found |>
    group_by(across(all_of(trial_grouping_cols(found)))) |>
    summarise(max_run_length = max(run_lengths(v)),
              n_runs = length(run_lengths(v)),
              .groups = "drop")
}

summarise_over_trials <- function(found, feature) {

  # by position in the sequence (n) and whatever grouping columns found has:
  # the proportion of each class (item classes) or the mean (distance etc.)
  grps <- intersect(names(found), c("n", "person", "condition", ".draw"))

  if (is.factor(found[[feature]])) {
    found |>
      ungroup() |>
      count(across(all_of(c(grps, feature))), name = "value") |>
      rename(class = all_of(feature)) |>
      mutate(value = value / sum(value), .by = all_of(grps))
  } else {
    found |>
      group_by(across(all_of(grps))) |>
      summarise(value = mean(.data[[feature]]), .groups = "drop") |>
      pivot_longer(-any_of(grps)) |>
      filter(is.finite(value))
  }
}

#' Chance baseline from the items in each trial
#'
#' Expected accuracy of picking uniformly at random among the items still
#' on the display, for exhaustive search: 1 / (items in the trial - n + 1)
#' at the nth selection, for every n up to the number of items, averaged
#' over trials. Chance doesn't depend on when a trial stops, so the curve
#' covers every position even if the selections stop earlier.
#'
#' @param stim Data frame with one row per item: columns \code{trial} and
#'   \code{condition}
#' @param n_items Instead of stim: the number of items, one row per trial
#'   or condition (columns \code{condition} and \code{n_items})
#'
#' @return Tibble with columns condition, n and value
compute_chance_baseline <- function(stim = NULL, n_items = NULL) {

  if (!is.null(stim)) n_items <- count(stim, condition, trial, name = "n_items")

  # one row per trial (or condition) and selection number n
  n_items |>
    uncount(n_items, .id = "n", .remove = FALSE) |>
    summarise(value = mean(1 / (n_items - n + 1)), .by = c(condition, n)) |>
    arrange(condition, n)
}

#' Compute nearest-neighbour (proximity) baseline predictions
#'
#' At every selection step from the 2nd onward, predicts the next item as
#' whichever not-yet-selected item in that trial is closest (Euclidean
#' distance in x, y) to the item just selected, then checks whether that
#' greedy prediction matches what was actually selected next. A stronger,
#' data-driven point of comparison than \code{compute_chance_baseline()}'s
#' uniform-chance baseline.
#'
#' @param stim Data frame with one row per grid item: columns \code{trial},
#'   \code{id}, \code{x}, \code{y} (the full candidate set, found or not).
#' @param found Data frame with one row per selection, in order: columns
#'   \code{trial}, \code{id}, \code{n}, \code{x}, \code{y}, plus any grouping
#'   columns to keep (e.g. \code{condition}, \code{person}).
#' @param accuracy_exclude Which selections to leave out, as for accuracy
#'   (keep_for_accuracy()): "first_and_forced" (default; the first pick has
#'   no nearest neighbour anyway, and a forced last pick is always "right")
#'   or "none", every selection from the 2nd onward.
#'
#' @return Tibble with one row per kept selection step, containing \code{trial}, \code{n}, any of \code{condition} /
#'   \code{person} present in \code{found}, and \code{correct}: 1 if the
#'   actual selection was the nearest remaining item, 1 / k if it was one of
#'   k items tied for nearest, else 0.
#'
#' @seealso \code{\link{summarise_nn_baseline}}, \code{\link{plot_post_accuracy}}
compute_nn_baseline <- function(stim, found, accuracy_exclude = "first_and_forced") {

  # the baseline is for the observed data (one sequence per trial); several
  # posterior draws would be mixed together within each trial
  if (".draw" %in% names(found)) {
    stop("compute_nn_baseline() takes observed data only; `found` has a .draw ",
         "column (simulated/posterior draws).")
  }

  found <- arrange(found, trial, n)
  keep_cols <- intersect(names(found), c("trial", "n", "condition", "person"))
  stim_by_trial <- split(stim, stim$trial)

  map(split(found, found$trial), function(fseq) {

    if (nrow(fseq) < 2) return(NULL)
    st <- stim_by_trial[[as.character(fseq$trial[1])]]

    # for each step, is the next item the nearest remaining one? Ties (items
    # at the same distance, common on grids) count 1 / (number tied): the
    # rule's expected accuracy if it picks among them at random
    correct <- map_dbl(seq_len(nrow(fseq) - 1), function(ii) {
      d2 <- (fseq$x[ii] - st$x)^2 + (fseq$y[ii] - st$y)^2
      d2[st$id %in% fseq$id[seq_len(ii)]] <- Inf
      nearest <- st$id[d2 <= min(d2) + 1e-9]
      (fseq$id[ii + 1] %in% nearest) / length(nearest)
    })
    mutate(fseq[-1, keep_cols], correct = correct)
  }) |>
    list_rbind() |>
    # the same selections as accuracy (flag_selections(), keep_for_accuracy())
    left_join(flag_selections(found, stim) |> select(trial, n, first, forced),
              by = c("trial", "n")) |>
    keep_for_accuracy(accuracy_exclude) |>
    select(-first, -forced)
}

#' Summarise nearest-neighbour baseline accuracy
#'
#' Aggregates \code{compute_nn_baseline()}'s per-selection predictions into
#' the (condition, n, value) shape expected by \code{plot_post_accuracy()}'s
#' \code{nn_baseline} argument.
#'
#' @param nn Output of \code{\link{compute_nn_baseline}}.
#' @param grouping Character vector of columns to group by before averaging
#'   \code{correct}. Defaults to \code{c("condition", "n")}.
#'
#' @return Tibble with the grouping columns plus \code{value} (mean accuracy).
#'
#' @seealso \code{\link{compute_nn_baseline}}, \code{\link{plot_post_accuracy}}
summarise_nn_baseline <- function(nn, grouping = c("condition", "n")) {

  stopifnot(all(grouping %in% names(nn)))

  nn |>
    group_by(across(all_of(grouping))) |>
    summarise(value = mean(correct), .groups = "drop")
}
