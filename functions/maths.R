# =============================================================================
# MATHS
#
# R versions of the models' maths, for plotting and checks: softmax, the von
# Mises density, reference-class log-odds for class transitions, and small
# numeric helpers for generate_stim().
# =============================================================================

softmax <- function(x) {
  z <- x - max(x)
  exp_z <- exp(z)
  exp_z / sum(exp_z)
}

vonmises_density <- function(x, mu, kappa) {
  # the von Mises density, written to avoid overflow at large kappa
  # (vectorised over x and mu, unlike circular::dvonmises())
  exp(kappa * (cos(x - mu) - 1)) / (2 * pi * besselI(kappa, 0, expon.scaled = TRUE))
}

#' Reference-class log-odds for one previous class
#'
#' R port of get_logratio_mean() in models/include/FoMo_functions.stan -
#' converts one column of a class-transition probability matrix into
#' log-odds relative to the last class. Keep in sync with the Stan version
#' if that ever changes.
#'
#' @param b_v_pr Square matrix of transition probabilities; b_v_pr[, c1] is
#'   the probability of transitioning to each class given previous class c1
#' @param c1 Which previous class (column) to compute log-odds for
#'
#' @return Numeric vector of length n_classes - 1
get_logratio_mean <- function(b_v_pr, c1) {
  p <- b_v_pr[, c1]
  p <- p / sum(p)
  log(p[-length(p)] / p[length(p)])
}

#' Prior means for the b_v reference-class log-odds parameterisation
#'
#' R port of compute_b_v_mean() in models/include/FoMo_functions.stan -
#' applies get_logratio_mean() to every (condition, previous-class) pair.
#' Returns a tidy table rather than the Stan version's array[K, n_classes]
#' vector[n_classes - 1], since R has no equivalent compact ragged-array
#' type; filter/join on `condition` and `prev_class` as needed.
#'
#' @param b_v_pr List of K transition-probability matrices, one per condition
#' @param n_classes Number of item classes
#'
#' @return Tibble with columns condition, prev_class, class (the "row" the
#'   log-odds is for, i.e. seq_len(n_classes - 1)), and mu
compute_b_v_mean <- function(b_v_pr, n_classes) {

  map(seq_along(b_v_pr), function(k) {
    map(seq_len(n_classes), function(c1) {
      mu <- get_logratio_mean(b_v_pr[[k]], c1)
      tibble(condition = k, prev_class = c1, class = seq_along(mu), mu = mu)
    }) |>
      list_rbind()
  }) |>
    list_rbind()
}
