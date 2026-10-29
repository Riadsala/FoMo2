# von Mises concentration of each directional component: the one default
# used for simulating, fitting and plotting
FOMO_KAPPA <- 30

# per-person effects sampled centred by default in multi-level fits (each
# model centres the ones it has; spatial models have no v). From an Oct
# 2026 comparison on the real datasets: centring delta and v halved the
# gradient evaluations with no divergences, and was the only setting that
# fitted one of them; centring direction or psi made funnels
FOMO_CENTRED_DEFAULT <- c("delta", "v")

prep_stan_list <- function(stim, sim_params = NULL, found = NULL, priors = NULL,
                           xy_scaling = 10, directional_model = "fourier",
                           kappa = FOMO_KAPPA, n_dircomps = 4, dir_offset = 0,
                           harmonics = c(1, 2, 4), centred_params = "default",
                           symmetric_v = FALSE, save_selections = FALSE) {

  # Builds the data list for fit_fomo() and sim_fomo(). Three uses:
  #   fitting:           prep_stan_list(stim, found = found)
  #   prior predictive:  prep_stan_list(stim), then sim_fomo(): simulates
  #                      from the fitting priors
  #   simulating from chosen values (as in examples/simulated/):
  #                      prep_stan_list(stim, sim_params), then sim_fomo()
  #
  # Fitting priors come from default_priors.csv. `priors` overrides entries
  # by their csv name, one value each, e.g. priors = list(rho_delta_sd = 0.5).
  # The csv has Fourier priors for harmonics 1, 2 and 4 only; any other
  # harmonic needs both of its priors here, e.g. harmonics = 1:4 with
  # priors = list(fourier_h3_sd = 0.4, fourier_h3_sd_lambda = 8).
  #
  # sim_params: the values to simulate from, named as in the Stan data
  # (see set_sim_params() below for what each one can be). Anything left
  # out is set to the fitting prior's value (after any `priors` overrides),
  # with a one-line message.
  # The simulate models still draw each fixed effect from normal(mu, sd),
  # so the *_sd entries default to 0.01, which pins it to mu.
  #
  # directional_model defaults to "fourier", with harmonics c(1, 2, 4),
  # Fourier-124: it beat von Mises (the papers' model) and Fourier-12 on
  # loo on every real dataset (examples/real_data/model_comparison/), and
  # fitted faster than von Mises.
  #
  # centred_params (multi-level models only): which groups of per-person
  # effects are sampled centred rather than non-centred: "default"
  # (FOMO_CENTRED_DEFAULT, delta and v, as far as the model has them),
  # "none", "all", or any of "delta", "psi", "direction" and "v" (cat
  # models), e.g. c("delta", "v"). Same model either way; only the
  # sampling geometry changes. Centring suits effects each person's data
  # pin down well, non-centring weakly informed ones (small simulated
  # datasets may sample better with "none").
  #
  # symmetric_v = TRUE (cat items, multi-level Fourier fits only, so far)
  # selects cat_sym.stan: each condition's transition matrix has one
  # parameter, the stay log-odds (b_stay), with every other class equally
  # likely, so the classes have no weights of their own. A simpler
  # alternative to the full matrix, for model comparison. No simulate model.
  #
  # save_selections (fitting only): FALSE (default)
  # leaves P, W and loglik out of the fit (size-0 arrays) and keeps only
  # the per-trial log_lik that loo needs; predict_fomo() gets them
  # afterwards from a subset of draws (accuracy: extract_post(pred, ...),
  # plot_post_accuracy()). TRUE makes the fit models output them for every
  # selection and draw, which is too big to read back into R on large
  # datasets (hughes2024rsos at iter 1000 needed more than 16 GB).

  stopifnot("directional_model must be 'von_mises' or 'fourier'" =
            directional_model %in% c("von_mises", "fourier"))
  # the data format the models rely on (data/README.md); import_data() also
  # prints notes on it
  verify_data(stim, if (is.data.frame(found)) found, "prep_stan_list()", show_notes = FALSE)
  if (!is.null(sim_params)) {
    stopifnot("sim_params is for simulating: to fit, pass found and no sim_params" =
              !is.data.frame(found))
  }

  ########################################################################
  # generic inputs for Stan

  K <- nlevels(stim$condition)
  stim$condition <- as.numeric(stim$condition)

  item_features <- if ("v" %in% names(stim)) "cat" else "none"

  stim |>
    group_by(trial) |>
    summarise(s = n(), person = unique(person), condition = unique(condition)) -> trials

  stanlist <- list(
    N = nrow(stim),
    K = K,
    T = nrow(trials),
    L = max(stim$person),
    z = trials$person,
    s = trials$s,
    X = trials$condition,
    x = stim$x * xy_scaling,
    y = stim$y * xy_scaling,
    # kept for prep_stan_list_like() (not read by the models)
    xy_scaling = xy_scaling)

  # categorical item features: classes per condition
  if (item_features == "cat") {
    n_classes <- stim |> summarise(n = n_distinct(v), .by = condition) |> pull(n) |> max()
    stanlist <- append(stanlist, list(v = stim$v, n_classes = n_classes))
  }

  # behaviour data, if present
  if (is.data.frame(found)) {
    stanlist$Y <- found$id
    stanlist$trial <- found$trial
    stanlist$new_trl_flag <- as.numeric(found$n)
    stanlist$Ns <- nrow(found)
    stanlist$f <- count(found, trial)$n
    stanlist$save_selections <- as.integer(save_selections)
  }

  ########################################################################
  # priors, or the values to simulate from

  dims <- list(K = K, n_dircomps = n_dircomps, n_coef = 2 * length(harmonics),
               n_classes = stanlist$n_classes,
               multi_level = length(unique(stim$person)) > 1)
  if (symmetric_v) {
    stopifnot("symmetric_v needs categorical items (stim$v a factor)" = item_features == "cat",
              "symmetric_v is multi-level only (so far)" = dims$multi_level,
              "symmetric_v is Fourier only (so far)" = directional_model == "fourier")
  }

  pr <- default_prior_list(read_default_priors(priors), dims, directional_model,
                           harmonics, item_features)
  if (!is.null(sim_params)) {
    pr <- set_sim_params(pr, sim_params, dims, directional_model, item_features)
  }
  stanlist <- append(stanlist, pr)

  ########################################################################
  # directional hyper-parameters. dir_offset can be a single value (used
  # for every condition) or one value per condition. NB: these must be
  # supplied to *both* the simulation and the fit call (found = df doesn't
  # carry them).

  dir_offset <- as.array(expand_vec(dir_offset, K, "dir_offset"))

  if (directional_model == "von_mises") {
    stanlist <- append(stanlist, list(kappa = kappa, n_dircomps = n_dircomps,
                                      dir_offset = dir_offset))
  } else {
    # n_harmonics alongside harmonics: Stan's data block needs an explicit
    # size for arrays (no bare array[]), so the fourier models declare
    # `array[n_harmonics] int harmonics`
    stanlist <- append(stanlist, list(n_harmonics = length(harmonics),
                                      harmonics = as.array(as.integer(harmonics)),
                                      dir_offset = dir_offset))
  }

  # the symmetric transition matrix (cat_sym.stan, chosen in load_stan_model())
  if (symmetric_v) stanlist$symmetric_v <- 1L

  ########################################################################
  # multi-level models: which per-person effects are centred

  if (length(unique(stim$person)) > 1) {
    stanlist <- append(stanlist, centring_flags(centred_params, item_features))
  } else if (!(identical(centred_params, "none") || identical(centred_params, "default"))) {
    stop("centred_params only applies to multi-level models (data with more than one person)")
  }

  return(stanlist)

}

centring_flags <- function(centred_params, item_features) {
  # centred_params ("default", "none", "all" or group names) -> the
  # multi-level fit models' data flags centre_<group> (1 = centred, 0 =
  # non-centred)
  groups <- c("delta", "psi", "direction", if (item_features == "cat") "v")

  centred <- if (identical(centred_params, "all")) groups
             else if (identical(centred_params, "default")) intersect(FOMO_CENTRED_DEFAULT, groups)
             else setdiff(centred_params, "none")
  unknown <- setdiff(centred, groups)
  if (length(unknown) > 0) {
    stop("centred_params: no group ", paste(unknown, collapse = ", "), " in this model; its groups are ",
         paste(groups, collapse = ", "), " (or \"default\" / \"none\" / \"all\")")
  }

  set_names(as.list(as.integer(groups %in% centred)), paste0("centre_", groups))
}

default_prior_list <- function(dp, dims, directional_model, harmonics, item_features) {
  # The fitting priors in Stan's shapes: default_priors.csv holds one base
  # value per prior (dp, from read_default_priors()), replicated here to K
  # conditions, n_dircomps directions, Fourier coefficients or item
  # classes, since those depend on the stimulus.

  K <- dims$K

  pr <- list(
    rho_delta_mu = rep(dp$rho_delta_mu, K),
    rho_delta_sd = rep(dp$rho_delta_sd, K),
    rho_psi_mu = rep(dp$rho_psi_mu, K),
    rho_psi_sd = rep(dp$rho_psi_sd, K),
    # random-effect scales: sigma ~ exponential(lambda)
    sd_delta_lambda = rep(dp$sd_delta_lambda, K),
    sd_psi_lambda = rep(dp$sd_psi_lambda, K))

  if (directional_model == "von_mises") {
    pr <- append(pr, list(
      log_theta_mu = matrix(dp$log_theta_mu, nrow = K, ncol = dims$n_dircomps),
      log_theta_sd = dp$log_theta_sd,
      sd_log_theta_lambda = rep(dp$sd_log_theta_lambda, dims$n_dircomps)))
  } else {
    # {cos, sin} coefficients per harmonic, mean 0 (no preferred direction)
    pr <- append(pr, list(
      fourier_coef_mu = matrix(0, nrow = K, ncol = dims$n_coef),
      fourier_coef_sd = get_fourier_coef_sd(harmonics, dp),
      sd_fourier_coef_lambda = get_fourier_coef_sd(harmonics, dp, "_sd_lambda")))
  }

  if (item_features == "cat") {
    # each row's log-odds (b_v_logodds, at both levels) ~
    # normal(log-ratios within each column of b_v_pr, b_v_sd), so equal
    # values give mean 0 and only b_v_sd sets the width
    mat <- matrix(dp$cat_b_v_pr, nrow = dims$n_classes, ncol = dims$n_classes)
    pr <- append(pr, list(b_v_pr = rep(list(mat), K),
                          b_v_sd = dp$cat_b_v_sd,
                          sd_v_lambda = dp$cat_sd_v_lambda))
    # (each trial's first selection has no class parameters of its own:
    # its class weights are the stationary distribution of b_v, see
    # stationary_class_weights() in FoMo_functions.stan)
  }

  pr
}

set_sim_params <- function(pr, sp, dims, directional_model, item_features) {
  # Replaces the fitting priors pr with the simulation values in sp
  # (prep_stan_list()'s sim_params). What each entry can be:
  #   rho_delta_mu, rho_psi_mu: one value, or one per condition
  #   log_theta_mu (von Mises) / fourier_coef_mu (Fourier):
  #     one value, one row (used for every condition) or a K-row matrix
  #   b_v_pr (cat): one class-transition matrix, or a list of K; column c1
  #     holds the probabilities of moving to each class after class c1
  #   sd_lambda: how much people differ from each other (multi-level models
  #     only). Each random-effect SD is drawn from exponential(sd_lambda), so
  #     a LARGER value means LESS person-to-person variation (mean SD =
  #     1 / sd_lambda). By parameter: rho_delta, rho_psi (one or K values),
  #     log_theta (one or n_dircomps), fourier_coef (one, one per harmonic
  #     or one per coefficient), b_v (cat: one value)
  #   *_sd (rho_delta_sd, log_theta_sd, ...): NOT between-person
  #     variation. The simulate model draws each population-level value from
  #     normal(*_mu, *_sd); the default 0.01 pins it to *_mu. There's rarely
  #     a reason to set these. They share names with the fitting priors
  #     (default_priors.csv) but don't affect the fit
  # Anything left out keeps the fitting prior's value, with a message.

  K <- dims$K

  dir_names <- if (directional_model == "von_mises") c("log_theta_mu", "log_theta_sd") else
    c("fourier_coef_mu", "fourier_coef_sd")
  item_names <- if (item_features == "cat") c("b_v_pr", "b_v_sd")
  allowed <- c("rho_delta_mu", "rho_delta_sd", "rho_psi_mu", "rho_psi_sd",
               dir_names, item_names, "sd_lambda")
  unknown <- setdiff(names(sp), allowed)
  if (length(unknown) > 0) {
    stop("sim_params entries not used by this model (", directional_model, ", ", item_features,
         " item features): ", paste(unknown, collapse = ", "))
  }

  # sp[[name]] if given, otherwise the default, with a message
  pick <- function(name, default, label = name) {
    if (!is.null(sp[[name]])) return(sp[[name]])
    message("sim_params: no ", label, ", using ", describe_default(default))
    default
  }
  pin <- function(name) if (is.null(sp[[name]])) 0.01 else sp[[name]]

  pr$rho_delta_mu <- expand_vec(pick("rho_delta_mu", pr$rho_delta_mu), K, "rho_delta_mu")
  pr$rho_delta_sd <- expand_vec(pin("rho_delta_sd"), K, "rho_delta_sd")
  pr$rho_psi_mu <- expand_vec(pick("rho_psi_mu", pr$rho_psi_mu), K, "rho_psi_mu")
  pr$rho_psi_sd <- expand_vec(pin("rho_psi_sd"), K, "rho_psi_sd")

  if (directional_model == "von_mises") {
    pr$log_theta_mu <- expand_mat(pick("log_theta_mu", pr$log_theta_mu), K, dims$n_dircomps, "log_theta_mu")
    pr$log_theta_sd <- pin("log_theta_sd")
  } else {
    pr$fourier_coef_mu <- expand_mat(pick("fourier_coef_mu", pr$fourier_coef_mu), K, dims$n_coef, "fourier_coef_mu")
    pr$fourier_coef_sd <- expand_vec(pin("fourier_coef_sd"), dims$n_coef, "fourier_coef_sd")
  }

  if (item_features == "cat") {
    b_v_pr <- pick("b_v_pr", pr$b_v_pr)
    if (is.matrix(b_v_pr)) b_v_pr <- rep(list(b_v_pr), K)
    stopifnot("b_v_pr must be one matrix or a list of one per condition" = length(b_v_pr) == K)
    # the models only use the ratios within each column of b_v_pr (as
    # log-odds prior means), and the pinned b_v_sd fixes the simulated
    # transition probabilities to them (the x 1000 changes nothing)
    pr$b_v_pr <- map(b_v_pr, ~ .x * 1000)
    pr$b_v_sd <- pin("b_v_sd")
    # zeros would give log(0)
    stopifnot("b_v_pr must be positive" = all(unlist(pr$b_v_pr) > 0))
  }

  # random-effect scales: only multi-level models use them, so only those
  # get messages about defaults
  sl <- sp$sd_lambda
  lam <- function(field, default) {
    if (!is.null(sl[[field]])) return(sl[[field]])
    if (dims$multi_level) message("sim_params: no sd_lambda$", field, ", using ", describe_default(default))
    default
  }

  pr$sd_delta_lambda <- expand_vec(lam("rho_delta", pr$sd_delta_lambda), K, "sd_lambda$rho_delta")
  pr$sd_psi_lambda <- expand_vec(lam("rho_psi", pr$sd_psi_lambda), K, "sd_lambda$rho_psi")
  if (directional_model == "von_mises") {
    pr$sd_log_theta_lambda <- expand_vec(lam("log_theta", pr$sd_log_theta_lambda), dims$n_dircomps,
                                         "sd_lambda$log_theta")
  } else {
    fc <- lam("fourier_coef", pr$sd_fourier_coef_lambda)
    if (length(fc) == dims$n_coef / 2) fc <- rep(fc, each = 2)  # one per harmonic
    pr$sd_fourier_coef_lambda <- expand_vec(fc, dims$n_coef, "sd_lambda$fourier_coef")
  }
  if (item_features == "cat") {
    pr$sd_v_lambda <- lam("b_v", pr$sd_v_lambda)
  }

  pr
}

describe_default <- function(x) {
  # short description of a default value for set_sim_params()'s messages
  if (is.list(x)) return("equal class-transition probabilities")
  vals <- unique(round(as.numeric(x), 3))
  if (length(vals) == 1) as.character(vals) else paste(vals, collapse = ", ")
}

expand_vec <- function(x, n, name) {
  # a single value, or exactly n values, as a length-n numeric vector
  x <- as.numeric(unlist(x))
  if (length(x) == 1) x <- rep(x, n)
  if (length(x) != n) stop(name, " must have length 1 or ", n)
  x
}

expand_mat <- function(x, nr, nc, name) {
  # a single value, one row of nc values (used for every row) or an
  # nr x nc matrix, as an nr x nc matrix
  if (length(x) == 1) return(matrix(as.numeric(x), nr, nc))
  if (is.matrix(x) && all(dim(x) == c(nr, nc))) return(x)
  if (length(x) == nc && (!is.matrix(x) || nrow(x) == 1)) {
    return(matrix(as.numeric(x), nr, nc, byrow = TRUE))
  }
  stop(name, " must be a single value, one row of ", nc, " values, or a ", nr, " x ", nc, " matrix")
}

get_fourier_coef_sd <- function(harmonics, dp, suffix = "_sd") {
  # Looks up each requested harmonic's default prior value
  # (fourier_h<k><suffix> in default_priors.csv: "_sd" for the coefficient
  # SD, "_sd_lambda" for the random-effect scale) and expands it to a
  # length-2*length(harmonics) vector - cos and sin share the same value
  # within a harmonic, since a normal(0, sigma) prior on both is what gives
  # the rotation-invariant (phase-free) amplitude distribution.
  fields <- paste0("fourier_h", harmonics, suffix)
  missing <- harmonics[!(fields %in% names(dp))]
  if (length(missing) > 0) {
    stop("no prior for harmonic ", paste(missing, collapse = ", "),
         " (default_priors.csv has 1, 2 and 4): pass it in priors, e.g. priors = list(",
         "fourier_h", missing[1], "_sd = 0.4, fourier_h", missing[1], "_sd_lambda = 8)")
  }
  rep(as.numeric(unlist(dp[fields])), each = 2)
}

read_default_priors <- function(overrides = NULL) {
  # Reads functions/default_priors.csv (the repo-wide fitting priors; see
  # that file's header for what it holds and why) into a named list, so
  # callers can do dp$rho_delta_mu etc. overrides: a named list of single
  # values replacing csv entries (prep_stan_list()'s `priors`), or adding
  # Fourier priors for harmonics the csv doesn't have (fourier_h<k>_sd,
  # fourier_h<k>_sd_lambda). Uses rprojroot to find the file from any
  # working directory.

  root <- rprojroot::find_root(rprojroot::has_dir("models"))
  dp <- read_csv(file.path(root, "functions", "default_priors.csv"),
                 comment = "#", show_col_types = FALSE)
  dp <- setNames(as.list(dp$value), dp$parameter)

  if (!is.null(overrides)) {
    new_harmonic <- str_detect(names(overrides), "^fourier_h[0-9]+_sd(_lambda)?$")
    unknown <- setdiff(names(overrides)[!new_harmonic], names(dp))
    if (length(unknown) > 0) {
      stop("unknown prior(s): ", paste(unknown, collapse = ", "),
           " - use the names in functions/default_priors.csv")
    }
    stopifnot("each prior in `priors` must be a single value" = all(lengths(overrides) == 1))
    stopifnot("Fourier prior SDs and lambdas must be positive" =
              all(unlist(overrides[new_harmonic]) > 0))
    for (nm in intersect(names(overrides), "cat_b_v_pr")) {
      if (overrides[[nm]] <= 0) stop(nm, " must be positive")
    }
    dp[names(overrides)] <- overrides
  }

  dp
}
