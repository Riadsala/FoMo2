# =============================================================================
# FIT
#
# fit_fomo() fits the model a data list (prep_stan_list()) calls for; the
# model file is chosen from the data (see load_stan_model()).
# =============================================================================

fit_fomo <- function(sl, iter = 500, chains = 4, filename = "tmp",
                     save_dir = NULL, grainsize = 1, threads_per_chain = 2,
                     pathfinder_init = FALSE, seed = NULL, init = "prior", ...) {

  # the fit is saved as <save_dir>/<filename>.rds; save_dir defaults to
  # scratch/ at the repo root (gitignored), so fits don't land in whichever
  # folder the calling script runs from. seed = NULL: cmdstanr draws one
  # from R's generator, which set_up_FoMo.R seeds (FOMO_SEED), so the fit
  # still depends on the random draws made before it; set seed to fix it.

  mod <- load_stan_model(sl, model_task = "fit")

  # all fit models (single_level and multi_level) are reduce_sum-parallel
  # (see load_stan_model)
  sl$grainsize <- grainsize

  # pathfinder_init: run pathfinder() first and initialise the sampler
  # from it, instead of the starting values below. Off by default: on 5
  # hughes2024rsos participants (single-level cat model) it didn't shorten
  # warmup, so the whole fit was ~40% slower. Use it when a chain starts
  # badly (it fixed one stuck in early warmup on a real dataset);
  # fit_fomo() suggests it after a fit where one chain was much slower.
  #
  # init: "prior" (default) starts near the priors
  # (init_near_prior()). NULL is cmdstan's default: every parameter uniform
  # on (-2, 2) on its unconstrained scale, so random-effect SDs up to e^2
  # and person offsets up to +-2, which in the multi-level models starts far
  # from the posterior. A number is cmdstan's radius instead
  # of 2. A list or function goes straight to cmdstan. pathfinder_init =
  # TRUE overrides it. `...` goes to $sample() (e.g. save_warmup = TRUE).
  # Caveat: init_near_prior()'s jitter is small, so chains start close
  # together and R-hat, which relies on spread-out starts, may miss some
  # non-convergence
  if (identical(init, "prior") && !pathfinder_init) init <- init_near_prior(mod, sl, chains, seed)
  if (pathfinder_init) {
    pf <- mod$pathfinder(data = sl, num_threads = threads_per_chain, seed = seed)
    init <- pathfinder_inits(pf, mod, chains, seed)
  }

  m <- mod$sample(data = sl,
                  iter_warmup = iter,
                  iter_sampling = iter,
                  chains = chains,
                  threads_per_chain = threads_per_chain,
                  init = init,
                  seed = seed,
                  ...)

  # centring per-person effects (prep_stan_list(centred_params = ), the
  # default for multi-level fits) suits lots of data per person; with
  # little it makes a funnel, so point to the non-centred fit
  centred <- unlist(sl[startsWith(names(sl), "centre_")])
  n_divergent <- sum(m$diagnostic_summary(quiet = TRUE)$num_divergent)
  if (any(centred == 1) && n_divergent > 0) {
    message(n_divergent, " divergent transitions with centred person effects (",
            paste(sub("centre_", "", names(centred)[centred == 1]), collapse = ", "),
            "). With few trials per person, try prep_stan_list(..., centred_params = \"none\")")
  }

  # a chain much slower than the others usually started somewhere bad and
  # crawled through warmup (on one real-data fit, one chain was at warmup
  # iteration 100 when the others had finished); pathfinder initial
  # values fixed it
  chain_mins <- m$time()$chains$total / 60
  if (!pathfinder_init && length(chain_mins) > 1 && max(chain_mins) > 2 * median(chain_mins)) {
    message("Chain ", which.max(chain_mins), " took ", round(max(chain_mins)), " min, against a median of ",
            round(median(chain_mins)), " min: it probably started badly. Check its draws, ",
            "and try fit_fomo(..., pathfinder_init = TRUE)")
  }

  # the data list travels with the fit (attr(fit, "fomo_sl"), kept by
  # save_object()), so predict_fomo() and fomo_compare() can use it later
  attr(m, "fomo_sl") <- sl

  if (is.null(save_dir)) {
    save_dir <- file.path(rprojroot::find_root(rprojroot::has_dir("models")), "scratch")
  }
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  m$save_object(file.path(save_dir, paste0(filename, ".rds")))

  m
}

# A fit's draws of mod's parameters only (not its generated quantities,
# e.g. a fit's stored P/W/loglik: up to 77,000 columns), as a draws matrix
# for mod$generate_quantities() (predict_fomo(), sim_fomo()). Skips size-0
# parameters (no draws, and left out of cmdstanr's metadata). Only
# parameters are passed: every simplex in the models is a transformed
# parameter, recomputed exactly from them
parameter_draws <- function(fit, mod, sl) {
  sizes <- fit$metadata()$stan_variable_sizes
  pars <- keep(names(mod$variables()$parameters), \(p) prod(sizes[[p]] %||% 0) > 0)
  posterior::as_draws_matrix(fit$draws(variables = pars))
}

# f() with its standard output captured (cmdstanr's progress lines);
# messages, warnings and errors still show
run_quietly <- function(f) {
  out <- NULL
  utils::capture.output(out <- f())
  out
}

# One-step-ahead predictions for every observed selection, given the
# selections before it: the item sampled from the model's probabilities
# (P), the most probable item (W) and the log-likelihood of the actual
# choice (loglik), on n_draws of the posterior draws. For accuracy, and for
# fits made with prep_stan_list(save_selections = FALSE), which don't store
# them. (sim_fomo(sl, fit) differs: it simulates whole new selection
# sequences.) Runs the fit model's generated quantities with the flag on,
# reading only the fit's parameters (the model recomputes the rest), so it
# stays small however big the dataset. sl: the data list the fit was made
# with (by default the one fit_fomo() stored with it). Returns a
# CmdStanGQ object, which extract_post(pred, found, stim)
# reads like a fit.
predict_fomo <- function(fit, sl = attr(fit, "fomo_sl"), n_draws = 500, threads_per_chain = 2,
                         seed = NULL, selections = TRUE) {

  stopifnot("no data list: pass sl (fits from fit_fomo() carry their own)" = !is.null(sl))

  # selections = FALSE: only the per-trial log_lik (kfold_fomo())
  sl$save_selections <- as.integer(selections)
  sl$grainsize <- sl$grainsize %||% 1
  mod <- load_stan_model(sl, model_task = "fit", quiet = TRUE)
  draws <- parameter_draws(fit, mod, sl)
  # the seed only for choosing the draws: set.seed() here would reset R's
  # global random numbers, which cmdstanr also uses to name its output
  # files, so two calls within a minute could write (and later delete) the
  # same file
  pick <- function() sample(nrow(draws), min(n_draws, nrow(draws)))
  idx <- if (is.null(seed)) pick() else withr::with_seed(seed, pick())
  # cmdstanr prints its progress to standard output (not as messages):
  # captured, so it doesn't clutter pages and logs. Errors still show
  run_quietly(\() mod$generate_quantities(fitted_params = draws[idx, ], data = sl,
                                     threads_per_chain = threads_per_chain, seed = seed))
}

# The data entries of a data list (prep_stan_list()): everything sized by
# items, trials or selections. Everything else (K, L, the item classes or
# dimensions, priors, directional model, centring, ...) is a setting
STAN_LIST_DATA <- c("N", "T", "z", "s", "X", "x", "y", "v", "Y", "trial", "new_trl_flag", "Ns", "f")

# A data list for the data d (list(stim, found)) with the settings of sl, an
# existing data list: sl with its data entries replaced by d's. Used to
# refit or predict a model on part of its data (kfold_fomo()) without
# knowing the prep_stan_list() arguments that made it. d must hold the same
# people, numbered the same way, and conditions (K and L are kept): a person
# missing from d would leave the others' numbers pointing at the wrong
# person's effects
prep_stan_list_like <- function(sl, d) {
  stopifnot("sl has no data: make it with prep_stan_list(stim, found = )" = !is.null(sl$Y))
  new <- suppressMessages(prep_stan_list(d$stim, found = d$found, xy_scaling = sl$xy_scaling %||% 10))
  stopifnot("d has other conditions than sl" = new$K == sl$K,
            "d must hold every person of sl (people numbered 1..L as in sl)" =
              new$L == sl$L && n_distinct(d$stim$person) == sl$L)
  data <- intersect(STAN_LIST_DATA, names(new))   # v only with item features
  sl[data] <- new[data]
  sl
}

# K-fold cross-validation: the fallback for comparing models when loo's
# Pareto k diagnostics are poor (fomo_compare(kfold_cv = K) runs it for
# every model). Each stratum's trials (by default each person x condition)
# are dealt round the K folds in turn, so every fold holds every person and
# condition, and held-out trials are predicted with their person's own
# fitted effects (prediction of a known person's new trials). For each
# fold: fit_fomo() on the other folds (save_selections = FALSE), then
# predict_fomo() with every draw on the held-out trials, whose per-trial
# log_lik give each trial's elpd. So K folds cost K fits.
# d: list(stim, found) (import_data() or simulated). The model: sl, the
# data list of a fit to all of d (attr(fit, "fomo_sl")), or prep_args, the
# prep_stan_list() arguments (e.g. directional_model, harmonics,
# dir_offset). strata: "person" must be one, as held-out trials need their
# person's effects; c("person", "condition") (default) or "person". group:
# not yet supported (leaving out whole groups, e.g. people, needs
# predictions for new people, which the models don't make yet). quiet:
# no progress output. ...: fit_fomo() arguments (iter, chains,
# threads_per_chain, ...). Returns list(elpd, folds, fits): elpd a loo
# package elpd object, one value per trial in the order of
# sort(unique(d$stim$trial)), so loo::loo_compare() compares models fitted
# to the same data and folds; folds the fold of each trial; fits each
# fold's diagnostics.
kfold_fomo <- function(d, K = 3, sl = NULL, prep_args = list(), strata = c("person", "condition"),
                       group = NULL, quiet = FALSE, filename = "kfold", save_dir = NULL, ...) {

  if (!is.numeric(K) || length(K) != 1 || is.na(K) || K < 2 || K != round(K)) {
    stop("kfold_fomo(): K must be a whole number of at least 2 (K = 1 would leave no trials to ",
         "predict); got ", deparse(K), call. = FALSE)
  }
  if (!is.null(group)) {
    stop("kfold_fomo(): group (leaving out whole groups, e.g. people) isn't supported yet: it needs ",
         "predictions for new people, which the models don't make", call. = FALSE)
  }
  if (!("person" %in% strata) || !all(strata %in% c("person", "condition"))) {
    stop("kfold_fomo(): strata must be c(\"person\", \"condition\") or \"person\": held-out trials ",
         "are predicted with their person's effects, so every fold needs every person", call. = FALSE)
  }

  folds <- d$stim |>
    distinct(person, condition, trial) |>
    arrange(person, condition, trial) |>
    mutate(fold = (row_number() - 1) %% K + 1, .by = all_of(strata))
  short <- folds |> count(across(all_of(strata))) |> filter(n < K)
  if (nrow(short) > 0) {
    stop("kfold_fomo(): ", nrow(short), " ", paste(strata, collapse = " x "),
         " cells have fewer than K = ", K, " trials, so some fold would miss them; use a smaller K",
         call. = FALSE)
  }

  # each fold's data list: sl's settings (sl given) or prep_args
  prep <- function(x) {
    out <- if (!is.null(sl)) prep_stan_list_like(sl, x) else
      suppressMessages(do.call(prep_stan_list, c(list(x$stim, found = x$found), prep_args)))
    # only the held-out trials' log_lik is needed, not P, W and loglik
    out$save_selections <- 0L
    out
  }
  say <- function(...) if (!quiet) message(...)
  run <- function(f) if (quiet) suppressMessages(run_quietly(f)) else f()

  ll <- list(); fits <- list()
  for (k in seq_len(K)) {
    say("kfold_fomo(): fold ", k, " of ", K, " (", format(Sys.time(), "%H:%M"), ")")
    held_out <- filter(folds, fold == k)
    train <- subset_trials(d, filter(folds, fold != k))
    test  <- subset_trials(d, held_out)

    t0 <- Sys.time()
    fit <- run(\() fit_fomo(prep(train), filename = paste0(filename, "_fold", k), save_dir = save_dir, ...))
    pred <- run(\() predict_fomo(fit, prep(test), n_draws = .Machine$integer.max, selections = FALSE))

    # subset_trials() numbers the held-out trials in (person, condition,
    # trial) order, as held_out is sorted: map them back to the originals
    m <- unclass(pred$draws("log_lik", format = "matrix"))
    stopifnot("held-out log_lik doesn't match the held-out trials" = ncol(m) == nrow(held_out))
    colnames(m) <- held_out$trial
    ll[[k]] <- m

    diag <- fit$diagnostic_summary(quiet = TRUE)
    fits[[k]] <- tibble(fold = k, n_train_trials = n_distinct(train$stim$trial),
                        n_test_trials = nrow(held_out),
                        wall_mins = round(as.numeric(Sys.time() - t0, units = "mins"), 1),
                        max_rhat = round(max(fit$summary(fomo_variables(fit))$rhat, na.rm = TRUE), 3),
                        n_divergent = sum(diag$num_divergent))
  }

  # one column per trial, in trial order; each column's draws come from the
  # fit that held it out (every fit has the same number of draws)
  ll <- do.call(cbind, ll)
  ll <- ll[, order(as.numeric(colnames(ll))), drop = FALSE]
  list(elpd = loo::elpd(ll), folds = folds, fits = bind_rows(fits))
}

# one init list per chain, near the priors: population parameters at their
# prior means (sl's *_mu) with a little jitter, random-effect SDs small
# (0.1-0.3), class log-odds near 0
# (near-uniform transition probabilities). So every person starts close to the population
# value, as with weak data, and each chain gets its own jitter (R-hat
# needs dispersed starts). Works for any fit model: parameter shapes come
# from a zero-iteration run, and parameters not named here are left to
# cmdstan's random inits.
init_near_prior <- function(mod, sl, chains, seed = NULL) {

  sizes <- fomo_parameter_sizes(mod, sl)
  # seeded locally (see predict_fomo())
  if (!is.null(seed)) withr::local_seed(seed)

  jitter <- function(base, dims, sd) array(rep_len(base, prod(dims)) + rnorm(prod(dims), 0, sd), dim = dims)

  one_chain <- function() {
    out <- list()
    for (p in names(sizes)) {
      d <- sizes[[p]]
      out[[p]] <- switch(
        p,
        rho_delta = pmax(jitter(sl$rho_delta_mu, d, 0.1), 0.05),
        rho_psi = jitter(sl$rho_psi_mu, d, 0.1),
        log_theta = jitter(sl$log_theta_mu, d, 0.1),
        fourier_coef = jitter(sl$fourier_coef_mu, d, 0.05),
        # person values (multi-level models) start near the population's
        u_delta = pmax(jitter(sl$rho_delta_mu, d, 0.1), 0.05),
        u_psi = jitter(sl$rho_psi_mu, d, 0.1),
        u_log_theta = jitter(sl$log_theta_mu, d, 0.1),
        u_fourier_coef = jitter(sl$fourier_coef_mu, d, 0.05),
        u_logit_v = jitter(0, d, 0.1),
        # cat_sym.stan: the stay log-odds
        b_stay = jitter(0, d, 0.1),
        u_stay = jitter(0, d, 0.1),
        b_v_logodds = jitter(0, d, 0.1),
        eta_v = jitter(0, d, 0.1),
        if (str_starts(p, "sigma_")) array(runif(prod(d), 0.1, 0.3), dim = d)
        else if (str_detect(p, "^z")) jitter(0, d, 0.1)
        else NULL)
    }
    out <- out[!vapply(out, is.null, logical(1))]
    # scalars as plain numbers, not length-1 arrays
    scalar <- names(out)[vapply(names(out), \(p) mod$variables()$parameters[[p]]$dimensions == 0, logical(1))]
    out[scalar] <- lapply(out[scalar], as.numeric)
    out
  }

  lapply(seq_len(chains), function(i) one_chain())
}

# the declared size of every parameter of a fit model, from a one-iteration
# run with tree depth 1 (about two gradient evaluations, plus the model's
# transformed data), whose csv header cmdstanr reads the sizes from
fomo_parameter_sizes <- function(mod, sl) {
  sl$grainsize <- sl$grainsize %||% 1
  f <- mod$sample(data = sl, chains = 1, iter_warmup = 0, iter_sampling = 1, adapt_engaged = FALSE,
                  max_treedepth = 1, threads_per_chain = 1, refresh = 0,
                  show_messages = FALSE, show_exceptions = FALSE)
  # parameters of size 0 (e.g. a switched-off option's) are left out of
  # cmdstanr's metadata: skip them (they need no initial value)
  sizes <- f$metadata()$stan_variable_sizes
  sizes[intersect(names(mod$variables()$parameters), names(sizes))]
}

# one init list per chain from random pathfinder draws. Built by hand
# because cmdstanr's init = <pathfinder fit> drops size-1 dimensions (e.g.
# b_v with 2 classes is declared (K, 2, 1) but arrives as (K, 2), and every
# chain fails at initialisation). Stan's draw columns are column-major, as
# array() fills, so each parameter is reshaped to its declared sizes
pathfinder_inits <- function(pf, mod, chains, seed = NULL) {
  pars <- mod$variables()$parameters
  sizes <- pf$metadata()$stan_variable_sizes
  # parameters of size 0 have no draws and need no initial value
  params <- keep(names(pars), \(p) prod(sizes[[p]] %||% 0) > 0)
  draws <- posterior::as_draws_matrix(pf$draws(variables = params))
  # seeded locally (see predict_fomo())
  if (!is.null(seed)) withr::local_seed(seed)
  rows <- sample(nrow(draws), chains)
  lapply(rows, function(r) {
    setNames(lapply(params, function(p) {
      if (pars[[p]]$dimensions == 0) return(as.numeric(draws[r, p]))
      cols <- grep(paste0("^", p, "\\["), colnames(draws), value = TRUE)
      array(as.numeric(draws[r, cols]), dim = sizes[[p]])
    }), params)
  })
}

# quiet: no "Loading ... model" message (predict_fomo(), which runs inside
# fomo_compare() and kfold_fomo())
load_stan_model <- function(sl, model_task, quiet = FALSE) {

  root <- rprojroot::find_root(rprojroot::has_dir("models"))

  # chosen once, in prep_stan_list(directional_model = ...)
  directional_model <- get_directional_model(sl)

  # model file from the item features, and single- or multi-level from the
  # number of participants
  item_features <- if ("n_classes" %in% names(sl)) "cat" else "spatial"
  modfile <- paste0(item_features, ".stan")
  model_type <- if (length(unique(sl$z)) > 1) "multi_level" else "single_level"
  # the symmetric transition matrix (prep_stan_list(symmetric_v = TRUE)):
  # cat_sym.stan, a multi-level Fourier fit model only
  if (isTRUE(sl$symmetric_v == 1)) {
    if (model_task != "fit") stop("symmetric_v (cat_sym.stan) has no simulate model (so far)", call. = FALSE)
    modfile <- "cat_sym.stan"
  }

  if (!quiet) message("Loading ", model_type, " ", directional_model, " model: ", modfile)

  # fit models (single_level and multi_level) sum trial log-likelihoods via
  # reduce_sum and need STAN_THREADS=TRUE to compile; simulate models are
  # unchanged and don't need it. stanc's O1 optimizer gives a further ~10%
  # speedup on these same models (checked against baseline for the single_level
  # spatial fit model: identical log-lik, posterior means within
  # MCMC noise) with no downside, so it's on whenever we compile with
  # threading.
  is_threaded <- model_task == "fit"

  cmdstan_model(file.path(root, "models", model_type, model_task, directional_model, modfile),
                cpp_options = if (is_threaded) list(stan_threads = TRUE) else list(),
                stanc_options = if (is_threaded) list("O1") else list())
}

get_directional_model <- function(sl) {
  # prep_stan_list() adds n_harmonics for Fourier and n_dircomps for von Mises
  if ("n_harmonics" %in% names(sl)) return("fourier")
  if ("n_dircomps" %in% names(sl)) return("von_mises")
  stop("can't tell the directional model from sl - build it with prep_stan_list()")
}
