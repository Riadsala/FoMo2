// FoMo multi-level cat model with a symmetric transition matrix, Fourier
// directional model: fit. A copy of cat.stan (same directory), changed only
// in the class transitions: keep the two in step.
//
// As cat.stan, but each condition's transition matrix has one parameter,
// the stay log-odds b_stay: the log-odds of selecting an item of the class
// just selected against an item of any one other class. Every other class
// is equally likely, so for 2 classes the matrix is symmetric
// ([p, 1 - p; 1 - p, p]) and the classes have no weights of their own; the
// first selection's class weights (the stationary distribution) are then
// equal. Per person: u_stay (centred with centre_v). The transition matrix
// is output as in cat.stan (b_v_logodds: log-odds against the last class,
// b_v: probabilities; u_logit_v, u_v, u_first_v per person), so
// extract_post() and the plots read it unchanged. Selected by
// prep_stan_list(symmetric_v = TRUE). No simulate model.
//
// cat.stan's description, which applies here too:
//
// Each selection is a categorical choice among the items left, with
// log-weights from distance (rho_delta), absolute direction (fourier_coef:
// {cos, sin} pairs per harmonic) and relative direction (rho_psi), plus the
// item's class through the class-transition matrix b_v (row c: the next
// class after class c). A trial's first selection has class weights only
// (first_item_weight()): the stationary distribution of b_v.
//
// Multi-level: per-person effects on every term, each group centred or
// non-centred (centre_*: prep_stan_list(centred_params)).
//
// The likelihood uses the flat layout (../../../include/flat_functions.stan):
// every selection after a trial's first has one column per item left in F,
// a matrix of data with rows distance, Fourier basis (H2 rows) and turning
// angle, so an item's log-weight is
//   beta . F[, item] + log b_v[previous class, item's class]
// with beta = (-rho_delta, fourier_coef, -rho_psi) for that person and
// condition: one matrix product per trial, plus a lookup for the classes.
//
// Data-only quantities (F and its indices) are computed once in transformed
// data; trials are summed with reduce_sum. Generated quantities: prior
// draws, per-trial log_lik (loo) and, with save_selections = 1, per
// selection P (sampled item), W (most likely item) and loglik.
// The simulate model for cat.stan (../../simulate/fourier/cat.stan)
// computes each selection separately: the readable reference for this
// likelihood.

functions {
  #include ../../../include/FoMo_functions.stan
  #include ../../../include/fourier_functions.stan
  #include ../../../include/flat_functions.stan

  vector stay_to_logodds(real stay, int c1, int n) {
    // row c1 of a symmetric transition matrix from its stay log-odds, as
    // log-odds against the last class (cat.stan's scale): every other
    // class has log-weight 0, class c1 has stay
    vector[n] w = rep_vector(0, n);
    w[c1] = stay;
    return w[1:(n - 1)] - w[n];
  }

  // a person's coefficients, in F's row order: distance, Fourier basis,
  // turning angle (distance and turn enter with a minus sign)
  row_vector person_beta(array[] vector u_delta, array[,] vector u_fourier_coef,
                         array[] vector u_psi, int k, int l) {
    return append_col(append_col(-u_delta[k, l], get_person_fourier_coef(u_fourier_coef, k, l)),
                      -u_psi[k, l]);
  }

  // unnormalised log-weights of trial t's block of F: every item left at
  // each selection after the first. lv: the person's flattened log
  // transition matrix (flat_log_v()). Only for trials with f[t] > 1
  row_vector trial_weights(int t, array[] int flat_start, array[] int flat_len,
                           array[] int cls_flat, vector lv, row_vector beta, matrix F) {
    int st = flat_start[t];
    int len = flat_len[t];
    return beta * block(F, 1, st, rows(F), len)
           + to_row_vector(lv[cls_flat[st:(st + len - 1)]]);
  }

  // log-likelihood of a chunk of trials (reduce_sum)
  real partial_sum_lpmf(
      array[] int trial_slice, int start, int end,
      array[] int f, array[] int s, array[] int X, array[] int z,
      array[] int sel_start_idx, array[] int p_at_step,
      array[] int flat_start, array[] int flat_len, array[] int cls_flat,
      array[] int seg_rel, array[] int chosen_rel,
      array[] vector u_delta, array[,] vector u_fourier_coef, array[] vector u_psi,
      array[,,] vector u_v, array[,] vector u_first_v,
      matrix F, array[,] int vt_first) {
    real lp = 0;

    for (t in trial_slice) {
      int k = X[t];   // condition
      int l = z[t];   // person
      int idx = sel_start_idx[t];

      // first selection: class weights only (not in F)
      lp += first_item_weight(s[t], u_first_v[k, l], vt_first[t, 1:s[t]])[p_at_step[idx]];

      // the rest: every selection's log-weights at once
      if (f[t] > 1) {
        row_vector[flat_len[t]] w = trial_weights(t, flat_start, flat_len, cls_flat,
          flat_log_v(u_v[k, l]), person_beta(u_delta, u_fourier_coef, u_psi, k, l), F);
        lp += flat_trial_loglik(w, idx, s[t], f[t], 2, seg_rel, chosen_rel);
      }
    }

    return lp;
  }

}

data {
  int <lower = 1> N; // items, over all trials
  int <lower = 1> K; // conditions
  int <lower = 1> T; // trials
  int <lower = 1> L; // people

  // person of each trial
  array[T] int<lower = 1, upper = L> z;

  int <lower = 1, upper = N> Ns; // selections, over all trials
  array[Ns] int<lower = 1> Y; // selected item (index within its trial)
  array[Ns] int<lower = 1, upper = T> trial; // trial of each selection
  array[T] int<lower = 1> f; // selections per trial

  array[Ns] int<lower = 0> new_trl_flag; // selection's position in its trial (not a 0/1 flag)
  array[T] int<lower = 1> s; // items per trial

  // condition
  array[T] int<lower = 1, upper = K> X;

  // item positions
  vector[N] x;
  vector[N] y;

  // item classes
  int<lower = 1> n_classes;
  array[N] int<lower = 1, upper = n_classes> v;

  // priors: transition matrices
  array[K] matrix[n_classes, n_classes] b_v_pr;
  real<lower = 0> b_v_sd;
  // priors: distance
  vector<lower = 0>[K] rho_delta_mu;
  vector<lower = 0>[K] rho_delta_sd;
  // priors: relative direction
  vector[K] rho_psi_mu;
  vector<lower = 0>[K] rho_psi_sd;

  // harmonics to fit, e.g. {1, 2, 4} (n_harmonics: Stan needs data array
  // sizes)
  int<lower = 1> n_harmonics;
  array[n_harmonics] int harmonics;
  vector[K] dir_offset; // absolute-direction offset per condition (radians)

  // priors: absolute direction
  matrix[K, 2 * n_harmonics] fourier_coef_mu;
  vector[2 * n_harmonics] fourier_coef_sd;

  // priors: between-person SDs
  real<lower = 0> sd_v_lambda;
  vector<lower = 0>[K] sd_delta_lambda;
  vector<lower = 0>[K] sd_psi_lambda;
  vector<lower = 0>[2 * n_harmonics] sd_fourier_coef_lambda;

  // per-person effects centred (1) or not (0): prep_stan_list(centred_params)
  int<lower = 0, upper = 1> centre_delta;
  int<lower = 0, upper = 1> centre_psi;
  int<lower = 0, upper = 1> centre_direction;
  int<lower = 0, upper = 1> centre_v;

  int<lower = 1> grainsize; // reduce_sum grainsize (1: let TBB decide)

  int<lower=0, upper=1> save_selections; // 1: output P, W and loglik per selection

}

transformed data {
  // P, W and loglik have size 0 unless save_selections
  int Ns_out = save_selections ? Ns : 0;

  int H2 = 2 * n_harmonics; // width of fourier_coef (2 per harmonic)
  int D = H2 + 2; // rows of F: distance, Fourier basis, turning angle

  // start of each trial's items (cumulative s) and selections (cumulative f)
  array[T] int trial_start_idx = create_trial_start_idx(s, T);
  array[T] int sel_start_idx = create_trial_start_idx(f, T);

  // prior mean of b_stay: the stay log-odds implied by b_v_pr, averaged
  // over the previous classes. As in cat.stan, b_v_pr's column c1 (scaled to
  // sum to 1) is the prior for the next class after c1; the stay log-odds
  // is log p_c1 against the mean of the other entries
  vector[K] b_stay_mean = rep_vector(0, K);
  for (k in 1:K)
    for (c1 in 1:n_classes) {
      vector[n_classes] p = col(b_v_pr[k], c1) / sum(col(b_v_pr[k], c1));
      b_stay_mean[k] += (log(p[c1]) - log((1 - p[c1]) / (n_classes - 1))) / n_classes;
    }

  // items left at each selection, and the selected item's position among them
  array[Ns, max(s)] int avail_order = compute_availability_order(Ns, T, s, f, Y);
  array[Ns] int p_at_step = compute_selected_position(Ns, T, s, f, Y);

  // the flat layout (flat_functions.stan): every selection after a trial's
  // first (first_sel = 2)
  int M = count_flat(s, f, T, 2);
  array[T] int flat_start = flat_trial_start(s, f, T, 2);
  array[T] int flat_len = flat_trial_len(s, f, T, 2);
  array[Ns] int seg_rel = flat_seg_rel(s, f, T, Ns, 2);
  array[Ns] int chosen_rel = flat_chosen_rel(seg_rel, p_at_step);
  matrix[D, M] F;          // each column's features
  array[M] int cls_flat;   // each column's index into flat_log_v()

  // each trial's item classes in availability order at its first selection
  // (unfilled beyond s[t])
  array[T, max(s)] int vt_first;

  {
    // distance, direction and turning angle by item id: local to this
    // block, so they are freed once copied into F
    matrix[Ns, max(s)] delta = compute_delta(Ns, Y, s, trial, x, y, new_trl_flag, trial_start_idx);
    matrix[Ns, max(s)] phi   = compute_phi(Ns, Y, s, trial, x, y, new_trl_flag, trial_start_idx);
    matrix[Ns, max(s)] psi   = compute_psi(Ns, Y, s, trial, x, y, new_trl_flag, trial_start_idx);
    // Fourier basis for every selection, once
    array[Ns] matrix[H2, max(s)] basis =
      compute_fourier_basis_by_cond(Ns, s, trial, X, phi, harmonics, dir_offset);

    // each column's selection and item
    array[M, 2] int cols = flat_cols(s, f, T, 2, avail_order, M);
    for (c in 1:M) {
      int sel = cols[c, 1];
      int item = cols[c, 2];
      F[1, c] = delta[sel, item];
      F[2:(H2 + 1), c] = basis[sel][, item];
      F[D, c] = psi[sel, item];
    }
    cls_flat = flat_class_index(cols, Y, v, trial, trial_start_idx, n_classes);
  }

  for (t in 1:T) {
    int v0 = trial_start_idx[t] - 1;
    for (j in 1:s[t]) vt_first[t, j] = v[v0 + avail_order[sel_start_idx[t], j]];
  }

  // trial ids, for reduce_sum to partition
  array[T] int trial_ids;
  for (t in 1:T) trial_ids[t] = t;

}

parameters {
  // fixed effects
  array[K] real<lower = 0> rho_delta; // distance
  array[K] real rho_psi; // relative direction
  matrix[K, H2] fourier_coef; // absolute direction: {cos, sin} per harmonic

  // class transitions: the stay log-odds, and its between-person SD
  array[K] real b_stay;
  vector<lower = 0>[K] sigma_u_v;

  // random effects

  // between-person SDs
  vector<lower = 0>[K] sigma_u_delta;
  vector<lower = 0>[K] sigma_u_psi;
  matrix<lower = 0>[K, H2] sigma_u_fourier_coef;

  // per-person values
  array[K] vector<offset = person_offset(to_vector(b_stay), L, centre_v),
                  multiplier = person_multiplier(sigma_u_v, L, centre_v)>[L] u_stay;
  array[K] vector<offset = person_offset(to_vector(rho_delta), L, centre_delta),
                  multiplier = person_multiplier(sigma_u_delta, L, centre_delta)>[L] u_delta;
  array[K] vector<offset = person_offset(to_vector(rho_psi), L, centre_psi),
                  multiplier = person_multiplier(sigma_u_psi, L, centre_psi)>[L] u_psi;
  array[K, H2] vector<offset = person_offset(fourier_coef, L, centre_direction),
                      multiplier = person_multiplier(sigma_u_fourier_coef, L, centre_direction)>[L] u_fourier_coef;

}

transformed parameters {

  // the transition matrices on cat.stan's scales: log-odds against the last
  // class, for the population (b_v_logodds) and per person (u_logit_v), and
  // the population's as probabilities (b_v)
  array[K, n_classes] vector[n_classes - 1] b_v_logodds;
  array[K, L, n_classes] vector[n_classes - 1] u_logit_v;
  for (k in 1:K)
    for (c1 in 1:n_classes) {
      b_v_logodds[k, c1] = stay_to_logodds(b_stay[k], c1, n_classes);
      for (l in 1:L) u_logit_v[k, l, c1] = stay_to_logodds(u_stay[k, l], c1, n_classes);
    }
  array[K, n_classes] simplex[n_classes] b_v = logodds_to_simplex(b_v_logodds);

  // per-person transition matrices
  array[K, L, n_classes] simplex[n_classes] u_v;
  u_v = class_probs(u_logit_v);

  // per-person first-selection class weights: the stationary distribution of
  // their transition matrix, which for a symmetric matrix is uniform (kept,
  // as in cat.stan, for extract_post() and the likelihood)
  array[K, L] simplex[n_classes] u_first_v;
  for (k in 1:K)
    for (l in 1:L)
      u_first_v[k, l] = rep_vector(inv(n_classes), n_classes);

}

model {
  // priors: transition matrices
  b_stay ~ normal(b_stay_mean, b_v_sd);
  sigma_u_v ~ exponential(sd_v_lambda);

  rho_delta ~ normal(rho_delta_mu, rho_delta_sd);
  rho_psi   ~ normal(rho_psi_mu,   rho_psi_sd);

  for (k in 1:K) {
    to_vector(fourier_coef[k,]) ~ normal(to_vector(fourier_coef_mu[k,]), fourier_coef_sd);
  }

  // priors: between-person SDs
  sigma_u_delta ~ exponential(sd_delta_lambda);
  sigma_u_psi   ~ exponential(sd_psi_lambda);

  for (k in 1:K) {
    for (d in 1:H2) {
      sigma_u_fourier_coef[k, d] ~ exponential(sd_fourier_coef_lambda[d]);
    }
  }

  // per-person values around the population values (centred or not: see their
  // offset/multiplier)
  for (k in 1:K) u_delta[k] ~ normal(rho_delta[k], sigma_u_delta[k]);
  for (k in 1:K) u_psi[k] ~ normal(rho_psi[k], sigma_u_psi[k]);
  for (k in 1:K)
    for (d in 1:H2)
      u_fourier_coef[k, d] ~ normal(fourier_coef[k, d], sigma_u_fourier_coef[k, d]);

  for (k in 1:K) u_stay[k] ~ normal(b_stay[k], sigma_u_v[k]);

  // likelihood, summed over trials in parallel (reduce_sum)

  target += reduce_sum(partial_sum_lpmf, trial_ids, grainsize,
    f, s, X, z, sel_start_idx, p_at_step,
    flat_start, flat_len, cls_flat, seg_rel, chosen_rel,
    u_delta, u_fourier_coef, u_psi, u_v, u_first_v, F, vt_first);

}

generated quantities {
  // prior draws

  array[K] real pr_rho_delta, pr_rho_psi;
  matrix[K, H2] pr_fourier_coef;
  matrix[K, H2] pr_u_fourier_coef;
  array[K, n_classes] simplex[n_classes] pr_b_v;
  array[K] real pr_u_delta, pr_u_psi, pr_u_v;

  for (k in 1:K) {
    // fixed effects
    pr_rho_delta[k] = normal_lower0_rng(rho_delta_mu[k], rho_delta_sd[k]);
    pr_rho_psi[k]   = normal_rng(rho_psi_mu[k], rho_psi_sd[k]);
    // between-person SDs
    pr_u_delta[k] = exponential_rng(sd_delta_lambda[k]);
    pr_u_psi[k] = exponential_rng(sd_psi_lambda[k]);
    pr_u_v[k] = exponential_rng(sd_v_lambda);
    for (d in 1:H2)
      pr_u_fourier_coef[k, d] = exponential_rng(sd_fourier_coef_lambda[d]);

    for (d in 1:H2)
      pr_fourier_coef[k, d] = normal_rng(fourier_coef_mu[k, d], fourier_coef_sd[d]);

    {
      real pr_stay = normal_rng(b_stay_mean[k], b_v_sd);
      for (c1 in 1:n_classes) pr_b_v[k, c1] = logodds_to_simplex(stay_to_logodds(pr_stay, c1, n_classes));
    }
  }

  // for extract_post()
  vector[K] directional_offset = dir_offset;
  array[n_harmonics] int fourier_harmonics = harmonics; // labels fourier_coef's pairs (extract_fourier())

  // log-likelihood and predictions

  array[Ns_out] int P, W;
  vector[Ns_out] loglik;
  vector[T] log_lik; // per-trial log-likelihood (loo)

  loglik = rep_vector(0, Ns_out);
  log_lik = rep_vector(0, T);

  for (t in 1:T) {
    int k = X[t]; // condition
    int l = z[t]; // person
    int idx = sel_start_idx[t];

    // log-weights of the trial's later selections (unnormalised), once
    row_vector[flat_len[t]] w;
    if (f[t] > 1)
      w = trial_weights(t, flat_start, flat_len, cls_flat, flat_log_v(u_v[k, l]),
                        person_beta(u_delta, u_fourier_coef, u_psi, k, l), F);

    for (ii in 1:f[t]) {
      int n_avail = s[t] - ii + 1; // items left
      int p = p_at_step[idx];
      row_vector[n_avail] weights = ii == 1
        ? first_item_weight(n_avail, u_first_v[k, l], vt_first[t, 1:n_avail])
        : flat_selection_weights(w, seg_rel[idx], n_avail);

      if (save_selections) loglik[idx] = weights[p];
      log_lik[t] += weights[p];
      // P: an item sampled from the weights; W: the most likely item
      // (indices within the trial)
      if (save_selections) P[idx] = avail_order[idx, categorical_rng(exp(to_vector(weights)))];
      if (save_selections) W[idx] = avail_order[idx, which_max(weights)];

      idx += 1;
    }
  }
}
