// FoMo multi-level cat model, Fourier directional model: fit.
//
// Each selection is a categorical choice among the items left, with
// log-weights from distance (rho_delta), absolute direction (fourier_coef:
// {cos, sin} per harmonic, or cos only) and relative direction (rho_psi), plus the
// item's class through the class-transition matrix b_v (row c: the next
// class after class c). A trial's first selection has class weights only
// (first_item_weight()): the stationary distribution of b_v.
//
// Multi-level: per-person effects on every term, each group centred or
// not (centre_*: prep_stan_list(centred_params)). Not centred: distance and
// class effects non-centred (offset/multiplier); relative and absolute
// direction in sum-to-zero form (see their parameters).
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
// Simulate model: ../../simulate/fourier/cat.stan (computes each selection
// separately: the readable reference for this likelihood).

functions {
  #include ../../../include/FoMo_functions.stan
  #include ../../../include/fourier_functions.stan
  #include ../../../include/flat_functions.stan

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
  // 1: the harmonic's sine term is fitted, 0: cosine only
  // (prep_stan_list(cos_only)); n_coef: the number of coefficients,
  // n_harmonics + sum(has_sin)
  array[n_harmonics] int<lower = 0, upper = 1> has_sin;
  int<lower = 1> n_coef;
  vector[K] dir_offset; // absolute-direction offset per condition (radians)

  // priors: absolute direction
  matrix[K, n_coef] fourier_coef_mu;
  vector[n_coef] fourier_coef_sd;

  // priors: between-person SDs
  real<lower = 0> sd_v_lambda;
  vector<lower = 0>[K] sd_delta_lambda;
  vector<lower = 0>[K] sd_psi_lambda;
  vector<lower = 0>[n_coef] sd_fourier_coef_lambda;

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

  int H2 = n_coef; // width of fourier_coef: cos (and sin) per harmonic
  if (n_coef != fourier_n_coef(has_sin)) reject("n_coef must be n_harmonics + sum(has_sin)");
  int D = H2 + 2; // rows of F: distance, Fourier basis, turning angle

  // start of each trial's items (cumulative s) and selections (cumulative f)
  array[T] int trial_start_idx = create_trial_start_idx(s, T);
  array[T] int sel_start_idx = create_trial_start_idx(f, T);

  // prior means of b_v_logodds: log-ratios within b_v_pr's columns
  array[K, n_classes] vector[n_classes - 1] b_v_mean;
  b_v_mean = compute_b_v_mean(b_v_pr, n_classes);

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
      compute_fourier_basis_by_cond(Ns, s, trial, X, phi, harmonics, has_sin, dir_offset);

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
  // relative direction (rho_psi) and absolute direction (fourier_coef:
  // {cos, sin} per harmonic, or cos only): with their per-person values below

  // class transitions: log-odds against the last class
  array[K, n_classes] vector[n_classes - 1] b_v_logodds;
  array[K, n_classes] vector<lower = 0>[n_classes - 1] sigma_u_v;

  // random effects

  // between-person SDs
  vector<lower = 0>[K] sigma_u_delta;
  vector<lower = 0>[K] sigma_u_psi;
  matrix<lower = 0>[K, H2] sigma_u_fourier_coef;

  // per-person values
  array[K, L, n_classes] vector<offset = person_offset(b_v_logodds, L, centre_v),
                                multiplier = person_multiplier(sigma_u_v, L, centre_v)>[n_classes - 1] u_logit_v;
  array[K] vector<offset = person_offset(to_vector(rho_delta), L, centre_delta),
                  multiplier = person_multiplier(sigma_u_delta, L, centre_delta)>[L] u_delta;

  // relative and absolute direction, population and per-person values
  // (rho_psi, u_psi, fourier_coef, u_fourier_coef: transformed parameters).
  // Centred (centre_psi / centre_direction = 1): sampled as they are
  // (centred_*). Otherwise, the default, in sum-to-zero form: for each
  // condition (and coefficient), with sigma the between-person SD,
  //   u_l = m + sigma * zt_l      mu = m - sigma * eta / sqrt(L)
  // zt a sum_to_zero_vector ~ normal(0, 1), eta ~ normal(0, 1). m is the
  // people's mean, which the data pin down, and eta how far the population
  // mean mu sits from it. The same model as non-centred u ~ normal(mu,
  // sigma), but without its ridge (mu trading off against the people's mean
  // standardised offset), so the population means mix far better (Nov
  // 2026); and unlike centring, no funnel when sigma is small. Unused sets
  // have size 0
  array[centre_psi ? K : 0] real centred_rho_psi;
  array[centre_psi ? K : 0] vector[L] centred_u_psi;
  array[centre_psi ? 0 : K] real m_psi;
  array[centre_psi ? 0 : K] real eta_psi;
  array[centre_psi ? 0 : K] sum_to_zero_vector[L] zt_psi;
  matrix[centre_direction ? K : 0, H2] centred_fourier_coef;
  array[centre_direction ? K : 0, H2] vector[L] centred_u_fourier_coef;
  matrix[centre_direction ? 0 : K, H2] m_fc;
  matrix[centre_direction ? 0 : K, H2] eta_fc;
  array[centre_direction ? 0 : K, H2] sum_to_zero_vector[L] zt_fc;

}

transformed parameters {
  // relative and absolute direction: population and per-person values
  // (see their parameters)
  array[K] real rho_psi;
  array[K] vector[L] u_psi;
  matrix[K, H2] fourier_coef;
  array[K, H2] vector[L] u_fourier_coef;
  for (k in 1:K) {
    if (centre_psi) {
      rho_psi[k] = centred_rho_psi[k];
      u_psi[k] = centred_u_psi[k];
    } else {
      rho_psi[k] = m_psi[k] - sigma_u_psi[k] * eta_psi[k] / sqrt(L);
      u_psi[k] = m_psi[k] + sigma_u_psi[k] * zt_psi[k];
    }
    for (d in 1:H2) {
      if (centre_direction) {
        fourier_coef[k, d] = centred_fourier_coef[k, d];
        u_fourier_coef[k, d] = centred_u_fourier_coef[k, d];
      } else {
        fourier_coef[k, d] = m_fc[k, d] - sigma_u_fourier_coef[k, d] * eta_fc[k, d] / sqrt(L);
        u_fourier_coef[k, d] = m_fc[k, d] + sigma_u_fourier_coef[k, d] * zt_fc[k, d];
      }
    }
  }

  // class-transition probabilities, as in the single-level models: row c1
  // is a simplex over the next class, after an item of class c1
  array[K, n_classes] simplex[n_classes] b_v = logodds_to_simplex(b_v_logodds);

  // per-person transition matrices
  array[K, L, n_classes] simplex[n_classes] u_v;
  u_v = class_probs(u_logit_v);

  // per-person first-selection class weights: the stationary distribution of
  // their transition matrix
  array[K, L] simplex[n_classes] u_first_v;
  for (k in 1:K)
    for (l in 1:L)
      u_first_v[k, l] = stationary_class_weights(u_v[k, l]);

}

model {
  // priors: transition matrices
  for (k in 1:K)
    for (c1 in 1:n_classes) {
      b_v_logodds[k, c1] ~ normal(b_v_mean[k, c1], b_v_sd);
      sigma_u_v[k, c1] ~ exponential(sd_v_lambda);
    }

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
  // relative and absolute direction: centred, or in sum-to-zero form (the
  // same model: the change of variables from (mu, u) has a constant Jacobian)
  if (centre_psi) {
    for (k in 1:K) u_psi[k] ~ normal(rho_psi[k], sigma_u_psi[k]);
  } else {
    eta_psi ~ std_normal();
    for (k in 1:K) zt_psi[k] ~ std_normal();
  }
  if (centre_direction) {
    for (k in 1:K)
      for (d in 1:H2)
        u_fourier_coef[k, d] ~ normal(fourier_coef[k, d], sigma_u_fourier_coef[k, d]);
  } else {
    to_vector(eta_fc) ~ std_normal();
    for (k in 1:K)
      for (d in 1:H2) zt_fc[k, d] ~ std_normal();
  }

  for (k in 1:K)
    for (l in 1:L)
      for (c1 in 1:n_classes)
        u_logit_v[k, l, c1] ~ normal(b_v_logodds[k, c1], sigma_u_v[k, c1]);

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

    for (c1 in 1:n_classes)
      pr_b_v[k, c1] =
        logodds_to_simplex(to_vector(normal_rng(b_v_mean[k, c1], b_v_sd)));
  }

  // for extract_post()
  vector[K] directional_offset = dir_offset;
  // label fourier_coef's columns (extract_fourier())
  array[n_harmonics] int fourier_harmonics = harmonics;
  array[n_harmonics] int fourier_has_sin = has_sin;

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
