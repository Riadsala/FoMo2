// FoMo single-level spatial model, Fourier directional model: fit.
//
// Each selection is a categorical choice among the items left, with
// log-weights from distance (rho_delta), absolute direction (fourier_coef:
// {cos, sin} pairs per harmonic) and relative direction (rho_psi).
//
// Data-only quantities (distances, directions, the items left) are computed
// once in transformed data; trials are summed with reduce_sum. Generated
// quantities: prior draws, per-trial log_lik (loo) and, with save_selections
// = 1, per selection P (sampled item), W (most likely item) and loglik.
// Simulate model: ../../simulate/fourier/spatial.stan.

functions {
  #include ../../../include/FoMo_functions.stan
  #include ../../../include/fourier_functions.stan

  row_vector compute_weights(
    int n_avail,
    real rho_delta, row_vector coef, real rho_psi,
    row_vector d, matrix basis, row_vector psi) {
      row_vector[n_avail] weights;

      weights = to_row_vector(-rho_delta * d);
      weights += compute_fourier_absdir_weights(coef, basis);
      weights += -rho_psi * psi;

      // every item passed in is available: normalise
      weights -= log_sum_exp(weights);

      return(weights);
  }

  // log-likelihood of a chunk of trials (reduce_sum); sel_start_idx[t]: trial
  // t's first row in the per-selection arrays
  real partial_sum_lpmf(
      array[] int trial_slice, int start, int end,
      array[] int f, array[] int s, array[] int X, array[] int sel_start_idx,
      array[] int p_at_step,
      array[] real rho_delta, matrix fourier_coef, array[] real rho_psi,
      matrix delta_g, array[] matrix basis_g, matrix psi_g) {
    real lp = 0;

    for (t in trial_slice) {
      int k = X[t];
      int n_avail = s[t];
      int idx = sel_start_idx[t];

      for (ii in 1:f[t]) {
        row_vector[n_avail] weights = compute_weights(n_avail,
          rho_delta[k], fourier_coef[k,], rho_psi[k],
          delta_g[idx, 1:n_avail], basis_g[idx][, 1:n_avail], psi_g[idx, 1:n_avail]);

        lp += weights[p_at_step[idx]];

        n_avail -= 1;
        idx += 1;

      }
    }

    return lp;
  }
}

data {
  int <lower = 1> N; // items, over all trials
  int <lower = 1> K; // conditions
  int <lower = 1> T; // trials

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

  // harmonics to fit, e.g. {1, 2, 4} (n_harmonics: Stan needs data array
  // sizes)
  int<lower = 1> n_harmonics;
  array[n_harmonics] int harmonics;
  vector[K] dir_offset; // absolute-direction offset per condition (radians)

  // priors: distance
  vector<lower = 0>[K] rho_delta_mu;
  vector<lower = 0>[K] rho_delta_sd;
  // priors: relative direction
  vector[K] rho_psi_mu;
  vector<lower = 0>[K] rho_psi_sd;
  // priors: absolute direction
  matrix[K, 2 * n_harmonics] fourier_coef_mu;
  vector[2 * n_harmonics] fourier_coef_sd;

  int<lower = 1> grainsize; // reduce_sum grainsize (1: let TBB decide)

  int<lower=0, upper=1> save_selections; // 1: output P, W and loglik per selection
}

transformed data {
  // P, W and loglik have size 0 unless save_selections
  int Ns_out = save_selections ? Ns : 0;

  int H2 = 2 * n_harmonics; // width of fourier_coef (2 per harmonic)

  // start of each trial's items (cumulative s) and selections (cumulative f)
  array[T] int trial_start_idx = create_trial_start_idx(s, T);
  array[T] int sel_start_idx = create_trial_start_idx(f, T);

  // items left at each selection, and the selected item's position among them
  array[Ns, max(s)] int avail_order = compute_availability_order(Ns, T, s, f, Y);
  array[Ns] int p_at_step = compute_selected_position(Ns, T, s, f, Y);

  // per selection, delta, direction and psi in availability order (unfilled
  // beyond n_avail)
  matrix[Ns, max(s)] delta_g;
  array[Ns] matrix[H2, max(s)] basis_g;
  matrix[Ns, max(s)] psi_g;

  {
    // distance, direction and turning angle by item id: local to this
    // block, so they are freed once gathered into availability order
    matrix[Ns, max(s)] delta = compute_delta(Ns, Y, s, trial, x, y, new_trl_flag, trial_start_idx);
    matrix[Ns, max(s)] phi   =   compute_phi(Ns, Y, s, trial, x, y, new_trl_flag, trial_start_idx);
    matrix[Ns, max(s)] psi   =   compute_psi(Ns, Y, s, trial, x, y, new_trl_flag, trial_start_idx);
    // Fourier basis for every selection, once
    array[Ns] matrix[H2, max(s)] basis =
      compute_fourier_basis_by_cond(Ns, s, trial, X, phi, harmonics, dir_offset);

    int idx = 1;
    for (t in 1:T) {
      int n_avail = s[t];
      for (ii in 1:f[t]) {
        array[n_avail] int avail = avail_order[idx, 1:n_avail];
        delta_g[idx, 1:n_avail] = delta[idx, avail];
        basis_g[idx][, 1:n_avail] = basis[idx][, avail];
        psi_g[idx, 1:n_avail] = psi[idx, avail];
        n_avail -= 1;
        idx += 1;
      }
    }
  }

  // trial ids, for reduce_sum to partition
  array[T] int trial_ids;
  for (t in 1:T) trial_ids[t] = t;

}

parameters {
  array[K] real<lower = 0> rho_delta; // distance
  array[K] real rho_psi; // relative direction
  matrix[K, H2] fourier_coef; // absolute direction: {cos, sin} per harmonic

}

model {
  rho_delta ~ normal(rho_delta_mu, rho_delta_sd);
  rho_psi ~ normal(rho_psi_mu, rho_psi_sd);
  for (k in 1:K) {
    to_vector(fourier_coef[k,]) ~ normal(to_vector(fourier_coef_mu[k,]), fourier_coef_sd);
  }

  // likelihood, summed over trials in parallel (reduce_sum)

  target += reduce_sum(partial_sum_lpmf, trial_ids, grainsize,
    f, s, X, sel_start_idx, p_at_step,
    rho_delta, fourier_coef, rho_psi, delta_g, basis_g, psi_g);

}

generated quantities {
  // prior draws

  array[K] real pr_rho_delta, pr_rho_psi;

  matrix[K, H2] pr_fourier_coef;

  for (k in 1:K) {
    pr_rho_delta[k] = normal_lower0_rng(rho_delta_mu[k], rho_delta_sd[k]);
    pr_rho_psi[k]   = normal_rng(rho_psi_mu[k], rho_psi_sd[k]);
    for (d in 1:H2) {
      pr_fourier_coef[k, d] = normal_rng(fourier_coef_mu[k, d], fourier_coef_sd[d]);
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

{
  int idx = 1;

  for (t in 1:T) {
    int k = X[t]; // condition
    int n_avail = s[t]; // items left

    for (ii in 1:f[t]) {
      row_vector[n_avail] weights = compute_weights(n_avail,
        rho_delta[k], fourier_coef[k,], rho_psi[k],
        delta_g[idx, 1:n_avail], basis_g[idx][, 1:n_avail], psi_g[idx, 1:n_avail]);

      int p = p_at_step[idx];

      if (save_selections) loglik[idx] = weights[p];
      log_lik[t] += weights[p];
      // P: an item sampled from the weights; W: the most likely item (indices
      // within the trial)
      if (save_selections) P[idx] = avail_order[idx, categorical_rng(exp(to_vector(weights)))];
      if (save_selections) W[idx] = avail_order[idx, which_max(weights)];

      n_avail -= 1;
      idx += 1;
      }
    }
  }
}
