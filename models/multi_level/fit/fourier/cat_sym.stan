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
// b_v: probabilities; u_logit_v, u_v, u_first_v per person), so extract_post() and the
// plots read it unchanged. Selected by prep_stan_list(symmetric_v = TRUE).
// No simulate model.

functions {
  #include ../../../include/FoMo_functions.stan
  #include ../../../include/fourier_functions.stan

  vector stay_to_logodds(real stay, int c1, int n) {
    // row c1 of a symmetric transition matrix from its stay log-odds, as
    // log-odds against the last class (cat.stan's scale): every other
    // class has log-weight 0, class c1 has stay
    vector[n] w = rep_vector(0, n);
    w[c1] = stay;
    return w[1:(n - 1)] - w[n];
  }

  row_vector compute_weights(
    int n_avail,
    real rho_delta, row_vector coef, real rho_psi, row_vector log_b_v,
    row_vector d, matrix basis, row_vector psi, array[] int vt) {
      row_vector[n_avail] weights;

      weights = to_row_vector(-rho_delta * d);
      weights += compute_fourier_absdir_weights(coef, basis);
      weights +=  -rho_psi * psi;
      weights += log_b_v[vt];

      // every item passed in is available: normalise
      weights -= log_sum_exp(weights);

      return(weights);
  }

  // log-likelihood of a chunk of trials (reduce_sum); sel_start_idx[t]: trial
  // t's first row in the per-selection arrays
  real partial_sum_lpmf(
      array[] int trial_slice, int start, int end,
      array[] int f, array[] int s, array[] int X, array[] int z,
      array[] int sel_start_idx, array[] int p_at_step, int n_classes,
      array[] int Y, array[,] int vt_full,
      array[] vector u_delta, array[,] vector u_fourier_coef, array[] vector u_psi,
      array[,,] vector u_v, array[,] vector u_first_v,
      matrix delta_g, array[] matrix basis_g, matrix psi_g, array[,] int vt_g) {
    real lp = 0;

    for (t in trial_slice) {
      int k = X[t];
      int n_avail = s[t];
      int idx = sel_start_idx[t];

      row_vector[size(u_fourier_coef[k])] person_coef = get_person_fourier_coef(u_fourier_coef, k, z[t]);
      // log transition probabilities, once per trial
      matrix[n_classes, n_classes] log_person_v = log(get_person_v(u_v, k, z[t], n_classes));

      for (ii in 1:f[t]) {
        // previous selection's class: by item id, so from the ungathered
        // vt_full
        row_vector[n_avail] weights;
        if (ii == 1) {
          weights = first_item_weight(n_avail, u_first_v[k, z[t]], vt_g[idx, 1:n_avail]);
        } else {
          weights = compute_weights(n_avail,
            u_delta[k, z[t]], person_coef, u_psi[k, z[t]],
            log_person_v[vt_full[t, Y[idx-1]]],
            delta_g[idx, 1:n_avail], basis_g[idx][, 1:n_avail], psi_g[idx, 1:n_avail], vt_g[idx, 1:n_avail]);
        }

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

  // each trial's item classes by item id (vt_full) and, per selection, delta,
  // direction, psi and classes in availability order (unfilled beyond
  // n_avail)
  array[T, max(s)] int vt_full;
  matrix[Ns, max(s)] delta_g;
  array[Ns] matrix[H2, max(s)] basis_g;
  matrix[Ns, max(s)] psi_g;
  array[Ns, max(s)] int vt_g;

  {
    // distance, direction and turning angle by item id: local to this
    // block, so they are freed once gathered into availability order
    matrix[Ns, max(s)] delta = compute_delta(Ns, Y, s, trial, x, y, new_trl_flag, trial_start_idx);
    matrix[Ns, max(s)] phi   = compute_phi(Ns, Y, s, trial, x, y, new_trl_flag, trial_start_idx);
    matrix[Ns, max(s)] psi   = compute_psi(Ns, Y, s, trial, x, y, new_trl_flag, trial_start_idx);
    // Fourier basis for every selection, once
    array[Ns] matrix[H2, max(s)] basis =
      compute_fourier_basis_by_cond(Ns, s, trial, X, phi, harmonics, dir_offset);

    int pos = 1;
    int idx = 1;
    for (t in 1:T) {
      int nt = s[t];
      int n_avail = nt;

      array[nt] int vt = v[pos:(pos + nt - 1)];
      vt_full[t, 1:nt] = vt;
      pos = pos + nt;

      for (ii in 1:f[t]) {
        array[n_avail] int avail = avail_order[idx, 1:n_avail];
        delta_g[idx, 1:n_avail] = delta[idx, avail];
        basis_g[idx][, 1:n_avail] = basis[idx][, avail];
        psi_g[idx, 1:n_avail] = psi[idx, avail];
        vt_g[idx, 1:n_avail] = vt[avail];
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
    f, s, X, z, sel_start_idx, p_at_step, n_classes, Y, vt_full,
    u_delta, u_fourier_coef, u_psi, u_v, u_first_v, delta_g, basis_g, psi_g, vt_g);

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

  {
    int idx = 1;

    for (t in 1:T) {
      int k = X[t]; // condition
      int n_avail = s[t]; // items left

      // this person's values
      row_vector[H2] person_coef = get_person_fourier_coef(u_fourier_coef, k, z[t]);
      // log transition probabilities, once per trial
      matrix[n_classes, n_classes] log_person_v = log(get_person_v(u_v, k, z[t], n_classes));

      for (ii in 1:f[t]) {
        row_vector[n_avail] weights;
        if (ii == 1) {
          weights = first_item_weight(n_avail, u_first_v[k, z[t]], vt_g[idx, 1:n_avail]);
        } else {
          weights = compute_weights(n_avail,
            u_delta[k, z[t]], person_coef, u_psi[k, z[t]],
            log_person_v[vt_full[t, Y[idx-1]]],
            delta_g[idx, 1:n_avail], basis_g[idx][, 1:n_avail], psi_g[idx, 1:n_avail], vt_g[idx, 1:n_avail]);
        }

        int p = p_at_step[idx];

        if (save_selections) loglik[idx] = weights[p];
        log_lik[t] += weights[p];
        // P: an item sampled from the weights; W: the most likely item
        // (indices within the trial)
        if (save_selections) P[idx] = avail_order[idx, categorical_rng(exp(to_vector(weights)))];
        if (save_selections) W[idx] = avail_order[idx, which_max(weights)];

        n_avail -= 1;
        idx += 1;

      }
    }
  }
}
