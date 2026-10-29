// FoMo multi-level cat model, von Mises directional model: fit.
//
// Each selection is a categorical choice among the items left, with
// log-weights from distance (rho_delta), absolute direction (theta:
// n_dircomps von Mises lobes, kappa fixed) and relative direction (rho_psi),
// plus the item's class through the class-transition matrix b_v (row c: the
// next class after class c). A trial's first selection has class weights
// only (first_item_weight()): the stationary distribution of b_v.
//
// Multi-level: per-person effects on every term, each group centred or
// non-centred (centre_*: prep_stan_list(centred_params)).
//
// Data-only quantities (distances, directions, the items left) are computed
// once in transformed data; trials are summed with reduce_sum. Generated
// quantities: prior draws, per-trial log_lik (loo) and, with save_selections
// = 1, per selection P (sampled item), W (most likely item) and loglik.
// Simulate model: ../../simulate/von_mises/cat.stan.

functions {
  #include ../../../include/FoMo_functions.stan
  #include ../../../include/von_mises_functions.stan

   row_vector compute_weights(
    int n_avail,
    real rho_delta, row_vector theta, real rho_psi, row_vector log_b_v,
    row_vector d, matrix dircomps, row_vector psi, array[] int vt) {
      row_vector[n_avail] weights;

      weights = to_row_vector(-rho_delta * d);
      weights += compute_absdir_weights_fixed_kappa(theta, dircomps);
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
      array[] vector u_delta, array[,] vector u_theta, array[] vector u_psi,
      array[,,] vector u_v, array[,] vector u_first_v,
      matrix delta_g, array[] matrix dircomps_g, matrix psi_g, array[,] int vt_g) {
    real lp = 0;

    for (t in trial_slice) {
      int k = X[t];
      int n_avail = s[t];
      int idx = sel_start_idx[t];

      row_vector[size(u_theta[k])] person_theta = get_person_theta(u_theta, k, z[t]);
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
            u_delta[k, z[t]], person_theta, u_psi[k, z[t]],
            log_person_v[vt_full[t, Y[idx-1]]],
            delta_g[idx, 1:n_avail], dircomps_g[idx][, 1:n_avail], psi_g[idx, 1:n_avail], vt_g[idx, 1:n_avail]);
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
  // absolute direction
  real<lower = 0> kappa; // von Mises concentration (FOMO_KAPPA)
  int<lower = 1> n_dircomps; // directions (4: E, N, W, S)
  vector[K] dir_offset; // absolute-direction offset per condition (radians)

  // priors: absolute direction
  matrix[K, n_dircomps] log_theta_mu;
  real log_theta_sd;

  // priors: between-person SDs
  real<lower = 0> sd_v_lambda;
  vector<lower = 0>[K] sd_delta_lambda;
  vector<lower = 0>[K] sd_psi_lambda;
  vector<lower = 0>[n_dircomps] sd_log_theta_lambda;

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

  // start of each trial's items (cumulative s) and selections (cumulative f)
  array[T] int trial_start_idx = create_trial_start_idx(s, T);
  array[T] int sel_start_idx = create_trial_start_idx(f, T);

  // prior means of b_v_logodds: log-ratios within b_v_pr's columns
  array[K, n_classes] vector[n_classes - 1] b_v_mean;
  b_v_mean = compute_b_v_mean(b_v_pr, n_classes);

  // items left at each selection, and the selected item's position among them
  array[Ns, max(s)] int avail_order = compute_availability_order(Ns, T, s, f, Y);
  array[Ns] int p_at_step = compute_selected_position(Ns, T, s, f, Y);

  // each trial's item classes by item id (vt_full) and, per selection, delta,
  // direction, psi and classes in availability order (unfilled beyond
  // n_avail)
  array[T, max(s)] int vt_full;
  matrix[Ns, max(s)] delta_g;
  array[Ns] matrix[n_dircomps, max(s)] dircomps_g;
  matrix[Ns, max(s)] psi_g;
  array[Ns, max(s)] int vt_g;

  {
    // distance, direction and turning angle by item id: local to this
    // block, so they are freed once gathered into availability order
    matrix[Ns, max(s)] delta = compute_delta(Ns, Y, s, trial, x, y, new_trl_flag, trial_start_idx);
    matrix[Ns, max(s)] phi   = compute_phi(Ns, Y, s, trial, x, y, new_trl_flag, trial_start_idx);
    matrix[Ns, max(s)] psi   = compute_psi(Ns, Y, s, trial, x, y, new_trl_flag, trial_start_idx);
    // direction components for every selection, once
    array[Ns] matrix[n_dircomps, max(s)] dircomps =
      compute_dir_components_by_cond(Ns, s, trial, X, phi, kappa, n_dircomps, dir_offset);

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
        dircomps_g[idx][, 1:n_avail] = dircomps[idx][, avail];
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
  matrix[K, n_dircomps] log_theta; // absolute direction

  // class transitions: log-odds against the last class
  array[K, n_classes] vector[n_classes - 1] b_v_logodds;
  array[K, n_classes] vector<lower = 0>[n_classes - 1] sigma_u_v;

  // random effects

  // between-person SDs
  vector<lower = 0>[K] sigma_u_delta;
  vector<lower = 0>[K] sigma_u_psi;
  matrix<lower = 0>[K, n_dircomps] sigma_u_log_theta;

  // per-person values
  array[K, L, n_classes] vector<offset = person_offset(b_v_logodds, L, centre_v),
                                multiplier = person_multiplier(sigma_u_v, L, centre_v)>[n_classes - 1] u_logit_v;
  array[K] vector<offset = person_offset(to_vector(rho_delta), L, centre_delta),
                  multiplier = person_multiplier(sigma_u_delta, L, centre_delta)>[L] u_delta;
  array[K] vector<offset = person_offset(to_vector(rho_psi), L, centre_psi),
                  multiplier = person_multiplier(sigma_u_psi, L, centre_psi)>[L] u_psi;
  array[K, n_dircomps] vector<offset = person_offset(log_theta, L, centre_direction),
                              multiplier = person_multiplier(sigma_u_log_theta, L, centre_direction)>[L] u_log_theta;

}

transformed parameters {
  // class-transition probabilities, as in the single-level models: row c1
  // is a simplex over the next class, after an item of class c1
  array[K, n_classes] simplex[n_classes] b_v = logodds_to_simplex(b_v_logodds);
  matrix[K, n_dircomps] theta = exp(log_theta);

  array[K, n_dircomps] vector[L] u_theta;
  for (k in 1:K) {
    for (d in 1:n_dircomps) u_theta[k, d] = exp(u_log_theta[k, d]);
  }

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

  to_vector(log_theta) ~ normal(to_vector(log_theta_mu), log_theta_sd);

  // priors: between-person SDs
  sigma_u_delta ~ exponential(sd_delta_lambda);
  sigma_u_psi   ~ exponential(sd_psi_lambda);

  for (d in 1:n_dircomps) sigma_u_log_theta[, d] ~ exponential(sd_log_theta_lambda[d]);

  // per-person values around the population values (centred or not: see their
  // offset/multiplier)
  for (k in 1:K) u_delta[k] ~ normal(rho_delta[k], sigma_u_delta[k]);
  for (k in 1:K) u_psi[k] ~ normal(rho_psi[k], sigma_u_psi[k]);

  for (k in 1:K)
    for (d in 1:n_dircomps)
      u_log_theta[k, d] ~ normal(log_theta[k, d], sigma_u_log_theta[k, d]);

  for (k in 1:K)
    for (l in 1:L)
      for (c1 in 1:n_classes)
        u_logit_v[k, l, c1] ~ normal(b_v_logodds[k, c1], sigma_u_v[k, c1]);

  // likelihood, summed over trials in parallel (reduce_sum)

  target += reduce_sum(partial_sum_lpmf, trial_ids, grainsize,
    f, s, X, z, sel_start_idx, p_at_step, n_classes, Y, vt_full,
    u_delta, u_theta, u_psi, u_v, u_first_v, delta_g, dircomps_g, psi_g, vt_g);

}

generated quantities {
  // prior draws

  array[K] real pr_rho_delta, pr_rho_psi;
  matrix[K, n_dircomps] pr_theta;
  matrix[K, n_dircomps] pr_u_log_theta;
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
    for (d in 1:n_dircomps)
      pr_u_log_theta[k, d] = exponential_rng(sd_log_theta_lambda[d]);

    for (d in 1:n_dircomps)
      pr_theta[k, d] =
        lognormal_rng(log_theta_mu[k, d], log_theta_sd);

    for (c1 in 1:n_classes)
      pr_b_v[k, c1] =
        logodds_to_simplex(to_vector(normal_rng(b_v_mean[k, c1], b_v_sd)));
  }

  // for extract_post()
  vector[K] directional_offset = dir_offset;

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
      row_vector[n_dircomps] person_theta = get_person_theta(u_theta, k, z[t]);
      // log transition probabilities, once per trial
      matrix[n_classes, n_classes] log_person_v = log(get_person_v(u_v, k, z[t], n_classes));

      for (ii in 1:f[t]) {
        row_vector[n_avail] weights;
        if (ii == 1) {
          weights = first_item_weight(n_avail, u_first_v[k, z[t]], vt_g[idx, 1:n_avail]);
        } else {
          weights = compute_weights(n_avail,
            u_delta[k, z[t]], person_theta, u_psi[k, z[t]],
            log_person_v[vt_full[t, Y[idx-1]]],
            delta_g[idx, 1:n_avail], dircomps_g[idx][, 1:n_avail], psi_g[idx, 1:n_avail], vt_g[idx, 1:n_avail]);
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
