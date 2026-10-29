// FoMo multi-level cat model, von Mises directional model: simulate.
//
// Each selection is a categorical choice among the items left, with
// log-weights from distance (rho_delta), absolute direction (theta:
// n_dircomps von Mises lobes, kappa fixed) and relative direction (rho_psi),
// plus the item's class through the class-transition matrix b_v (row c: the
// next class after class c). A trial's first selection has class weights
// only (first_item_weight()): the stationary distribution of b_v.
//
// Multi-level: per-person effects on every term (always non-centred here).
//
// Simulates selection sequences Q from the priors (model block), or from a
// fit's draws (generate_quantities()). Parameters are declared in the same
// order as in the fit model, which generate_quantities() relies on.

functions {
  #include ../../../include/FoMo_functions.stan
  #include ../../../include/von_mises_functions.stan

}

data {

  int <lower = 1> N; // items, over all trials
  int <lower = 1> K; // conditions
  int <lower = 1> T; // trials
  int <lower = 1> L; // people

  array[T] int<lower = 1, upper = L> z; // person of each trial
  array[T] int<lower = 1, upper = K> X; // condition
  array[T] int<lower = 1> s; // items per trial

  // item positions
  vector[N] x;
  vector[N] y;
  // item classes
  int<lower = 1> n_classes;
  array[N] int<lower = 1, upper = n_classes> v;

  // priors (or values to simulate from)

  // priors: transition matrices
  array[K] matrix[n_classes, n_classes] b_v_pr;
  real<lower = 0> b_v_sd;
  // distance
  vector<lower = 0>[K] rho_delta_mu;
  vector<lower = 0>[K] rho_delta_sd;
  // relative direction
  vector[K] rho_psi_mu;
  vector<lower = 0>[K] rho_psi_sd;
  // absolute direction
  real<lower = 0> kappa; // von Mises concentration (FOMO_KAPPA)
  int<lower = 1> n_dircomps; // directions (4: E, N, W, S)
  vector[K] dir_offset; // absolute-direction offset per condition (radians)

  // absolute direction
  matrix[K, n_dircomps] log_theta_mu;
  real log_theta_sd;

  // priors: between-person SDs
  real<lower = 0> sd_v_lambda;
  vector<lower = 0>[K] sd_delta_lambda;
  vector<lower = 0>[K] sd_psi_lambda;
  vector<lower = 0>[n_dircomps] sd_log_theta_lambda;

}

transformed data {

  array[K, n_classes] vector[n_classes - 1] b_v_mean;
  b_v_mean = compute_b_v_mean(b_v_pr, n_classes);

}

parameters {
  // fixed effects
  array[K] real<lower = 0> rho_delta; // distance
  array[K] real rho_psi; // relative direction
  matrix[K, n_dircomps] log_theta; // absolute direction
  array[K, n_classes] vector[n_classes - 1] b_v_logodds; // class transitions: log-odds against the last class

  // same order as the fit model's parameters: generate_quantities() reads
  // them by position
  array[K, n_classes] vector<lower = 0>[n_classes - 1] sigma_u_v;

  // random effects

  // between-person SDs
  vector<lower = 0>[K]    sigma_u_delta;
  vector<lower = 0>[K]    sigma_u_psi;
  matrix<lower = 0>[K, n_dircomps] sigma_u_log_theta;

  // per-person values
  array[K, L, n_classes] vector<offset = person_offset(b_v_logodds, L, 0),
                                multiplier = person_multiplier(sigma_u_v, L, 0)>[n_classes - 1] u_logit_v;
  array[K] vector<offset = person_offset(to_vector(rho_delta), L, 0),
                  multiplier = person_multiplier(sigma_u_delta, L, 0)>[L] u_delta;
  array[K] vector<offset = person_offset(to_vector(rho_psi), L, 0),
                  multiplier = person_multiplier(sigma_u_psi, L, 0)>[L] u_psi;
  array[K, n_dircomps] vector<offset = person_offset(log_theta, L, 0),
                              multiplier = person_multiplier(sigma_u_log_theta, L, 0)>[L] u_log_theta;

}

transformed parameters {
  // class-transition probabilities, as in the single-level models: row c1
  // is a simplex over the next class, after an item of class c1
  array[K, n_classes] simplex[n_classes] b_v = logodds_to_simplex(b_v_logodds);
  // theta, for the output
  matrix[K, n_dircomps] theta = exp(log_theta);

  // per-person values

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

}

generated quantities {
  array[N] int Q = rep_array(0, N);

  // for extract_post()
  vector[K] directional_offset = dir_offset;

  {
    int pos = 1;

    for (t in 1:T) {
      int nt = s[t];
      int k = X[t];
      row_vector[nt] item_available = rep_row_vector(1, nt);
      // this trial's items
      vector[nt] xt = segment(x, pos, nt);
      vector[nt] yt = segment(y, pos, nt);
      array[nt] int vt = v[pos:(pos + nt - 1)];

      // item-by-item distance and direction
      matrix[nt, nt] delta; // distance
      matrix[nt, nt] phi; // absolute direction
      row_vector[nt] psi; // relative direction
      array[nt] matrix[n_dircomps, nt] dircomps;

      for (i in 1:nt) {
        delta[i,] = to_row_vector(sqrt((xt - xt[i])^2 + (yt - yt[i])^2));
        phi[i,] =   to_row_vector(atan2((yt - yt[i]), (xt - xt[i])));
        dircomps[i] = compute_dir_components_row(phi[i,], kappa, n_dircomps, dir_offset[k], nt);
      }

      // this person's values
      row_vector[n_dircomps] person_theta = get_person_theta(u_theta, k, z[t]);
      matrix[n_classes, n_classes] person_v = get_person_v(u_v, k, z[t], n_classes);

      for (ii in 1:nt) {
        row_vector[nt] weights = rep_row_vector(0, nt);

        // first selection: class weights only (see first_item_weight())
        if (ii == 1)
          weights = first_item_weight(nt, u_first_v[k, z[t]], vt);

        if (ii > 1) {
          int qm1 = Q[pos+ii-2]; // previous selection

          // transition weights from previously selected item's class
          weights += to_row_vector(log(person_v[vt[qm1],vt]));
          weights += -u_delta[k, z[t]] * delta[qm1,];
          weights += compute_absdir_weights_fixed_kappa(person_theta, dircomps[qm1]);

          if (ii > 2) {
            // from the third: relative direction
            psi = compute_psi_from_phi(phi, qm1, Q[pos+ii-3], nt);
            weights += -u_psi[k, z[t]] * psi;
          }
        }

      // log-probabilities over the items left
      weights = normalise_weights(weights, item_available);

      int q = categorical_rng(exp(to_vector(weights)));
      Q[pos+ii-1] = q;
      item_available[q] = 0;

    }

    pos = pos + s[t];

    }
  }
}
