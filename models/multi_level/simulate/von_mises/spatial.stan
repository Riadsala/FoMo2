// FoMo multi-level spatial model, von Mises directional model: simulate.
//
// Each selection is a categorical choice among the items left, with
// log-weights from distance (rho_delta), absolute direction (theta:
// n_dircomps von Mises lobes, kappa fixed) and relative direction (rho_psi).
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

  // priors (or values to simulate from)

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
  vector<lower = 0>[K] sd_delta_lambda;
  vector<lower = 0>[K] sd_psi_lambda;
  vector<lower = 0>[n_dircomps] sd_log_theta_lambda;

}

transformed data {

}

parameters {
  // fixed effects
  array[K] real<lower = 0> rho_delta; // distance
  array[K] real rho_psi; // relative direction
  matrix[K, n_dircomps] log_theta; // absolute direction

  // random effects

  // between-person SDs
  vector<lower = 0>[K]    sigma_u_delta;
  vector<lower = 0>[K]    sigma_u_psi;
  matrix<lower = 0>[K, n_dircomps] sigma_u_log_theta;

  // per-person values
  array[K] vector<offset = person_offset(to_vector(rho_delta), L, 0),
                  multiplier = person_multiplier(sigma_u_delta, L, 0)>[L] u_delta;
  array[K] vector<offset = person_offset(to_vector(rho_psi), L, 0),
                  multiplier = person_multiplier(sigma_u_psi, L, 0)>[L] u_psi;
  array[K, n_dircomps] vector<offset = person_offset(log_theta, L, 0),
                              multiplier = person_multiplier(sigma_u_log_theta, L, 0)>[L] u_log_theta;
}

transformed parameters {
  // theta, for the output
  matrix[K, n_dircomps] theta = exp(log_theta);

  // per-person values

  array[K, n_dircomps] vector[L] u_theta;

  for (k in 1:K) {
    for (d in 1:n_dircomps) u_theta[k, d] = exp(u_log_theta[k, d]);
  }

}

model {
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

      for (ii in 1:nt) {
        row_vector[nt] weights = rep_row_vector(0, nt);

        if (ii > 1) {
          weights += -u_delta[k, z[t]] * delta[Q[pos+ii-2],];
          weights += compute_absdir_weights_fixed_kappa(person_theta, dircomps[Q[pos+ii-2]]);

        if (ii > 2) {
            // from the third: relative direction
            psi = compute_psi_from_phi(phi, Q[pos+ii-2], Q[pos+ii-3], nt);
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
