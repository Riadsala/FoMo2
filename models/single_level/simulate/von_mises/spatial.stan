// FoMo single-level spatial model, von Mises directional model: simulate.
//
// Each selection is a categorical choice among the items left, with
// log-weights from distance (rho_delta), absolute direction (theta:
// n_dircomps von Mises lobes, kappa fixed) and relative direction (rho_psi).
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

  array[T] int<lower = 1> s; // items per trial

  // condition
  array[T] int<lower = 1, upper = K> X;

  // item positions
  vector[N] x;
  vector[N] y;

    // absolute direction
  real<lower = 0> kappa; // von Mises concentration (FOMO_KAPPA)
  int<lower = 1> n_dircomps; // directions (4: E, N, W, S)
  vector[K] dir_offset; // absolute-direction offset per condition (radians)

  // priors (or values to simulate from)
  vector<lower = 0>[K] rho_delta_mu;
  vector<lower = 0>[K] rho_delta_sd;
  vector[K] rho_psi_mu;
  vector<lower = 0>[K] rho_psi_sd;
  matrix[K, n_dircomps] log_theta_mu;
  real log_theta_sd;

}

transformed data {

}

parameters {
  array[K] real<lower = 0> rho_delta; // distance
  array[K] real rho_psi; // relative direction
  matrix[K, n_dircomps] log_theta; // direction weights

}

transformed parameters {
  matrix[K, n_dircomps] theta = exp(log_theta);

}

model {
  rho_delta ~ normal(rho_delta_mu, rho_delta_sd);
  rho_psi ~ normal(rho_psi_mu, rho_psi_sd);
  to_vector(log_theta) ~ normal(to_vector(log_theta_mu), log_theta_sd);

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
      matrix[nt, nt] delta; // distance
      matrix[nt, nt] phi; // absolute direction
      row_vector[nt] psi; // relative direction

      array[nt] matrix[n_dircomps, nt] dircomps;
      for (i in 1:nt) {
        delta[i,] = to_row_vector(sqrt((xt - xt[i])^2 + (yt - yt[i])^2));
        phi[i,] =   to_row_vector(atan2((yt - yt[i]), (xt - xt[i])));
        dircomps[i] = compute_dir_components_row(phi[i,], kappa, n_dircomps, dir_offset[k], nt);
      }

      for (ii in 1:nt) {
        row_vector[nt] weights = rep_row_vector(0, nt);

        if (ii > 1) {
          // from the second selection: distance and absolute direction
          weights += -rho_delta[k] * delta[Q[pos+ii-2],];
          weights += compute_absdir_weights_fixed_kappa(theta[k,], dircomps[Q[pos+ii-2]]);

        if (ii > 2) {
            // from the third: relative direction
            psi = compute_psi_from_phi(phi, Q[pos+ii-2], Q[pos+ii-3], nt);
            weights += -rho_psi[k] * psi;
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
