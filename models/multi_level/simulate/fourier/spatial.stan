// FoMo multi-level spatial model, Fourier directional model: simulate.
//
// Each selection is a categorical choice among the items left, with
// log-weights from distance (rho_delta), absolute direction (fourier_coef:
// {cos, sin} pairs per harmonic) and relative direction (rho_psi).
//
// Multi-level: per-person effects on every term (always non-centred here).
//
// Simulates selection sequences Q from the priors (model block), or from a
// fit's draws (generate_quantities()). Parameters are declared in the same
// order as in the fit model, which generate_quantities() relies on.

functions {
  #include ../../../include/FoMo_functions.stan
  #include ../../../include/fourier_functions.stan

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

  // harmonics to fit, e.g. {1, 2, 4} (n_harmonics: Stan needs data array
  // sizes)
  int<lower = 1> n_harmonics;
  array[n_harmonics] int harmonics;
  vector[K] dir_offset; // absolute-direction offset per condition (radians)

  // absolute direction
  matrix[K, 2 * n_harmonics] fourier_coef_mu;
  vector[2 * n_harmonics] fourier_coef_sd;

  // priors: between-person SDs
  vector<lower = 0>[K] sd_delta_lambda;
  vector<lower = 0>[K] sd_psi_lambda;
  vector<lower = 0>[2 * n_harmonics] sd_fourier_coef_lambda;
}

transformed data {
  int H2 = 2 * n_harmonics; // width of fourier_coef (2 per harmonic)

}

parameters {
  // fixed effects
  array[K] real<lower = 0> rho_delta; // distance
  array[K] real rho_psi; // relative direction
  matrix[K, H2] fourier_coef; // absolute direction

  // random effects

  // between-person SDs
  vector<lower = 0>[K] sigma_u_delta;
  vector<lower = 0>[K] sigma_u_psi;
  matrix<lower = 0>[K, H2] sigma_u_fourier_coef;

  // per-person values
  array[K] vector<offset = person_offset(to_vector(rho_delta), L, 0),
                  multiplier = person_multiplier(sigma_u_delta, L, 0)>[L] u_delta;
  array[K] vector<offset = person_offset(to_vector(rho_psi), L, 0),
                  multiplier = person_multiplier(sigma_u_psi, L, 0)>[L] u_psi;
  array[K, H2] vector<offset = person_offset(fourier_coef, L, 0),
                      multiplier = person_multiplier(sigma_u_fourier_coef, L, 0)>[L] u_fourier_coef;
}

model {
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

}

generated quantities {
  array[N] int Q = rep_array(0, N);

  // for extract_post()
  vector[K] directional_offset = dir_offset;
  array[n_harmonics] int fourier_harmonics = harmonics; // labels fourier_coef's pairs (extract_fourier())

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
      array[nt] matrix[H2, nt] basis;

      for (i in 1:nt) {
        delta[i,] = to_row_vector(sqrt((xt - xt[i])^2 + (yt - yt[i])^2));
        phi[i,] =   to_row_vector(atan2((yt - yt[i]), (xt - xt[i])));
        basis[i] = compute_fourier_basis_row(phi[i,], harmonics, dir_offset[k], nt);
      }

      // this person's values
      row_vector[H2] person_coef = get_person_fourier_coef(u_fourier_coef, k, z[t]);

      for (ii in 1:nt) {
        row_vector[nt] weights = rep_row_vector(0, nt);

        if (ii > 1) {
          weights += -u_delta[k, z[t]] * delta[Q[pos+ii-2],];
          weights += compute_fourier_absdir_weights(person_coef, basis[Q[pos+ii-2]]);

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
