// FoMo single-level cat model, Fourier directional model: simulate.
//
// Each selection is a categorical choice among the items left, with
// log-weights from distance (rho_delta), absolute direction (fourier_coef:
// {cos, sin} pairs per harmonic) and relative direction (rho_psi), plus the
// item's class through the class-transition matrix b_v (row c: the next
// class after class c). A trial's first selection has class weights only
// (first_item_weight()): the stationary distribution of b_v.
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

  array[T] int<lower = 1> s; // items per trial

  // condition
  array[T] int<lower = 1, upper = K> X;

  // item positions
  vector[N] x;
  vector[N] y;

  // item classes
  int<lower = 1> n_classes;
  array[N] int<lower = 1, upper = n_classes> v;

  // harmonics to fit, e.g. {1, 2, 4} (n_harmonics: Stan needs data array
  // sizes)
  int<lower = 1> n_harmonics;
  array[n_harmonics] int harmonics;
  vector[K] dir_offset; // absolute-direction offset per condition (radians)

  // priors (or values to simulate from)
  array[K] matrix[n_classes, n_classes] b_v_pr; // b_v_logodds' prior means are its columns' log-ratios
  real<lower = 0> b_v_sd; // prior SD of b_v_logodds
  vector<lower = 0>[K] rho_delta_mu;
  vector<lower = 0>[K] rho_delta_sd;
  vector[K] rho_psi_mu;
  vector<lower = 0>[K] rho_psi_sd;
  matrix[K, 2 * n_harmonics] fourier_coef_mu;
  vector[2 * n_harmonics] fourier_coef_sd;

}

transformed data {
  // prior means of b_v_logodds: the log-ratios within b_v_pr's columns
  array[K, n_classes] vector[n_classes - 1] b_v_mean = compute_b_v_mean(b_v_pr, n_classes);

  int H2 = 2 * n_harmonics; // width of fourier_coef (2 per harmonic)

}

parameters {
  array[K] real<lower = 0> rho_delta; // distance
  array[K] real rho_psi; // relative direction
  matrix[K, H2] fourier_coef; // absolute direction: {cos, sin} per harmonic
  // class transitions: row c1's log-odds of each class against the last
  // class (the transformed parameter b_v holds the probabilities)
  array[K, n_classes] vector[n_classes - 1] b_v_logodds;

}

transformed parameters {
  // class-transition probabilities: row c1 is a simplex over the next
  // class after class c1
  array[K, n_classes] simplex[n_classes] b_v = logodds_to_simplex(b_v_logodds);
  // the first selection's class weights: the stationary distribution of
  // b_v (see stationary_class_weights())
  array[K] simplex[n_classes] first_v;
  for (k in 1:K) first_v[k] = stationary_class_weights(b_v[k]);

}

model {
  // b_v: a normal prior on each row's log-odds against the last class
  for (k in 1:K)
    for (c1 in 1:n_classes)
      b_v_logodds[k, c1] ~ normal(b_v_mean[k, c1], b_v_sd);

  rho_delta ~ normal(rho_delta_mu, rho_delta_sd);
  rho_psi   ~ normal(rho_psi_mu,   rho_psi_sd);

  for (k in 1:K) {
    to_vector(fourier_coef[k,]) ~ normal(to_vector(fourier_coef_mu[k,]), fourier_coef_sd);
  }

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
      array[nt] int vt = v[pos:(pos + nt - 1)];
      matrix[nt, nt] delta; // distance
      matrix[nt, nt] phi; // absolute direction
      row_vector[nt] psi; // relative direction

      array[nt] matrix[H2, nt] basis;
      for (i in 1:nt) {
        delta[i,] = to_row_vector(sqrt((xt - xt[i])^2 + (yt - yt[i])^2));
        phi[i,] = to_row_vector(atan2((yt - yt[i]), (xt - xt[i])));
        basis[i] = compute_fourier_basis_row(phi[i,], harmonics, dir_offset[k], nt);
      }

      for (ii in 1:nt) {
        row_vector[nt] weights = rep_row_vector(0, nt);

        // first selection: class weights only (first_item_weight())
        if (ii == 1) {weights = first_item_weight(nt, first_v[k], vt);}

        if (ii > 1) {
          weights += to_row_vector(log(b_v[k, vt[Q[pos+ii-2]]][vt]));
          // from the second selection: distance and absolute direction
          weights += -rho_delta[k] * delta[Q[pos+ii-2],];
          weights += compute_fourier_absdir_weights(fourier_coef[k,], basis[Q[pos+ii-2]]);

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
