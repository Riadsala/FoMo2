// von Mises directional model: absolute direction as a mixture of
// n_dircomps fixed-kappa von Mises lobes at evenly spaced directions,
// weighted by theta. Included with FoMo_functions.stan.
//
// Contents:
//   compute_absdir_weights_fixed_kappa(): log(1 + theta x components)
//   compute_dir_components_by_cond(): the components for every selection (fit models)
//   compute_dir_components_row(): the components for one step (simulate models)
//   get_person_theta(): one person's theta

row_vector compute_absdir_weights_fixed_kappa(row_vector theta, matrix dircomps) {
  // absolute-direction log-weights, log(1 + theta x dircomps): multiplies
  // each item's weight by the mixture. theta >= 0, so always finite. One
  // matrix product rather than a loop over directions (faster autodiff)
  return log1p(theta * dircomps);
}


array[] matrix compute_dir_components_by_cond(int Ns, array[] int s,
  array[] int trial, array[] int X, matrix phi, real kappa,
  int n_dircomps, vector dir_offset) {

  // per selection, the n_dircomps von Mises densities at the angle to each
  // item, centred on 2 pi (dir - 1) / n_dircomps + dir_offset[condition].
  // Only the first s[trial[i]] columns are filled

  real kappa_norm_factor = 1 / (2 * pi() * modified_bessel_first_kind(0, kappa));

  array[Ns] matrix[n_dircomps, max(s)] dircomps;

  for (i in 1:Ns) {
    int nt = s[trial[i]];
    real cond_offset = dir_offset[X[trial[i]]];
    for (dir in 1:n_dircomps) {
      real preferred_dir = (2 * pi() * (dir - 1) / n_dircomps) + cond_offset;
      dircomps[i][dir, 1:nt] = exp(kappa * cos(phi[i, 1:nt] - preferred_dir)) * kappa_norm_factor;
    }
  }

  return dircomps;
}

row_vector get_person_theta(
  array[,] vector u_theta,
  int k,
  int person) {
  // one person's theta in condition k

  int D = size(u_theta[k]);
  row_vector[D] theta;
  for (d in 1:D) theta[d] = u_theta[k, d][person];
  return theta;
}

matrix compute_dir_components_row(row_vector phi, real kappa, int n_dircomps,
  real cond_offset, int nt) {
  // compute_dir_components_by_cond() for one step (phi: angle to each item)

  real kappa_norm_factor = 1 / (2 * pi() * modified_bessel_first_kind(0, kappa));
  matrix[n_dircomps, nt] dircomps;

  for (dir in 1:n_dircomps) {
    real preferred_dir = (2 * pi() * (dir - 1) / n_dircomps) + cond_offset;
    dircomps[dir, ] = exp(kappa * cos(phi[1:nt] - preferred_dir)) * kappa_norm_factor;
  }

  return dircomps;
}
