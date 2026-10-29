// Fourier directional model: absolute direction as a truncated Fourier
// series in angle, added to the log-weights. Included with
// FoMo_functions.stan.
//
// Contents:
//   compute_fourier_basis_by_cond(): the basis for every selection (fit models)
//   compute_fourier_basis_row(): the basis for one step (simulate models)
//   compute_fourier_absdir_weights(): coefficients x basis
//   get_person_fourier_coef(): one person's coefficients (multi-level)

array[] matrix compute_fourier_basis_by_cond(int Ns, array[] int s,
  array[] int trial, array[] int X, matrix phi, array[] int harmonics,
  vector dir_offset) {
  // per selection, a {cos(h a), sin(h a)} pair of rows per harmonic h in
  // harmonics, with a = angle to the item minus the condition's dir_offset.
  // Harmonic 3 can be skipped ({1, 2, 4}): only multiples of 4 give a
  // 4-peaked pattern. Only the first s[trial[i]] columns are filled

  int H = size(harmonics);
  array[Ns] matrix[2 * H, max(s)] basis;

  for (i in 1:Ns) {
    int nt = s[trial[i]];
    real cond_offset = dir_offset[X[trial[i]]];
    row_vector[nt] a = phi[i, 1:nt] - cond_offset;
    for (hh in 1:H) {
      int h = harmonics[hh];
      basis[i][2*hh - 1, 1:nt] = cos(h * a);
      basis[i][2*hh,     1:nt] = sin(h * a);
    }
  }

  return basis;
}

row_vector compute_fourier_absdir_weights(row_vector coef, matrix basis) {
  // absolute-direction log-weights (any sign, so no log(1 + ...) as in
  // the von Mises model)
  return coef * basis;
}

matrix compute_fourier_basis_row(row_vector phi, array[] int harmonics,
  real cond_offset, int nt) {
  // compute_fourier_basis_by_cond() for one step (phi: angle to each item)

  int H = size(harmonics);
  matrix[2 * H, nt] basis;
  row_vector[nt] a = phi[1:nt] - cond_offset;

  for (hh in 1:H) {
    int h = harmonics[hh];
    basis[2*hh - 1,] = cos(h * a);
    basis[2*hh,]     = sin(h * a);
  }

  return basis;
}

row_vector get_person_fourier_coef(
  array[,] vector u_fourier_coef,
  int k,
  int person) {
  // one person's coefficients in condition k
  int D = size(u_fourier_coef[k]);
  row_vector[D] coef;
  for (d in 1:D) coef[d] = u_fourier_coef[k, d][person];
  return coef;
}
