// Fourier directional model: absolute direction as a truncated Fourier
// series in angle, added to the log-weights. Included with
// FoMo_functions.stan.
//
// The coefficients (fourier_coef's columns, and the basis rows) are, for
// each harmonic h in harmonics: cos(h a), then sin(h a) if has_sin for
// that harmonic (prep_stan_list(cos_only) drops a harmonic's sine term,
// e.g. harmonic 4's, which would rotate the 4-way pattern off the screen's
// axes). n_coef = n_harmonics + sum(has_sin).
//
// Contents:
//   compute_fourier_basis_by_cond(): the basis for every selection (fit models)
//   compute_fourier_basis_row(): the basis for one step (simulate models)
//   compute_fourier_absdir_weights(): coefficients x basis
//   get_person_fourier_coef(): one person's coefficients (multi-level)

int fourier_n_coef(array[] int has_sin) {
  // number of coefficients: a cos per harmonic, and the sines kept
  return size(has_sin) + sum(has_sin);
}

array[] matrix compute_fourier_basis_by_cond(int Ns, array[] int s,
  array[] int trial, array[] int X, matrix phi, array[] int harmonics,
  array[] int has_sin, vector dir_offset) {
  // per selection, a cos(h a) row (and a sin(h a) row, if has_sin) per
  // harmonic h in harmonics, with a = angle to the item minus the
  // condition's dir_offset. Harmonic 3 can be skipped ({1, 2, 4}): only
  // multiples of 4 give a 4-peaked pattern. Only the first s[trial[i]]
  // columns are filled

  int H = size(harmonics);
  array[Ns] matrix[fourier_n_coef(has_sin), max(s)] basis;

  for (i in 1:Ns) {
    int nt = s[trial[i]];
    real cond_offset = dir_offset[X[trial[i]]];
    row_vector[nt] a = phi[i, 1:nt] - cond_offset;
    int r = 1;
    for (hh in 1:H) {
      int h = harmonics[hh];
      basis[i][r, 1:nt] = cos(h * a);
      r += 1;
      if (has_sin[hh]) {
        basis[i][r, 1:nt] = sin(h * a);
        r += 1;
      }
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
  array[] int has_sin, real cond_offset, int nt) {
  // compute_fourier_basis_by_cond() for one step (phi: angle to each item)

  int H = size(harmonics);
  matrix[fourier_n_coef(has_sin), nt] basis;
  row_vector[nt] a = phi[1:nt] - cond_offset;
  int r = 1;

  for (hh in 1:H) {
    int h = harmonics[hh];
    basis[r,] = cos(h * a);
    r += 1;
    if (has_sin[hh]) {
      basis[r,] = sin(h * a);
      r += 1;
    }
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
