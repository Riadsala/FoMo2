// FoMo shared functions, for every model whatever its directional model.
// Each model includes this file plus von_mises_functions.stan or
// fourier_functions.stan.
//
// Conventions:
//   - Y holds trial-local item indices (1..s[t]), not indices into the
//     N-length x, y, v; slice a trial's items first (from
//     trial_start_idx[t]), then index with Y.
//   - Per-selection matrices (compute_delta() etc.) have max(s) columns but
//     only the first s[t] are filled: always slice to 1:nt.
//   - new_trl_flag (`nts` here) is a selection's position within its
//     trial (1, 2, ...), not a 0/1 flag; it stops Y[i-1] / Y[i-2] reading
//     across a trial boundary.
//
// Contents:
//   which_max(), normalise_weights()
//   normal_lower0_rng(): prior draws of rho_delta (truncated at 0)
//   compute_delta(), compute_phi(), compute_psi(), transform_psi(),
//     compute_psi_from_phi(): distance and direction
//   compute_availability_order(), compute_selected_position(): which items
//     are left at each selection (fit models' transformed data)
//   get_logratio_mean(), compute_b_v_mean(), logodds_to_simplex(),
//     stationary_class_weights(): item classes
//   first_item_weight(): a trial's first selection
//   get_person_v(), create_trial_start_idx()
//   person_offset(), person_multiplier(), class_probs(): per-person
//     effects, centred or non-centred (multi-level)

real normal_lower0_rng(real mu, real sigma) {
  // a draw from normal(mu, sigma) truncated below at 0, as rho_delta's
  // prior is (its lower bound): the inverse CDF of a uniform draw above
  // the probability below 0
  real p0 = Phi(-mu / sigma);
  return mu + sigma * inv_Phi(uniform_rng(p0, 1));
}

int which_max(vector x) {
  // index of the largest element; ties go to the first
  int W = 1;
  for (j in 2:dims(x)[1]) {
    if (x[j] > x[W]) {
      W = j;
    }
  }

  return(W);

}

int which_max(row_vector x) {
  int W = 1;
  for (j in 2:dims(x)[2]) {
    if (x[j] > x[W]) {
      W = j;
    }
  }

  return(W);

}

row_vector normalise_weights(row_vector w, row_vector a) {
  // log-weights w -> log-probabilities over the available items
  // (a = 1 available, 0 selected; at least one must be available)

  row_vector[num_elements(w)] ws = w;

  ws += log(a); // -inf for selected items
  ws -= log_sum_exp(ws);

  return(ws);

}

matrix compute_delta(int Ns, array [] int Y, array [] int s, array [] int trial,
  vector x, vector y, array[] int nts, array[] int trial_start_idx) {
  // distance from the previous selection to every item in the trial, per
  // selection (row of zeros for a trial's first selection)

  matrix[Ns, max(s)] delta;

  for (i in 1:Ns) {

    int t = trial[i];
    int nt = s[t];

    if (nts[i] > 1) {

      int pos = trial_start_idx[t];
      vector[nt] xt = segment(x, pos, nt);
      vector[nt] yt = segment(y, pos, nt);
      delta[i, 1:nt] = to_row_vector(sqrt((xt - xt[Y[i-1]])^2 + (yt - yt[Y[i-1]])^2));

    } else {

      delta[i, 1:nt] = rep_row_vector(0, nt);

    }

  }

  return delta;

}

matrix compute_phi(int Ns, array [] int Y, array [] int s, array [] int trial,
  vector x, vector y, array [] int nts, array[] int trial_start_idx) {
  // absolute direction (radians, atan2) from the previous selection to every
  // item, per selection (zeros for a trial's first selection)

  matrix[Ns, max(s)] phi;

  for (i in 1:Ns) {

    int t = trial[i];
    int nt = s[t];

    if (nts[i] > 1) {

      int pos = trial_start_idx[t];
      vector[nt] xt = segment(x, pos, nt);
      vector[nt] yt = segment(y, pos, nt);
      phi[i, 1:nt] = to_row_vector(atan2( (yt - yt[Y[i-1]] ), (xt - xt[Y[i-1]]) ));

    } else {

      phi[i, 1:nt] = rep_row_vector(0, nt);

    }
  }

  return phi;

}

matrix compute_psi(int Ns, array [] int Y, array [] int s, array [] int trial,
  vector x, vector y, array [] int nts, array[] int trial_start_idx) {
  // relative direction, per selection: the angle to each item minus the
  // previous move's heading, folded to [0, 1] (0 = straight on, 1 = U-turn).
  // Needs two previous selections, so zeros for a trial's first two

  matrix[Ns, max(s)] psi;

  for (i in 1:Ns) {
    int t = trial[i];
    int nt = s[t];

    if (nts[i] > 2) {
      int pos = trial_start_idx[t];

      vector[nt] xt = segment(x, pos, nt);
      vector[nt] yt = segment(y, pos, nt);
      psi[i, 1:nt] = to_row_vector(atan2( yt - yt[Y[i-1]] , xt - xt[Y[i-1]] ));

      // heading of the previous move
      real prev_angle = atan2( (yt[Y[i-1]] - yt[Y[i-2]]), (xt[Y[i-1]] - xt[Y[i-2]]) ) ;
      psi[i, 1:nt] -= prev_angle;
      psi[i, 1:nt] = transform_psi(psi[i, 1:nt]);

      } else {

        psi[i, 1:nt] = rep_row_vector(0, nt);
    }
  }

  return psi;

}

row_vector transform_psi(row_vector psi) {
  // angle difference (radians) -> shortest arc, rescaled to [0, 1]

  return (pi() - abs(fmod(abs(psi), 2*pi()) - pi())) / pi();

}

row_vector compute_psi_from_phi(matrix phi, int prev_item, int prev_prev_item, int n) {
  // compute_psi() for one step of a simulation, from the trial's item-by-item
  // angle matrix (phi[a, b]: angle from item a to item b)

  row_vector[n] psi = phi[prev_item,] - phi[prev_prev_item, prev_item];

  return transform_psi(psi);

}

// compute_availability_order() and compute_selected_position() replay each
// trial's selections, removing the selected item by swapping it with the
// last available one. Data only, so the fit models run them once in
// transformed data and gather delta, phi, psi, ... into availability order
// there. Two functions rather than one returning both, as Stan's tuples
// can't hold unsized 2-D arrays.

array[,] int compute_availability_order(
    int Ns, int T, array[] int s, array[] int f, array[] int Y) {
  // avail_order[idx, 1:n_avail]: the items still available at selection idx

  int max_s = max(s);
  array[Ns, max_s] int avail_order;

  int idx = 1;
  for (t in 1:T) {

    int nt = s[t];
    array[nt] int avail_idx;
    array[nt] int pos_of;
    int n_avail = nt;

    for (i in 1:nt) {
      avail_idx[i] = i;
      pos_of[i] = i;
    }

    for (ii in 1:f[t]) {

      avail_order[idx, 1:n_avail] = avail_idx[1:n_avail];

      // remove the selected item: swap it with the last available one
      int p = pos_of[Y[idx]];
      int last_item = avail_idx[n_avail];
      avail_idx[p] = last_item;
      pos_of[last_item] = p;
      n_avail -= 1;

      idx += 1;

    }
  }

  return avail_order;
}

array[] int compute_selected_position(
    int Ns, int T, array[] int s, array[] int f, array[] int Y) {
  // p_at_step[idx]: the selected item's position in avail_order[idx, ]

  array[Ns] int p_at_step;

  int idx = 1;
  for (t in 1:T) {

    int nt = s[t];
    array[nt] int avail_idx;
    array[nt] int pos_of;
    int n_avail = nt;

    for (i in 1:nt) {
      avail_idx[i] = i;
      pos_of[i] = i;
    }

    for (ii in 1:f[t]) {

      int p = pos_of[Y[idx]];
      p_at_step[idx] = p;

      int last_item = avail_idx[n_avail];
      avail_idx[p] = last_item;
      pos_of[last_item] = p;
      n_avail -= 1;

      idx += 1;

    }
  }

  return p_at_step;
}

vector get_logratio_mean(matrix b_v_pr, int c1) {
  // column c1 of a transition-probability matrix -> log-odds against the
  // last class (the prior mean of that row's log-odds)

    int n_classes = rows(b_v_pr);

    vector[n_classes] p = col(b_v_pr, c1);
    p /= sum(p);

    return log(p[1:(n_classes - 1)]) - rep_vector(log(p[n_classes]), n_classes - 1);
  }


  array[,] vector compute_b_v_mean(
  array[] matrix b_v_pr,
  int n_classes) {
  // get_logratio_mean() for every [condition, previous class]

    int K = size(b_v_pr);

    array[K, n_classes] vector[n_classes - 1] out;

    for (k in 1:K)
      for (c1 in 1:n_classes)
        out[k, c1] =
          get_logratio_mean(b_v_pr[k], c1);

    return out;
}

vector logodds_to_simplex(vector lo) {
  // log-odds against the last class -> simplex (inverse of get_logratio_mean())

  return softmax(append_row(lo, 0));
}

array[,] vector logodds_to_simplex(array[,] vector lo) {
  // for every [condition, previous class] (single-level b_v_logodds -> b_v)

  int K = dims(lo)[1];
  int n = dims(lo)[2];
  array[K, n] vector[n] out;

  for (k in 1:K)
    for (c1 in 1:n)
      out[k, c1] = logodds_to_simplex(lo[k, c1]);

  return out;
}

vector stationary_class_weights(array[] vector rows) {
  // stationary distribution of a transition matrix (rows[c1]: the simplex
  // over the next class after c1): the first selection's class weights.
  // Chosen by loo over separate first-selection class weights and the
  // rows' geometric mean.
  // Solves (P' - I) pi = 0 with the last equation replaced by sum(pi) = 1

  int n = size(rows);
  matrix[n, n] A;
  vector[n] e = rep_vector(0, n);

  for (c1 in 1:n) A[, c1] = rows[c1];
  A -= diag_matrix(rep_vector(1, n));
  A[n] = rep_row_vector(1, n);
  e[n] = 1;

  return A \ e;
}

row_vector first_item_weight(int n_avail, vector class_w, array[] int vt) {
  // log-probabilities for a trial's first selection: no previous item, so
  // only the item's class weight (class_w; the place for any other
  // first-selection term, e.g. a start-position bias). class_w is a weight
  // per class, not P(starting with that class): that also depends on how
  // many items each class has

  row_vector[n_avail] weights = to_row_vector(log(class_w[vt]));

  return(weights - log_sum_exp(weights));
}

matrix get_person_v(
    array[,,] vector u_v,
    int k,
    int person,
    int n_classes) {
  // one person's transition matrix in condition k (row c: the simplex
  // after class c)

    matrix[n_classes, n_classes] out;

    for (d in 1:n_classes)
        out[d, ] = to_row_vector(u_v[k, person, d]);

    return out;
}

array[] int create_trial_start_idx(array[] int s, int T) {
  // position of each trial's first entry in an array concatenated over
  // trials (cumulative sum of s)

  array[T] int tsi;

  tsi[1] = 1;
  for (t in 2:T) {
    tsi[t] = tsi[t-1] + s[t-1];
  }

  return(tsi);
}


////////////////////////////////////////////////////////////////////////////////
// Per-person effects, centred or non-centred (multi-level models)
//
// Person-level parameters are declared with offset/multiplier, e.g.
//   vector<offset = person_offset(..., centre_delta),
//          multiplier = person_multiplier(..., centre_delta)>[L] u_delta;
//   u_delta[k] ~ normal(rho_delta[k], sigma_u_delta[k]);
// Offset = population value, multiplier = SD: non-centred. Offset 0,
// multiplier 1: centred. Same model either way, different geometry.
// person_offset() / person_multiplier() repeat those over L people;
// centred: 1 = centred, 0 = non-centred.
//
// Values are copied element by element, with if/else rather than ?:,
// because stanc O1 stores some parameter vectors in a layout that can't be
// assigned whole into a local array of vectors.
////////////////////////////////////////////////////////////////////////////////

// vector[K] -> array[K] vector[L]
array[] vector person_offset(vector pop, int L, int centred) {
  int K = num_elements(pop);
  array[K] vector[L] out;
  for (k in 1:K) {
    if (centred) out[k] = rep_vector(0, L);
    else out[k] = rep_vector(pop[k], L);
  }
  return out;
}
array[] vector person_multiplier(vector sd, int L, int centred) {
  int K = num_elements(sd);
  array[K] vector[L] out;
  for (k in 1:K) {
    if (centred) out[k] = rep_vector(1, L);
    else out[k] = rep_vector(sd[k], L);
  }
  return out;
}

// matrix[K, D] (direction components, item-feature dimensions)
// -> array[K, D] vector[L]
array[,] vector person_offset(matrix pop, int L, int centred) {
  int K = rows(pop);
  int D = cols(pop);
  array[K, D] vector[L] out;
  for (k in 1:K) for (d in 1:D) {
    if (centred) out[k, d] = rep_vector(0, L);
    else out[k, d] = rep_vector(pop[k, d], L);
  }
  return out;
}
array[,] vector person_multiplier(matrix sd, int L, int centred) {
  int K = rows(sd);
  int D = cols(sd);
  array[K, D] vector[L] out;
  for (k in 1:K) for (d in 1:D) {
    if (centred) out[k, d] = rep_vector(1, L);
    else out[k, d] = rep_vector(sd[k, d], L);
  }
  return out;
}

// array[K, n] vector[m] (transition log-odds) -> array[K, L, n] vector[m]
array[,,] vector person_offset(array[,] vector pop, int L, int centred) {
  int K = size(pop);
  int n = size(pop[1]);
  int m = num_elements(pop[1, 1]);
  array[K, L, n] vector[m] out;
  for (k in 1:K) for (l in 1:L) for (c in 1:n) {
    if (centred) out[k, l, c] = rep_vector(0, m);
    else for (j in 1:m) out[k, l, c, j] = pop[k, c, j];
  }
  return out;
}
array[,,] vector person_multiplier(array[,] vector sd, int L, int centred) {
  int K = size(sd);
  int n = size(sd[1]);
  int m = num_elements(sd[1, 1]);
  array[K, L, n] vector[m] out;
  for (k in 1:K) for (l in 1:L) for (c in 1:n) {
    if (centred) out[k, l, c] = rep_vector(1, m);
    else for (j in 1:m) out[k, l, c, j] = sd[k, c, j];
  }
  return out;
}

// per-person log-odds against the last class -> probabilities
array[,,] vector class_probs(array[,,] vector logits) {
  int K = size(logits);
  int L = size(logits[1]);
  int n = size(logits[1, 1]);
  array[K, L, n] vector[n] probs;
  for (k in 1:K) for (l in 1:L) for (c in 1:n) probs[k, l, c] = softmax(append_row(logits[k, l, c], 0));
  return probs;
}
