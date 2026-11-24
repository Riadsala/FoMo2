// The "flat" likelihood layout of the multi-level fit models: the index
// bookkeeping they share. Included with FoMo_functions.stan.
//
// Each selection is a categorical choice among the items left, with
// log-weights that are, apart from the class transitions, sums of
// (parameter x data) terms: distance, turn, and the direction term
// (Fourier basis, or von Mises lobes inside log1p()). The data side of
// every term is computed once, in transformed data, as one column per
// item left at each selection: one matrix (F) with a row per feature. A
// trial's selections are consecutive, so each trial owns one contiguous
// block of columns, and the log-weights of all its selections come from
// one matrix product per trial (the person's coefficients times the
// trial's block), instead of one small product per selection. Within the
// block, each selection's items are a segment; the selection adds
//   w[chosen item] - log_sum_exp(w[its segment])
// to the log-likelihood.
//
//   trial t:    | sel a: s-a+1 items | sel a+1: s-a items | ... | sel f |
//   F columns:  ^ flat_start[t]                                        ^
//               <---------------------- flat_len[t] ------------------->
//   seg_rel:    where each selection's segment starts, counted within the
//               trial's block; chosen_rel: where its chosen item is
//
// a = first_sel: 1 if every selection is in F (spatial models), 2 if a
// trial's first selection isn't (cat models: it has class weights only).
// Items within a segment are in availability order
// (compute_availability_order()), so the chosen item is at position
// p_at_step within its segment.
//
// The class-transition term is a lookup: the transition matrix is
// flattened into a vector of n_classes^2 log-probabilities
// (flat_log_v()), and each column's index into it is
// (previous class - 1) * n_classes + item class (flat_class_index()).
//
// Same log density as computing each selection separately (identical to
// ~1e-10), 3-4 times faster (benchmarks, November 2026). The simulate
// models compute each selection separately: the readable reference.
//
// Contents:
//   count_flat(), flat_trial_start(), flat_trial_len(): the columns
//   flat_seg_rel(), flat_chosen_rel(): each selection's segment and choice
//   flat_cols(): each column's selection and item (to fill F)
//   flat_class_index(), flat_log_v(): the class-transition lookup
//   flat_trial_loglik(): a trial's log-likelihood from its log-weights
//   flat_selection_weights(): one selection's normalised log-weights (GQ)

int count_flat(array[] int s, array[] int f, int T, int first_sel) {
  // number of columns: over trials, the items left at each selection from
  // first_sel on
  int M = 0;
  for (t in 1:T)
    for (ii in first_sel:f[t]) M += s[t] - ii + 1;
  return M;
}

array[] int flat_trial_start(array[] int s, array[] int f, int T, int first_sel) {
  // each trial's first column
  array[T] int st;
  int col = 1;
  for (t in 1:T) {
    st[t] = col;
    for (ii in first_sel:f[t]) col += s[t] - ii + 1;
  }
  return st;
}

array[] int flat_trial_len(array[] int s, array[] int f, int T, int first_sel) {
  // each trial's number of columns (0 if it has no selection from first_sel on)
  array[T] int len;
  for (t in 1:T) {
    len[t] = 0;
    for (ii in first_sel:f[t]) len[t] += s[t] - ii + 1;
  }
  return len;
}

array[] int flat_seg_rel(array[] int s, array[] int f, int T, int Ns, int first_sel) {
  // each selection's first column, counted within its trial's block (0 for
  // selections before first_sel, which aren't in F)
  array[Ns] int seg = rep_array(0, Ns);
  int idx = 1;
  for (t in 1:T) {
    int o = 1;
    for (ii in 1:f[t]) {
      if (ii >= first_sel) {
        seg[idx] = o;
        o += s[t] - ii + 1;
      }
      idx += 1;
    }
  }
  return seg;
}

array[] int flat_chosen_rel(array[] int seg_rel, array[] int p_at_step) {
  // each selection's chosen item's column, counted within its trial's block
  // (0 where seg_rel is)
  int Ns = size(seg_rel);
  array[Ns] int ch;
  for (i in 1:Ns) ch[i] = seg_rel[i] > 0 ? seg_rel[i] + p_at_step[i] - 1 : 0;
  return ch;
}

array[,] int flat_cols(array[] int s, array[] int f, int T, int first_sel,
                       array[,] int avail_order, int M) {
  // for each column: its selection (row of the per-selection arrays) and
  // its item (index within the trial), in availability order
  array[M, 2] int cols;
  int idx = 1;
  int col = 1;
  for (t in 1:T) {
    for (ii in 1:f[t]) {
      if (ii >= first_sel) {
        for (j in 1:(s[t] - ii + 1)) {
          cols[col, 1] = idx;
          cols[col, 2] = avail_order[idx, j];
          col += 1;
        }
      }
      idx += 1;
    }
  }
  return cols;
}

array[] int flat_class_index(array[,] int cols, array[] int Y, array[] int v,
                             array[] int trial, array[] int trial_start_idx,
                             int n_classes) {
  // each column's index into flat_log_v(): (previous class - 1) * n_classes
  // + item class. The previous selection is the one before the column's,
  // so cols must not include first selections (first_sel = 2)
  int M = dims(cols)[1];
  array[M] int cls;
  for (c in 1:M) {
    int sel = cols[c, 1];
    int v0 = trial_start_idx[trial[sel]] - 1;   // v[v0 + i]: the trial's item i
    int prev_class = v[v0 + Y[sel - 1]];
    cls[c] = (prev_class - 1) * n_classes + v[v0 + cols[c, 2]];
  }
  return cls;
}

vector flat_log_v(array[] vector person_v) {
  // a transition matrix (one simplex per previous class) as a vector of
  // log-probabilities: element (previous - 1) * n_classes + next
  int n = size(person_v);
  vector[n * n] lv;
  for (c1 in 1:n) lv[((c1 - 1) * n + 1):(c1 * n)] = log(person_v[c1]);
  return lv;
}

real flat_trial_loglik(row_vector w, int idx, int s_t, int f_t, int first_sel,
                       array[] int seg_rel, array[] int chosen_rel) {
  // log-likelihood of a trial's selections from first_sel on, given w (the
  // log-weights of the trial's block, unnormalised). idx: the trial's
  // first selection's row; s_t, f_t: its items and selections
  real lp = 0;
  for (ii in first_sel:f_t) {
    int sidx = idx + ii - 1;
    int o = seg_rel[sidx];
    int n = s_t - ii + 1;   // items left
    lp += w[chosen_rel[sidx]] - log_sum_exp(w[o:(o + n - 1)]);
  }
  return lp;
}

row_vector flat_selection_weights(row_vector w, int seg_start, int n_avail) {
  // one selection's normalised log-weights (availability order), from its
  // segment of the trial's log-weights w (generated quantities)
  row_vector[n_avail] weights = w[seg_start:(seg_start + n_avail - 1)];
  return weights - log_sum_exp(weights);
}
