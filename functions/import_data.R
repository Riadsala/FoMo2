# These functions import each dataset and do some initial processing to get
# it into our stim and found format (see data/README.md).
#
# import_data() returns a list with two data frames:
#   stim:  every item in every trial
#   found: the items that were selected, in selection order; found$n is
#          each selection's position in its trial (1 = first item selected,
#          2 = second, ...)
#
# (0, 0) is the bottom-left corner of the display, which is 1 unit wide
# (its height depends on the aspect ratio).

# shorthand dataset names (first-author initial + year; see CLAUDE.md),
# accepted by import_data() in any case, e.g. import_data("K2014")
DATASET_SHORTHAND <- c(
  K2014 = "kristjansson2014plos",
  H2024 = "hughes2024rsos",
  C2022 = "clarke2022qjep",
  T2022 = "tagu2022cog")

import_data <- function(dataset, ...) {

  # ...: passed to the importer (e.g. scarcity for hughes2024rsos)
  importers <- list(clarke2022qjep       = import_clarke2022qjep,
                    tagu2022cog          = import_tagu2022cog,
                    tagu2022cog2items    = import_tagu2022cog2items,
                    kristjansson2014plos = import_kristjansson2014plos,
                    hughes2024rsos       = import_hughes2024rsos)

  if (toupper(dataset) %in% names(DATASET_SHORTHAND)) {
    dataset <- DATASET_SHORTHAND[[toupper(dataset)]]
  }
  if (!(dataset %in% names(importers))) {
    stop("unknown dataset '", dataset, "'; known datasets: ",
         paste(names(importers), collapse = ", "), " (or shorthand: ",
         paste(names(DATASET_SHORTHAND), collapse = ", "), ")")
  }

  data_path <- file.path(rprojroot::find_root(rprojroot::has_dir("models")), "data")
  d <- importers[[dataset]](data_path, ...)

  # every dataset should end up with the same column names:
  #   n - position of each selection within its trial (1 = first item selected,
  #       2 = second, ...); read by prep_stan_list() as found$n
  #   v - item class: a factor (cat models)
  # the importers above build these directly, but the ones that read
  # pre-processed csv file (hughes2024rsos) may still carry the
  # old names, so catch those here.
  for (nm in c("stim", "found")) {
    if ("item_class" %in% names(d[[nm]]) && !("v" %in% names(d[[nm]]))) {
      d[[nm]] <- d[[nm]] |> rename(v = item_class) |> mutate(v = as.factor(v))
    }
  }

  if (!("n" %in% names(d$found))) {
    d$found <- if ("found" %in% names(d$found)) {
      rename(d$found, n = found)
    } else {
      mutate(d$found, n = row_number(), .by = trial)
    }
  }

  verify_data(d$stim, d$found, dataset)

  list(dataset = dataset, stim = d$stim, found = d$found)
}

verify_data <- function(stim, found = NULL, label = "data", show_notes = TRUE) {

  # Checks stim and found against the format in data/README.md, which the
  # Stan models rely on: most breaches would give silently misaligned data
  # rather than an error. Stops with every problem found. With show_notes,
  # also prints a note for things that are allowed but worth knowing (a
  # person missing from a condition, trials with one selection, items at
  # the same position). found = NULL checks stim alone (e.g. to simulate).
  # Run by import_data() (with notes) and by prep_stan_list() (without, so
  # every data list is checked, simulated or subsetted data included).
  # Returns TRUE invisibly.

  problems <- character(0)
  add <- function(ok, msg) if (!isTRUE(ok)) problems <<- c(problems, msg)
  report <- function() {
    if (length(problems) > 0) {
      stop(label, " fails verify_data():\n  - ", paste(problems, collapse = "\n  - "), call. = FALSE)
    }
  }
  consecutive <- function(x) all(x == seq_along(x))

  # --- columns: stop here if any are missing, since the rest need them
  missing_cols <- c(sprintf("stim$%s", setdiff(c("person", "condition", "trial", "id", "x", "y"), names(stim))),
                    if (!is.null(found)) sprintf("found$%s", setdiff(c("trial", "id", "n"), names(found))))
  add(length(missing_cols) == 0, paste("missing columns:", paste(missing_cols, collapse = ", ")))
  report()

  # --- stim
  add(!anyNA(select(stim, person, condition, trial, id, x, y)), "NAs in stim")
  add(all(is.finite(stim$x) & is.finite(stim$y)), "non-finite stim$x or stim$y")
  add(is.factor(stim$condition), "stim$condition is not a factor")
  add(all(levels(stim$condition) %in% stim$condition), "stim$condition has unused levels")
  add(setequal(stim$person, seq_len(n_distinct(stim$person))), "stim$person is not 1..L with no gaps")
  add(setequal(stim$trial, seq_len(n_distinct(stim$trial))), "stim$trial is not 1..T with no gaps")
  add(!is.unsorted(stim$trial), "stim is not sorted by trial")
  add(all(stim$x >= 0 & stim$x <= 1, na.rm = TRUE) && all(stim$y >= 0, na.rm = TRUE), "stim$x outside [0, 1] or stim$y negative")

  stim_trials <- stim |>
    summarise(n_items = n(),
              n_person = n_distinct(person),
              n_condition = n_distinct(condition),
              ids_ok = consecutive(id),
              dup_xy = anyDuplicated(paste(x, y)) > 0,
              .by = trial)
  add(all(stim_trials$n_person == 1), "a trial has more than one person")
  add(all(stim_trials$n_condition == 1), "a trial has more than one condition")
  add(all(stim_trials$ids_ok), "stim$id is not 1..nt in row order within every trial")

  # item features: prep_stan_list() passes a factor's codes to Stan as
  # classes 1..n_classes, with n_classes the most classes in any condition
  if ("v" %in% names(stim)) {
    add(is.factor(stim$v), "stim$v is not a factor")
    add(!anyNA(stim$v), "NAs in stim$v")
    n_classes <- max(summarise(stim, n = n_distinct(v), .by = condition)$n)
    add(max(as.integer(stim$v), na.rm = TRUE) <= n_classes,
        paste0("stim$v has factor codes above n_classes (", n_classes,
               "): unused levels, or different classes in different conditions"))
  }

  # --- found
  found_trials <- NULL
  if (!is.null(found)) {
    add(!anyNA(select(found, trial, id, n)), "NAs in found$trial, found$id or found$n")
    add(!is.unsorted(found$trial), "found is not sorted by trial")

    found_trials <- found |>
      summarise(f = n(), n_ok = consecutive(n), repeats = anyDuplicated(id) > 0,
                max_id = max(id), min_id = min(id), .by = trial) |>
      left_join(stim_trials, by = "trial")
    empty <- setdiff(stim_trials$trial, found_trials$trial)
    add(length(empty) == 0,
        paste0(length(empty), " trial(s) with no selections (e.g. trial ", empty[1], ")"))
    add(all(found_trials$trial %in% stim_trials$trial), "found has trials that aren't in stim")
    add(all(found_trials$n_ok), "found$n is not 1..f in order within every trial")
    add(!any(found_trials$repeats), "an item is selected more than once in a trial")
    add(all(found_trials$min_id >= 1 & found_trials$max_id <= found_trials$n_items, na.rm = TRUE),
        "found$id outside 1..(items in its trial)")
  }

  # columns carried over from stim must agree with the item they point to
  shared <- intersect(c("person", "condition", "x", "y", "v"), names(found))
  if (length(shared) > 0 && length(problems) == 0) {
    rows <- match(paste(found$trial, found$id), paste(stim$trial, stim$id))
    for (col in shared) {
      same <- if (is.numeric(stim[[col]])) {
        all(abs(found[[col]] - stim[[col]][rows]) < 1e-9)
      } else {
        all(as.character(found[[col]]) == as.character(stim[[col]][rows]))
      }
      add(isTRUE(same), paste0("found$", col, " disagrees with stim for the selected item"))
    }
    # extract_post() takes its condition labels from found
    if ("condition" %in% shared) {
      add(identical(levels(found$condition), levels(stim$condition)),
          "found$condition and stim$condition have different levels")
    }
  }
  report()

  if (!show_notes) return(invisible(TRUE))

  # --- allowed, but worth knowing
  notes <- character(0)
  note <- function(count, msg) if (count > 0) notes <<- c(notes, paste(count, msg))
  n_missing_cells <- n_distinct(stim$person) * nlevels(stim$condition) - nrow(distinct(stim, person, condition))
  note(n_missing_cells, "person x condition combination(s) with no trials")
  note(sum(found_trials$f == 1), "trial(s) with only one selection")
  note(sum(stim_trials$dup_xy), "trial(s) with two items at the same position")
  if (length(notes) > 0) message("verify_data() notes for ", label, ": ", paste(notes, collapse = "; "))

  invisible(TRUE)
}

########################################################################
# shared steps

filter_one_person <- function(d, pp) {

  # one person's data, as person 1 with trials numbered 1, 2, ...
  map(list(found = d$found, stim = d$stim), function(x) {
    x |>
      filter(person == pp) |>
      mutate(person = 1, trial = dense_rank(trial))
  })
}

subset_trials <- function(d, trials) {

  # the trials of d (list(stim, found)) listed in `trials` (a table with
  # person, condition and trial columns), with people renumbered 1..L and
  # trials 1..T from one key, so stim and found stay matched (running
  # fix_person_and_trial() on each table separately only lines up if both
  # keep exactly the same trials). trial_p is left as it was
  key <- trials |>
    distinct(person, condition, trial) |>
    arrange(person, condition, trial) |>
    mutate(new_person = consecutive_id(person), new_trial = row_number())

  renumber <- function(x) {
    x |>
      inner_join(key, by = c("person", "condition", "trial")) |>
      mutate(person = new_person, trial = new_trial) |>
      select(-new_person, -new_trial)
  }
  list(stim = renumber(d$stim) |> arrange(trial, id),
       found = renumber(d$found) |> arrange(trial, n))
}

fix_person_and_trial <- function(d) {

  # sort by person, condition and trial; number the people 1..L with no
  # gaps, and the trials 1..T over people and conditions (the original
  # trial number is kept as trial_p)
  d |>
    arrange(person, condition, trial) |>
    mutate(person = consecutive_id(person),
           trial_p = trial,
           trial = consecutive_id(person, condition, trial))
}

read_foraging_csv <- function(file, ...) {

  # the csv format shared by the tagu2022cog, clarke2022qjep and
  # kristjansson2014plos data: one row per item, with `found` its position
  # in the selection order (0 if not selected). ...: types of any extra
  # columns
  read_csv(file, col_types = cols(condition = col_double(),
                                  trial = col_double(),
                                  observer = col_double(),
                                  targ_type = col_double(),
                                  id = col_double(),
                                  x = col_double(),
                                  y = col_double(),
                                  found = col_double(),
                                  ...))
}

drop_na_trials <- function(d) {
  # remove every trial with a missing item class
  filter(d, !any(is.na(targ_type)), .by = c(condition, trial, observer))
}

rescale_xy <- function(d) {
  # shift to (0, 0) and scale by the width: x in [0, 1], y in [0, a],
  # where a is the aspect ratio
  width <- max(d$x) - min(d$x)
  mutate(d, x = (x - min(x)) / width, y = (y - min(y)) / width)
}

flip_y <- function(d) {
  # screen coordinates (y down) -> y up, so that (0, 0) is the bottom left
  mutate(d, y = 1 - y, y = y - min(y))
}

split_stim_found <- function(d) {
  # stim: every item; found: the selected items, in selection order
  list(stim = d |>
         select(person, condition, trial, id, x, y, v, trial_p) |>
         arrange(person, condition, trial, id),
       found = d |>
         filter(n > 0) |>
         arrange(person, condition, trial, n))
}

########################################################################
# one importer per dataset

import_tagu2022cog2items <- function(data_path) {

  # item classes: targ_type 0 (high value) is class 2, the rest class 1
  read_foraging_csv(file.path(data_path, "tagu2022cog", "tagu_2020_prox_mouse.csv")) |>
    drop_na_trials() |>
    select(person = observer, condition, trial, id, n = found, v = targ_type, x, y) |>
    mutate(v = as.factor(if_else(v == 0, 2, 1)),
           condition = fct_recode(as.factor(condition), value = "1", control = "2")) |>
    fix_person_and_trial() |>
    rescale_xy() |>
    flip_y() |>
    split_stim_found()
}

import_tagu2022cog <- function(data_path) {

  # like import_tagu2022cog2items(), but keeps all 3 target colours per block
  # as distinct item classes instead of collapsing to a binary
  # high/low-value flag.
  #
  # raw colour identity (1-6) can't be used as the item class directly:
  #   - it isn't comparable across conditions - target/distractor identities
  #     are reversed between Value and NoValue blocks for a given
  #     participant, so colour "3" means something different in each
  #   - it isn't even Stan-valid on its own - each block only uses 3 of the
  #     6 possible colours, and they're never a contiguous 1..3 range (e.g.
  #     {1,3,5}), but prep_stan_list()'s cat branch declares v with
  #     upper = n_classes = 3
  # so targ_type (raw colour) is recoded here into a "value role": 1 = this
  # block's designated colour (high_value_colour - the actual high-value
  # category in Value blocks, the nominal patch category in NoValue
  # blocks), 2/3 = the other two colours, in a fixed within-block order
  # (sorted by raw colour number). high_value_colour is constant across all
  # 10 trials of a given (observer, condition) block (verified directly
  # against the raw data), so this recoding is well-defined and stable.

  read_foraging_csv(file.path(data_path, "tagu2022cog", "tagu_2020_prox_mouse_allcolours.csv"),
                    high_value_colour = col_double()) |>
    drop_na_trials() -> d

  # for each (observer, condition) block, the 2 non-high-value colours get
  # roles 2 and 3, in colour-number order
  d |>
    filter(targ_type != high_value_colour) |>
    distinct(observer, condition, targ_type) |>
    arrange(observer, condition, targ_type) |>
    mutate(v = row_number() + 1, .by = c(observer, condition)) -> role_lookup

  stopifnot("expected exactly 2 non-high-value colours per (observer, condition) block" =
              all(count(role_lookup, observer, condition)$n == 2))

  d |>
    left_join(role_lookup, by = c("observer", "condition", "targ_type")) |>
    mutate(v = if_else(targ_type == high_value_colour, 1, v)) |>
    select(person = observer, condition, trial, id, n = found, v, x, y) |>
    mutate(v = as.factor(v),
           condition = fct_recode(as.factor(condition), value = "1", control = "2")) |>
    fix_person_and_trial() |>
    rescale_xy() |>
    flip_y() |>
    split_stim_found()
}

import_clarke2022qjep <- function(data_path) {

  # x and y are already scaled
  read_foraging_csv(file.path(data_path, "clarke2022qjep", "clarke2022qjep_collected.csv"),
                    RT = col_double()) |>
    select(person = observer, condition, trial, id, n = found, v = targ_type, x, y) |>
    mutate(v = as.factor(v + 1),
           condition = fct_recode(as.factor(condition), feature = "1", conjunction = "2")) |>
    fix_person_and_trial() |>
    flip_y() |>
    split_stim_found()
}

import_hughes2024rsos <- function(data_path, scarcity = TRUE) {

  # from the processed csv files (one row per item, one per selection). The
  # condition labels (e.g. "cond_conj_AB") hold both the task (feature or
  # conjunction) and the scarcity level (scarceA, scarceB or equal)
  tidy_hughes <- function(d) {
    d |>
      mutate(condition = condition |>
               str_remove("cond_") |>
               str_replace("conj", "conjunction") |>
               str_replace("_", "_scarce") |>
               str_replace("scarceAB", "equal")) |>
      separate_wider_delim(condition, "_", names = c("condition", "scarcity")) |>
      # number trial_p through the whole session, rather than within each
      # scarcity block
      mutate(condition = as_factor(condition),
             blk = as.numeric(factor(scarcity)),
             trial_p = if_else(trial_p > 5,
                               10 + (blk - 1) * 5 + trial_p,
                               (blk - 1) * 5 + trial_p)) |>
      select(-blk) |>
      arrange(person, trial)
  }

  d <- list(stim  = read_csv(file.path(data_path, "hughes2024rsos", "hughes2024rsos_stim.csv"),
                             show_col_types = FALSE),
            found = read_csv(file.path(data_path, "hughes2024rsos", "hughes2024rsos_found.csv"),
                             show_col_types = FALSE)) |>
    map(tidy_hughes)

  if (scarcity) {
    # combine condition (conjunction/feature) and scarcity (scarceA/equal/
    # scarceB) into a single 6-level condition factor, so prep_stan_list()
    # treats all 6 task/scarcity combinations as separate conditions (K = 6)
    # rather than folding scarcity into a covariate outside the model
    map(d, \(x) x |>
          unite("condition", condition, scarcity, sep = "_") |>
          mutate(condition = as_factor(condition)))
  } else {
    # scarcity = FALSE: just the feature/conjunction comparison - drop the
    # scarcity column entirely (K = 2), pooling the scarceA/equal/scarceB
    # blocks together within each task condition rather than splitting on
    # scarcity at all
    map(d, \(x) select(x, -scarcity))
  }
}

import_kristjansson2014plos <- function(data_path) {

  # printed straight away (a warning() would only show at the end of a
  # script), on every load
  message("Warning: kristjansson2014plos x and y are flipped (180 degree rotation) ",
          "so that first selections start top-left like the other datasets. This ",
          "orientation is an assumption - treat absolute-direction results with caution.")

  read_foraging_csv(file.path(data_path, "kristjansson2014plos", "human_data_2014_prox.csv")) |>
    select(person = observer, condition, trial, id, n = found, v = targ_type, x, y) |>
    mutate(v = as.factor(v + 1),
           condition = fct_recode(as.factor(condition), feature = "1", conjunction = "2"),
           # rotate the display by 180 degrees. As stored, first
           # selections cluster bottom-right, the reverse of every other
           # dataset (top-left); the raw x/y orientation hasn't been checked
           # against the original iPad display code
           x = 1 - x,
           y = 1 - y) |>
    fix_person_and_trial() |>
    split_stim_found()
}

