##########################################################################
# Checks that simulated and subsetted data meet the format in
# data/README.md. Only import_data() runs verify_data() itself, so data
# built any other way are checked here: each item-feature type
# (spatial, cat) is simulated with generate_stim() -> prep_stan_list()
# -> sim_fomo() -> extract_pred() (a prior predictive draw), then
#   - the simulated data (exhaustive: every item selected)
#   - truncate_trials() of it (early stopping)
#   - filter_one_person() of it (multi-level only)
# must pass verify_data(), and prep_stan_list(stim, found = found) must give
# consistent trial counts. filter_one_person() of a real dataset is checked
# too. Compiles the von Mises simulate models if needed; no fitting.
#
#   Rscript tests/test_simulation_data.R
#
# Exits with status 1 if any check fails.
##########################################################################

root <- rprojroot::find_root(rprojroot::has_dir("models"))
setwd(file.path(root, "tests"))
source("../functions/set_up_FoMo.R")

# verify_data() plus prep_stan_list()'s trial counts; "" if all is well
check <- function(stim, found) {
  tryCatch({
    suppressMessages(verify_data(stim, found))
    sl <- suppressMessages(prep_stan_list(stim, found = found))
    stopifnot("prep_stan_list(): f and s don't have one entry per trial" =
                length(sl$f) == sl$T && length(sl$s) == sl$T,
              "prep_stan_list(): sum(f) != Ns" = sum(sl$f) == sl$Ns,
              "prep_stan_list(): a trial has more selections than items" = all(sl$f <= sl$s))
    ""
  }, error = conditionMessage)
}

cases <- list(
  list(name = "spatial, single-level", n_people = 1,
       spec = list(K = 2, nrow = 5, ncol = 5, prop_missing = 0.2, v_type = "none")),
  list(name = "cat, multi-level", n_people = 3,
       spec = list(K = 2, nrow = 5, ncol = 5, prop_missing = 0.2, v_type = "cat", n_classes = 2)))

results <- list()
record <- function(label, err) {
  cat(sprintf("%-48s %-6s %s\n", label, if (err == "") "OK" else "FAILED", err))
  results[[label]] <<- err == ""
}

for (cs in cases) {

  stim <- generate_stim(cs$spec, n_trials = 3, n_people = cs$n_people)
  found <- tryCatch({
    capture.output(sl <- suppressMessages(prep_stan_list(stim)))
    capture.output(m <- suppressMessages(sim_fomo(sl, iter = 1, chains = 1)))
    extract_pred(m, stim)
  }, error = conditionMessage)

  if (is.character(found)) {
    record(paste0(cs$name, ": simulate"), found)
    next
  }

  record(paste0(cs$name, ": simulated"), check(stim, found))
  record(paste0(cs$name, ": truncate_trials()"),
         check(stim, truncate_trials(found, lambda = 5, min_keep = 3)))

  if (cs$n_people > 1) {
    one <- filter_one_person(list(stim = stim, found = found), 2)
    record(paste0(cs$name, ": filter_one_person()"), check(one$stim, one$found))
  }
}

d <- suppressMessages(import_data("kristjansson2014plos"))
one <- filter_one_person(d, 3)
record("kristjansson2014plos: filter_one_person()", check(one$stim, one$found))

n_fail <- sum(!unlist(results))
cat(sprintf("\n%d OK, %d FAILED\n", length(results) - n_fail, n_fail))
if (n_fail > 0) quit(status = 1)
