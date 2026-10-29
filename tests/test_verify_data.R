##########################################################################
# Checks that verify_data() (functions/import_data.R) catches each way the
# data can break the format in data/README.md: every case below breaks a
# clean dataset (kristjansson2014plos) in one way, and verify_data() must
# stop with the expected message. The clean data and the note-only cases (allowed,
# but worth knowing) must pass. Also: stim alone passes (simulating), and
# prep_stan_list() stops on broken data, since it runs verify_data() too.
# No models are compiled or fitted.
#
#   Rscript tests/test_verify_data.R
#
# Exits with status 1 if any case fails.
##########################################################################

root <- rprojroot::find_root(rprojroot::has_dir("models"))
setwd(file.path(root, "tests"))
source("../functions/set_up_FoMo.R")

d <- suppressMessages(import_data("kristjansson2014plos"))
s <- d$stim
f <- d$found

# drop one person x condition, renumbering the trials so nothing else breaks
drop_cell <- function(x) {
  x |>
    filter(!(person == 1 & condition == "feature")) |>
    mutate(trial = consecutive_id(trial))
}

# name, stim, found, and the expected error (a regex), or NA to pass
cases <- tribble(
  ~name,                         ~stim,                                          ~found,                                   ~expect,
  "clean",                       s,                                              f,                                        NA,
  "note: person misses a condition", drop_cell(s),                               drop_cell(f),                             NA,
  "note: one-selection trial",   s,                                              filter(f, !(trial == 2 & n > 1)),         NA,
  "missing columns",             select(s, -id),                                 select(f, -n),                            "missing columns: stim\\$id, found\\$n",
  "NA in x",                     mutate(s, x = if_else(row_number() == 10, NA, x)), f,                                     "NAs in stim",
  "condition not a factor",      mutate(s, condition = as.character(condition)), f,                                        "not a factor",
  "person numbering gap",        mutate(s, person = if_else(person == 16, 17, person)),
                                 mutate(f, person = if_else(person == 16, 17, person)),                                    "1\\.\\.L with no gaps",
  "stim not sorted",             arrange(s, desc(trial)),                        f,                                        "stim is not sorted",
  "stim ids out of order",       mutate(s, id = if_else(trial == 1, rev(id), id), .by = trial), f,                         "1\\.\\.nt in row order",
  "x outside [0, 1]",            mutate(s, x = x + 1),                           mutate(f, x = x + 1),                     "outside \\[0, 1\\]",
  "empty trial",                 s,                                              filter(f, trial != 5),                    "no selections",
  "found not sorted",            s,                                              arrange(f, desc(trial)),                  "found is not sorted",
  "found$n has a gap",           s,                                              mutate(f, n = if_else(trial == 1 & n == 3, 99, n)), "found\\$n is not 1\\.\\.f",
  "item selected twice",         s,                                              mutate(f, id = if_else(trial == 3 & n == 2, id[trial == 3 & n == 1], id)),
                                                                                                                           "selected more than once",
  "selected item doesn't exist", s,                                              mutate(f, id = if_else(trial == 1 & n == 1, 999, id)), "found\\$id outside",
  "found person disagrees",      s,                                              mutate(f, person = if_else(trial == 1, 2, person)), "found\\$person disagrees",
  "found condition levels",      s,                                              mutate(f, condition = fct_rev(condition)), "different levels",
  "cat v codes too high",        mutate(s, v = factor(v, levels = c("0", "1", "2"))),
                                 mutate(f, v = factor(v, levels = c("0", "1", "2"))),                                      "factor codes above n_classes",
  "v not a factor",              mutate(s, v = as.numeric(v)),                   f,                                        "stim\\$v is not a factor",
  "found v disagrees",           s,                                              mutate(f, v = factor(levels(v)[3 - as.integer(v)], levels(v))), "found\\$v disagrees")

n_fail <- 0
for (i in seq_len(nrow(cases))) {
  cs <- cases[i, ]
  err <- tryCatch({
    suppressMessages(verify_data(cs$stim[[1]], cs$found[[1]], cs$name))
    NA_character_
  }, error = conditionMessage)

  ok <- if (is.na(cs$expect)) is.na(err) else !is.na(err) && grepl(cs$expect, err)
  detail <- if (ok) "" else if (is.na(err)) "passed, but should have stopped" else err
  cat(sprintf("%-34s %-6s %s\n", cs$name, if (ok) "OK" else "FAILED", detail))
  n_fail <- n_fail + !ok
}

# prep_stan_list() runs verify_data() (found optional, no notes)
extra <- list(
  "stim only (to simulate)" = function() verify_data(s, label = "stim only"),
  "prep_stan_list(): clean data" = function() suppressMessages(prep_stan_list(s, found = f)),
  "prep_stan_list(): stim only" = function() suppressMessages(prep_stan_list(s)))
for (name in names(extra)) {
  err <- tryCatch({ extra[[name]](); NA_character_ }, error = conditionMessage)
  ok <- is.na(err)
  cat(sprintf("%-34s %-6s %s\n", name, if (ok) "OK" else "FAILED", if (ok) "" else err))
  n_fail <- n_fail + !ok
}
err <- tryCatch({ prep_stan_list(s, found = filter(f, trial != 5)); NA_character_ }, error = conditionMessage)
ok <- !is.na(err) && grepl("prep_stan_list\\(\\) fails verify_data", err) && grepl("no selections", err)
cat(sprintf("%-34s %-6s %s\n", "prep_stan_list(): empty trial", if (ok) "OK" else "FAILED", if (ok) "" else err))
n_fail <- n_fail + !ok

cat(sprintf("\n%d OK, %d FAILED\n", nrow(cases) + length(extra) + 1 - n_fail, n_fail))
if (n_fail > 0) quit(status = 1)
