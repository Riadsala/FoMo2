##########################################################################
# Checks that every dataset importer in functions/import_data.R produces
# data in the format prep_stan_list() relies on (see data/README.md; checked
# by verify_data(), which import_data() runs), and that prep_stan_list()
# turns it into a consistent Stan data list for both directional models.
# verify_data() itself is tested by test_verify_data.R. No models are
# compiled or fitted.
#
#   Rscript tests/test_dataset_import.R
#
# Datasets whose data files aren't present are reported as SKIPPED, not
# failures. Exits with status 1 if any dataset fails.
##########################################################################

root <- rprojroot::find_root(rprojroot::has_dir("models"))
setwd(file.path(root, "tests"))
source("../functions/set_up_FoMo.R")

# every importer in import_data(), plus both Hughes variants
cases <- list(
  list(name = "clarke2022qjep"),
  list(name = "tagu2022cog2items"),
  list(name = "tagu2022cog"),
  list(name = "kristjansson2014plos"),
  list(name = "hughes2024rsos", args = list(scarcity = TRUE)),
  list(name = "hughes2024rsos", args = list(scarcity = FALSE)))

check <- function(cond, msg, problems) if (!isTRUE(cond)) c(problems, msg) else problems

check_dataset <- function(d) {

  # the data format itself is checked by verify_data(), which import_data()
  # runs (a breach shows up as an import_data() error); here, that
  # prep_stan_list() turns it into a consistent Stan data list
  stim <- d$stim
  found <- d$found
  p <- character(0)
  T_ <- n_distinct(stim$trial)
  K <- nlevels(stim$condition)
  L <- n_distinct(stim$person)

  # --- prep_stan_list, both directional models
  for (dm in c("von_mises", "fourier")) {
    sl <- tryCatch(suppressMessages(prep_stan_list(stim, found = found, directional_model = dm)),
                   error = function(e) { p <<- c(p, paste0("prep_stan_list(", dm, ") failed: ", conditionMessage(e))); NULL })
    if (is.null(sl)) next
    p <- check(sl$T == T_ && length(sl$s) == T_ && length(sl$f) == T_ &&
                 length(sl$X) == T_ && length(sl$z) == T_,
               paste0(dm, ": T / s / f / X / z lengths disagree"), p)
    p <- check(sum(sl$s) == sl$N, paste0(dm, ": sum(s) != N"), p)
    p <- check(sl$Ns == nrow(found) && length(sl$Y) == sl$Ns, paste0(dm, ": Ns / Y disagree with found"), p)
    p <- check(all(sl$f <= sl$s), paste0(dm, ": more selections than items in some trial"), p)
    p <- check(sl$K == K && sl$L == L, paste0(dm, ": K or L wrong"), p)
    p <- check(all(is.finite(sl$x)) && all(is.finite(sl$y)), paste0(dm, ": non-finite x/y"), p)
  }

  p
}

results <- list()

for (cs in cases) {
  label <- paste0(cs$name, if (!is.null(cs$args)) paste0(" (", paste(names(cs$args), cs$args, sep = "=", collapse = ", "), ")") else "")

  d <- tryCatch(suppressMessages(do.call(import_data, c(list(cs$name), cs$args))),
                error = function(e) e)

  if (inherits(d, "error")) {
    if (grepl("does\\s+not\\s+exist", conditionMessage(d))) {
      status <- "SKIPPED"; detail <- "data file not present"
    } else {
      status <- "FAILED"; detail <- paste("import_data() error:", conditionMessage(d))
    }
  } else {
    problems <- check_dataset(d)
    status <- if (length(problems) == 0) "OK" else "FAILED"
    detail <- if (length(problems) == 0) {
      sprintf("K=%d L=%d T=%d N=%d Ns=%d", nlevels(d$stim$condition), length(unique(d$stim$person)),
              length(unique(d$stim$trial)), nrow(d$stim), nrow(d$found))
    } else paste(problems, collapse = "; ")
  }

  cat(sprintf("%-38s %-8s %s\n", label, status, detail))
  results[[label]] <- status
}

n_fail <- sum(unlist(results) == "FAILED")
n_skip <- sum(unlist(results) == "SKIPPED")
cat(sprintf("\n%d OK, %d SKIPPED, %d FAILED\n", sum(unlist(results) == "OK"), n_skip, n_fail))
if (n_fail > 0) quit(status = 1)
