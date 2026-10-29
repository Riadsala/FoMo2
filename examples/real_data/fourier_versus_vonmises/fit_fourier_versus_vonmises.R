##########################################################################
# The von Mises fits and the model comparison that
# fourier_versus_vonmises.qmd reads. Run the four dataset scripts first
# (../<dataset>/fit_<dataset>.R: they make the Fourier-124 fits), then
# this one, from this folder:
#   Rscript fit_fourier_versus_vonmises.R
# It fits the von Mises model to every dataset (about 14 hours on a
# laptop in all; run it where the computer won't go to sleep), then
# compares it with each dataset's Fourier-124 fit with fomo_compare().
# Everything is saved in fits/ in this folder. A fit already in fits/ is
# not redone, so a rerun carries on where the last one stopped.
##########################################################################

source("../../../functions/set_up_FoMo.R")

# the four datasets: the name used in the fit files, and each dataset's
# folder (where its Fourier-124 fit is)
datasets <- c(kristjansson2014 = "kristjansson2014plos",
              hughes2024rsos   = "hughes2024rsos",
              clarke2022qjep   = "clarke2022qjep",
              tagu2022cog      = "tagu2022cog")

# each dataset as its page fits it (hughes2024rsos: the two tasks, with the
# scarcity levels pooled)
import_example <- function(name) {
  if (name == "hughes2024rsos") return(import_data("hughes2024rsos", scarcity = FALSE))
  import_data(datasets[[name]])
}

fourier_file <- function(name) file.path("..", datasets[[name]], "fits", paste0(name, "_fourier124.rds"))
vonmises_file <- function(name) file.path("fits", paste0(name, "_vonmises.rds"))

missing <- names(datasets)[!file.exists(map_chr(names(datasets), fourier_file))]
if (length(missing) > 0) stop("Fourier-124 fits missing for ", paste(missing, collapse = ", "),
                              ": run each dataset's fit_<dataset>.R first")

# 1. the von Mises fits, one dataset at a time (the same model otherwise:
# default priors and centring, 4 chains of 1000 warmup and 1000 sampling
# iterations)
for (name in names(datasets)) {
  if (file.exists(vonmises_file(name))) next
  d <- import_example(name)
  sl <- prep_stan_list(d$stim, found = d$found, directional_model = "von_mises")
  fit <- fit_fomo(sl, iter = 1000, chains = 4, threads_per_chain = 2,
                  save_dir = "fits", filename = paste0(name, "_vonmises"))
  rm(fit); gc(verbose = FALSE)
}

# 2. fomo_compare() on each dataset: accuracy, loo and model weights. Each
# fit holds all its draws, so read two at a time and free them after
comparison <- map(names(datasets), function(name) {
  d <- import_example(name)
  fits <- list(fourier124 = readRDS(fourier_file(name)), vonmises = readRDS(vonmises_file(name)))
  out <- mutate(fomo_compare(fits, d), dataset = name, .before = 1)
  rm(fits); gc(verbose = FALSE)
  out
}) |>
  list_rbind()

write_csv(comparison, "fits/comparison.csv")
