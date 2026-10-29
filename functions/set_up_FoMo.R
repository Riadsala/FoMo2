# import required packages and load all of our functions

# quietly: the attach/masking messages repeat on every run. Also needed,
# called with their package prefix (so not attached): circular,
# ForagingOrg, kableExtra (fomo_summary()'s tables in rendered pages),
# bayesplot, ggdist, posterior, withr, rprojroot, knitr
suppressPackageStartupMessages({
  library(tidyverse)
  library(cmdstanr)
  library(tidybayes)
  library(patchwork)
  library(paletteer)
  library(loo)
})

functions_dir <- file.path(rprojroot::find_root(rprojroot::has_dir("models")), "functions")

# one file per stage of the workflow: import_data() -> prep_stan_list() ->
# fit_fomo() / sim_fomo() -> extract_post() -> summaries and plots
for (f in c("import_data.R", "prep.R", "simulate.R", "fit.R", "extract.R",
            "summarise.R", "behaviour.R", "plot.R", "maths.R")) {
  source(file.path(functions_dir, f))
}

options(digits = 3, mc.cores = 4)

# one seed for everything that sources this file (examples, tests, analyses),
# so simulated data are the same every run. cmdstanr takes each sampler's
# seed from R's generator when none is given, so fits follow it too
FOMO_SEED <- 65
set.seed(FOMO_SEED)
print(paste0("set_up_FoMo.R: random seed set to ", FOMO_SEED, " (FOMO_SEED)"))
