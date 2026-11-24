##########################################################################
# Fit FoMo to Haggouch, Doré-Mazars & Tagu (2026): the fit
# haggouch2026acta.qmd reads. Run this first, from this folder:
#   Rscript fit_haggouch2026acta.R
# The fit is saved in fits/ in this folder, where the page reads it.
# By far the slowest of the real-data examples: 304 people, 8512 trials of
# 40 selections (340,480 selections, about 4 times C2022). Not yet timed:
# from benchmarks, about 5 hours on a laptop, with about 1.5 GB of memory
# per chain (6 GB for the 4). Run it where the computer won't go to sleep.
#
# The minimal pipeline: data -> data list -> fit, with the default model
# (multi-level, Fourier-124 directions, default priors and centring).
##########################################################################

source("../../../functions/set_up_FoMo.R")

# the data: stim (every item on every trial) and found (the selections)
d <- import_data("haggouch2026acta")

# the data list for Stan, with the default model: multi-level (one person
# per participant), Fourier-124 directions, the default priors
sl <- prep_stan_list(d$stim, found = d$found)

# fit: 4 chains of 1000 warmup and 1000 sampling iterations, saved as
# fits/haggouch2026acta_fourier124.rds
fit <- fit_fomo(sl, iter = 1000, chains = 4, threads_per_chain = 2,
                save_dir = "fits", filename = "haggouch2026acta_fourier124")
