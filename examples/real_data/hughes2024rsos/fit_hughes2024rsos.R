##########################################################################
# Fit FoMo to Hughes, Statham & Clarke (2024): every fit hughes2024rsos.qmd
# reads. Run this first, from this folder:
#   Rscript fit_hughes2024rsos.R
# The fits are saved in fits/ in this folder, where the page reads them.
# About 3.5 hours on a laptop for both fits; run it where the computer
# won't go to sleep.
#
# The minimal pipeline: data -> data list -> fit, with the default model
# (multi-level, Fourier-124 directions, default priors and centring).
##########################################################################

source("../../../functions/set_up_FoMo.R")

# the data: stim (every item on every trial) and found (the selections):
# feature and conjunction foraging (K = 2), pooled over the scarcity levels
d <- import_data("hughes2024rsos", scarcity = FALSE)

# the data list for Stan, with the default model: multi-level (one person
# per participant), Fourier-124 directions, the default priors
sl <- prep_stan_list(d$stim, found = d$found)

# fit: 4 chains of 1000 warmup and 1000 sampling iterations, saved as
# fits/hughes2024rsos_fourier124.rds
fit <- fit_fomo(sl, iter = 1000, chains = 4, threads_per_chain = 2,
                save_dir = "fits", filename = "hughes2024rsos_fourier124")

# for the page's comparison: task x scarcity as six conditions (K = 6)
d6 <- import_data("hughes2024rsos", scarcity = TRUE)
sl6 <- prep_stan_list(d6$stim, found = d6$found)
fit6 <- fit_fomo(sl6, iter = 1000, chains = 4, threads_per_chain = 2,
                 save_dir = "fits", filename = "hughes2024rsos_scarcity_fourier124")
