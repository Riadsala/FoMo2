##########################################################################
# Fit FoMo to Tagu et al. (2022): every fit tagu2022cog.qmd reads. Run
# this first, from this folder:
#   Rscript fit_tagu2022cog.R
# The fits are saved in fits/ in this folder, where the page reads them.
# About 2.5 hours on a laptop for both fits; run it where the computer
# won't go to sleep.
#
# The minimal pipeline: data -> data list -> fit, with the default model
# (multi-level, Fourier-124 directions, default priors and centring).
##########################################################################

source("../../../functions/set_up_FoMo.R")

# the data: stim (every item on every trial) and found (the selections)
d <- import_data("tagu2022cog")

# the data list for Stan, with the default model: multi-level (one person
# per participant), Fourier-124 directions, the default priors
sl <- prep_stan_list(d$stim, found = d$found)

# fit: 4 chains of 1000 warmup and 1000 sampling iterations, saved as
# fits/tagu2022cog_fourier124.rds
fit <- fit_fomo(sl, iter = 1000, chains = 4, threads_per_chain = 2,
                save_dir = "fits", filename = "tagu2022cog_fourier124")

# for the page's comparison: the same data with two item classes instead
# of three (the high-value colour, and the rest)
d2 <- import_data("tagu2022cog2items")
sl2 <- prep_stan_list(d2$stim, found = d2$found)
fit2 <- fit_fomo(sl2, iter = 1000, chains = 4, threads_per_chain = 2,
                 save_dir = "fits", filename = "tagu2022cog2items_fourier124")
