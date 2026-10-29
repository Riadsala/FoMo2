##########################################################################
# Fit FoMo to Clarke, Irons, James, Leber & Hunt (2022): every fit
# clarke2022qjep.qmd reads. Run this first, from this folder:
#   Rscript fit_clarke2022qjep.R
# The fit is saved in fits/ in this folder, where the page reads it.
# About 4 hours on a laptop; run it where the computer won't go to sleep.
#
# The minimal pipeline: data -> data list -> fit, with the default model
# (multi-level, Fourier-124 directions, default priors and centring).
##########################################################################

source("../../../functions/set_up_FoMo.R")

# the data: stim (every item on every trial) and found (the selections)
d <- import_data("clarke2022qjep")

# the data list for Stan, with the default model: multi-level (one person
# per participant), Fourier-124 directions, the default priors
sl <- prep_stan_list(d$stim, found = d$found)

# fit: 4 chains of 1000 warmup and 1000 sampling iterations, saved as
# fits/clarke2022qjep_fourier124.rds
fit <- fit_fomo(sl, iter = 1000, chains = 4, threads_per_chain = 2,
                save_dir = "fits", filename = "clarke2022qjep_fourier124")
