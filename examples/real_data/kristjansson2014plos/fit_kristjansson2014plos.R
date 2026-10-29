##########################################################################
# Fit FoMo to Kristjansson, Johannesson & Thornton (2014): every fit
# kristjansson2014plos.qmd reads. Run this first, from this folder:
#   Rscript fit_kristjansson2014plos.R
# The fits are saved in fits/ in this folder, where the page reads them.
# About 2.5 hours on a laptop for both fits; run it where the computer
# won't go to sleep.
#
# The minimal pipeline: data -> data list -> fit, with the default model
# (multi-level, Fourier-124 directions, default priors and centring).
##########################################################################

source("../../../functions/set_up_FoMo.R")

# the data: stim (every item on every trial) and found (the selections)
d <- import_data("kristjansson2014plos")

# the data list for Stan, with the default model: multi-level (one person
# per participant), Fourier-124 directions, the default priors
sl <- prep_stan_list(d$stim, found = d$found)

# fit: 4 chains of 1000 warmup and 1000 sampling iterations, saved as
# fits/kristjansson2014_fourier124.rds
fit <- fit_fomo(sl, iter = 1000, chains = 4, threads_per_chain = 2,
                save_dir = "fits", filename = "kristjansson2014_fourier124")

# for the page's comparison: the same model with a symmetric transition
# matrix (one stay probability per condition)
sl_sym <- prep_stan_list(d$stim, found = d$found, symmetric_v = TRUE)
fit_sym <- fit_fomo(sl_sym, iter = 1000, chains = 4, threads_per_chain = 2,
                    save_dir = "fits", filename = "kristjansson2014_fourier124_symmetric")
