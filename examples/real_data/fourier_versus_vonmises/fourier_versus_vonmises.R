##########################################################################
# Fourier-124 against von Mises directional models, on every real
# dataset. Reads the fits and the comparison that
# fit_fourier_versus_vonmises.R saves (and each dataset's Fourier-124 fit,
# in ../<dataset>/fits/); nothing is refitted.
# fourier_versus_vonmises.qmd renders this script; the "## ----" lines
# mark its chunks. The directional curves (reading every fit) are cached
# by the .qmd.
##########################################################################

## ---- load
source("../../../functions/set_up_FoMo.R")

# the four datasets: the name used in the fit files, and each dataset's
# folder (where its Fourier-124 fit is)
datasets <- c(kristjansson2014 = "kristjansson2014plos",
              hughes2024rsos   = "hughes2024rsos",
              clarke2022qjep   = "clarke2022qjep",
              tagu2022cog      = "tagu2022cog")
models <- c(fourier124 = "Fourier-124", vonmises = "von Mises")

fourier_file <- function(name) file.path("..", datasets[[name]], "fits", paste0(name, "_fourier124.rds"))
vonmises_file <- function(name) file.path("fits", paste0(name, "_vonmises.rds"))

if (!file.exists("fits/comparison.csv")) {
  stop("Run fit_fourier_versus_vonmises.R first: it makes the fits and the comparison this page reads")
}

# each dataset as its page fits it (hughes2024rsos: the two tasks, with the
# scarcity levels pooled)
import_example <- function(name) {
  if (name == "hughes2024rsos") return(import_data("hughes2024rsos", scarcity = FALSE))
  import_data(datasets[[name]])
}

# fomo_compare() for every dataset
comparison <- read_csv("fits/comparison.csv", show_col_types = FALSE) |>
  mutate(dataset = factor(dataset, levels = names(datasets)),
         model = factor(models[model], levels = models)) |>
  arrange(dataset, model)

## ---- table
format_compare(comparison, group = "dataset",
               columns = c("dataset", "model", "deterministic", "nn", "gain_nn", "elpd_diff",
                           "se_diff", "pct_high_k", "robust"))

## ---- fig-elpd
# the loo difference per trial, so datasets of different sizes are on one
# scale (loo leaves out one trial at a time)
n_trials <- map_int(set_names(names(datasets)), \(name) n_distinct(import_example(name)$found$trial))

comparison |>
  filter(model == "von Mises") |>
  mutate(n_trials = n_trials[as.character(dataset)],
         per_trial = elpd_diff / n_trials,
         se = se_diff / n_trials) |>
  ggplot(aes(per_trial, fct_rev(dataset))) +
  theme_fomo() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_pointrange(aes(xmin = per_trial - 2 * se, xmax = per_trial + 2 * se),
                  colour = MODEL_COLOURS[2]) +
  labs(x = "elpd difference per trial (von Mises - Fourier-124)", y = NULL)

## ---- fig-accuracy
# deterministic accuracy (the model's most likely item is the one
# selected), with the nearest-neighbour rule on the same selections
comparison |>
  select(dataset, model, deterministic) |>
  bind_rows(distinct(comparison, dataset, deterministic = nn) |> mutate(model = "nearest neighbour")) |>
  mutate(model = factor(model, levels = c(models, "nearest neighbour"))) |>
  ggplot(aes(deterministic, fct_rev(dataset), colour = model, shape = model)) +
  theme_fomo() +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c(MODEL_COLOURS[1:2], "grey50")) +
  labs(x = "accuracy", y = NULL, colour = NULL, shape = NULL)

## ---- fig-stacking
# stacking weights: the weights of the mixture of the two models' loo
# predictive distributions that predicts held-out trials best
comparison |>
  ggplot(aes(stacking, fct_rev(dataset), fill = model)) +
  theme_fomo() +
  geom_col(width = 0.7) +
  scale_fill_manual(values = MODEL_COLOURS[1:2]) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "stacking weight", y = NULL, fill = NULL)

## ---- curves
# each model's directional log-weight curve (centred), from 200 posterior
# draws of each fit: only the directional parameters are read
directional_curves <- function(name, model) {
  fit <- readRDS(if (model == "fourier124") fourier_file(name) else vonmises_file(name))
  cond_labels <- levels(import_example(name)$found$condition)
  # the names from the metadata: get_variables() turns every draw into a
  # data frame first, too big for large fits
  vars <- fit$metadata()$model_params
  post <- list(vonmises = extract_vonmises(fit, vars, cond_labels),
               fourier = extract_fourier(fit, vars, cond_labels))
  out <- compute_directional_curve(post, draws = sample(posterior::ndraws(fit$draws("lp__")), 200)) |>
    mutate(dataset = name, model = models[model])
  # each fit holds all its draws (GBs): free it before reading the next, or
  # the 8 fits reach R's memory limit
  rm(fit, post)
  gc(verbose = FALSE)
  out
}

curves <- crossing(name = names(datasets), model = names(models)) |>
  pmap(directional_curves) |>
  list_rbind()

## ---- fig-curves
curves |>
  mutate(condition = paste0(dataset, ": ", condition),
         model = factor(model, levels = models)) |>
  plot_post_directional_curves() +
  facet_wrap(~condition, ncol = 3)
