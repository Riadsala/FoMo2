# FoMo 2.0

A Bayesian model of visual foraging that predicts each target selection from
proximity, movement direction and item features. R + Stan (cmdstanr)
implementation of Clarke & Hughes (2026, *PLoS Comp Bio*), building on Clarke,
Hunt & Hughes (2022, *PLoS Comp Bio*).

## Installation

FoMo is a set of R functions and Stan models, not a package: clone or
download this repository and run scripts from inside it.

1. R (4.1 or later; developed on 4.6) and the R packages:

   ```r
   install.packages(c("tidyverse", "tidybayes", "patchwork", "paletteer", "loo",
                      "posterior", "bayesplot", "ggdist", "kableExtra", "withr",
                      "rprojroot", "knitr", "rlang", "scales", "circular",
                      "ForagingOrg"))
   ```

   ggplot2 needs to be version 4.0 or later (it comes with the tidyverse).
2. cmdstanr and CmdStan (version 2.39 or later), which compile and run
   the Stan models:

   ```r
   install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
   cmdstanr::check_cmdstan_toolchain(fix = TRUE)
   cmdstanr::install_cmdstan()
   ```

   See the [cmdstanr guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)
   if the toolchain check fails (on Windows, it needs RTools).
3. Optional: [Quarto](https://quarto.org), to render the example pages
   (`.qmd`) yourself. The rendered pages are online at
   <https://riadsala.github.io/FoMo2/>.

Then start any script with `source("<path to>/functions/set_up_FoMo.R")`.
Each model compiles the first time it is used (a minute or two each).

## Start here

1. [`examples/simulated/single_level/spatial/`](https://riadsala.github.io/FoMo2/examples/simulated/single_level/spatial/spatial.html): simulate data
   from the model, fit it, check the parameters are recovered, and read
   the results. The model is explained there.
2. The other simulated examples (`single_level/cat`, then `multi_level/`):
   item classes and per-person effects.
3. `examples/real_data/<dataset>/`: working with a model fitted to real
   data. These pages read saved fits, which you make once with each
   page's `fit_<dataset>.R` (see `examples/real_data/README.md`).

`examples/README.md` lists every page.

## Models

Every model predicts each selection from the distance to each remaining
item (`rho_delta`), the turn relative to the previous step (`rho_psi`) and
the absolute direction of the step, plus, in the `cat` models, the items'
classes.

| Model | Item features | Single-level | Multi-level | Selected with |
|---|---|---|---|---|
| `spatial` | none (position only) | ✓ | ✓ | data without `v` |
| `cat` | item classes: a class-transition matrix (`b_v`); each trial's first selection is weighted by the matrix's stationary distribution | ✓ | ✓ | `v` a factor |
| `cat_sym` | as `cat`, with a symmetric transition matrix: one stay probability per condition (for model comparison; Fourier, fit only) | – | ✓ | `prep_stan_list(symmetric_v = TRUE)` |

- **Directional term:** `fourier` (the default, Fourier-124: harmonics 1,
  2 and 4) or `von_mises` (the papers' model), chosen with
  `prep_stan_list(directional_model = )`. Fourier-124 fits better on every
  real dataset, and faster.
- **Single- or multi-level:** the multi-level version (per-person random
  effects) is used automatically when the data contains more than one
  person. By default it samples the distance and class effects centred
  and the rest non-centred (`prep_stan_list(centred_params = )`).
- **Fit and simulate:** each model except `cat_sym` has a fit version and
  a simulate version (for simulating data and posterior predictive
  checks).

In all: 17 Stan models in `models/` (2 item-feature types x 2
directional terms x single- or multi-level, each with fit and simulate
versions, plus `cat_sym`); shared functions in `models/include/`.

## Usage

Scripts can live anywhere inside the repository: the setup and import
functions find `functions/` and `data/` from the repository root. From a
folder two levels down:

```r
source("../../functions/set_up_FoMo.R")

d    <- import_data("kristjansson2014plos")      # or the shorthand "K2014"
sl   <- prep_stan_list(d$stim, found = d$found)   # Fourier-124, multi-level
fit  <- fit_fomo(sl, iter = 1000)
post <- extract_post(fit, d$found, d$stim)        # tidy posterior tables
plot_post(post)

pred <- predict_fomo(fit, sl, n_draws = 500)      # one-step-ahead predictions
acc  <- extract_post(pred, d$found, d$stim)$acc   # accuracy
```

`fomo_compare()` compares models fitted to the same data (loo, accuracy,
model weights). To use your own data, put it in the format described in
`data/README.md` (one table of items, one of selections).

*Note*: FoMo's functions call dplyr verbs such as `select()`, `mutate()`
and `summarise()` without a package prefix, so a package attached after
`set_up_FoMo.R` can replace them. Attach any other packages first, and
call plyr and MASS with a prefix (`plyr::`, `MASS::`) rather than
attaching them: MASS's `select()` makes FoMo's functions fail, and plyr's
`summarise()` ignores `group_by()`, giving wrong results without an error.

## Datasets

`kristjansson2014plos` (K2014), `hughes2024rsos` (H2024), `clarke2022qjep`
(C2022) and `tagu2022cog` (T2022; `tagu2022cog2items` with 2 classes).
Data format, caveats and sources: `data/README.md`.

## Repository

| Folder | Contents |
|---|---|
| `models/` | Stan models: `{single,multi}_level/{fit,simulate}/{von_mises,fourier}/<model>.stan`; shared functions in `include/` |
| `functions/` | R code: data import and checks, data prep, fitting, extraction, summaries, plotting; default priors |
| `data/` | Datasets; data format in `data/README.md` |
| `examples/` | Worked examples. `simulated/`: `{single,multi}_level/{spatial,cat}/` and `model_comparison/`; `real_data/`: one folder per dataset (a fit script, then the walkthrough), and the Fourier vs von Mises comparison |
| `tests/` | Model, data, prediction, comparison and k-fold tests |

## Tests

```sh
Rscript tests/smoke_test_models.R    # after changing models or fitting code
Rscript tests/test_dataset_import.R  # after changing importers or data prep
Rscript tests/test_predict_fomo.R    # after changing predictions or accuracy
```

## References

- Clarke, A. D. F., & Hughes, A. E. (2026). FoMo: A unifying theory of
  visual foraging. *PLOS Computational Biology*, 22(5), e1014266.
  doi:10.1371/journal.pcbi.1014266. The model this repository implements.
- Clarke, A. D. F., Hunt, A. R., & Hughes, A. E. (2022a). Foraging as
  sampling without replacement: A Bayesian statistical model for
  estimating biases in target selection. *PLoS Computational Biology*, 18,
  e1009813. doi:10.1371/journal.pcbi.1009813
- Clarke, A. D. F., Hunt, A. R., & Hughes, A. E. (2022b). A Bayesian
  statistical model is able to predict target-by-target selection
  behaviour in a human foraging task. *Vision*, 6(4), 66.
  doi:10.3390/vision6040066

The datasets' papers are listed in `data/README.md`.

## Citing FoMo

Please cite Clarke & Hughes (2026), above. `CITATION.cff` has the details
(GitHub's "Cite this repository" button reads it).

## Licence

The code (`functions/`, `models/`, `tests/` and the example scripts) is
under the MIT licence (`LICENSE`). The text and figures of the example
pages and documentation, and the `hughes2024rsos` data, are under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The other
datasets keep their original licences: see `data/README.md`.
