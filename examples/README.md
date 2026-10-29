# Examples

Worked examples of FoMo, each a rendered page with the code that made it.
Read the pages at <https://riadsala.github.io/FoMo2/>. They are not in the
repository: render a `.qmd` (below) to make the `.html` locally.

## Reading order

1. **`simulated/single_level/spatial/`**: start here. Simulate foraging
   data from the model with known parameters, fit the model, check it
   recovers them, and read the results. The model is explained here.
2. **`simulated/single_level/cat/`**: add item classes (a
   class-transition matrix).
3. **`simulated/multi_level/`** (`spatial`, `cat`): several people,
   each with their own parameters around the population values.
4. **`real_data/<dataset>/`**: a model fitted to a real dataset, from
   checking convergence to posterior predictive checks, most ending with an
   example model comparison.

## Simulated data (`simulated/`)

Small datasets that fit in a few minutes. Each page simulates data, fits
the model, checks parameter recovery (are the true values inside the 97%
HPDIs?), then shows the posterior, prediction accuracy and a posterior
predictive check. All use the default directional model, Fourier-124.

| Page | Shows |
|---|---|
| `single_level/spatial/` | distance, turn and direction terms; the model's equations |
| `single_level/cat/` | item classes and the class-transition matrix |
| `multi_level/spatial/`, `cat/` | per-person effects, shrinkage |
| `model_comparison/fourier_v_vonmises/` | the Fourier and von Mises directional models compared |
| `model_comparison/condition_difference/` | do two conditions differ (a stick bias in one condition, and a null dataset)? the posterior difference, and separate vs pooled conditions on loo |

## Real data (`real_data/`)

| Page | Dataset | Comparison |
|---|---|---|
| `kristjansson2014plos/` | Kristjánsson et al. (2014): feature and conjunction foraging | full vs symmetric class-transition matrix |
| `tagu2022cog/` | Tagu et al. (2022): three target colours, with and without value | three item classes vs two |
| `hughes2024rsos/` | Hughes et al. (2024): target scarcity | does scarcity matter? (individuals change, the average doesn't) |
| `clarke2022qjep/` | Clarke et al. (2022): 58 people, feature and conjunction | none yet |
| `fourier_versus_vonmises/` | every dataset | Fourier vs von Mises directional models |

These pages read saved fits rather than fitting live: the fits take hours
and are not in the repository. In each folder, run `fit_<name>.R` first:
it makes every fit the page needs and saves them in `fits/` in that
folder, where the page reads them (see `real_data/README.md`). Without
the fits a page won't run, but its online version still shows the results.

## Running the examples

Each `.qmd` page runs its `.R` script chunk by chunk, so the `.R` file is
the code to read or copy. To re-render a page, from its folder:

```sh
quarto render spatial.qmd
```

- Slow steps (fits, posterior simulations) are cached, so a second render
  takes seconds.
- Everything is seeded in `functions/set_up_FoMo.R`.
- To update the online pages after re-rendering, run
  `sh examples/publish_pages.sh` (it replaces the `gh-pages` branch with
  the current pages, so their history isn't kept).
