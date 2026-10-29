# Real-data examples

Each page here works through FoMo fitted to one published dataset. The
fits take hours, so the pages don't fit anything themselves: they read
saved fits. The fits are **not in the repository** (each is 40-200 MB), so
you make them once, and then the pages can run.

The rendered pages are online at <https://riadsala.github.io/FoMo2/>, so
you can read them without any fits.

## How to run a page

Every folder works the same way:

1. Run its `fit_<name>.R` script, from that folder. It makes every fit the
   page needs and saves them in `fits/` in the same folder.
2. Render the page (`quarto render <name>.qmd`), or run its `<name>.R`
   script line by line in RStudio. Both read the fits from `fits/`.

For example, from `examples/real_data/tagu2022cog/`:

```sh
Rscript fit_tagu2022cog.R
quarto render tagu2022cog.qmd
```

| Folder | Fit script | Time on a laptop |
|---|---|---|
| `kristjansson2014plos/` | `fit_kristjansson2014plos.R` (2 fits) | about 2.5 h |
| `tagu2022cog/` | `fit_tagu2022cog.R` (2 fits) | about 2.5 h |
| `hughes2024rsos/` | `fit_hughes2024rsos.R` (2 fits) | about 3.5 h |
| `clarke2022qjep/` | `fit_clarke2022qjep.R` (1 fit) | about 4 h |
| `fourier_versus_vonmises/` | `fit_fourier_versus_vonmises.R` (4 fits, then the comparison); run the four scripts above first | about 14 h |

Run the scripts where the computer won't go to sleep (plugged in, lid
open), or from a terminal with `nohup`.

## The fits

A fit file (`fits/<dataset>_<model>.rds`) holds the whole posterior and
the data list it was made from, so `readRDS()` gives you everything
`fomo_summary()`, `predict_fomo()`, `sim_fomo()` and `fomo_compare()`
need. For example, `kristjansson2014_fourier124.rds` is the default model
(multi-level, Fourier-124 directions) fitted to Kristjánsson et al. (2014).

The fit scripts are also the simplest example of fitting FoMo to your own
data: import, `prep_stan_list()`, `fit_fomo()`.
