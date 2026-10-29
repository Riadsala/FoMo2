# Data format

Every importer in `functions/import_data.R` returns
`list(dataset, stim, found)`, and `prep_stan_list(stim, found = found)` turns
that into Stan data. The models rely on the rules below: breaking one would
usually give silently misaligned data rather than an error, so
`verify_data(stim, found)` (in `functions/import_data.R`) checks them and
stops on any breach. Both `import_data()` and `prep_stan_list()` run it (the
first also prints notes on allowed oddities), so data built any other way,
e.g. simulated or subsetted, is checked too.
`tests/test_dataset_import.R` also checks them for every dataset.

## Sources

| Dataset | Paper | Data source | Licence |
|---|---|---|---|
| `clarke2022qjep` | Clarke, Irons, James, Leber & Hunt (2022). Stable individual differences in strategies within, but not between, visual search tasks. *Quarterly Journal of Experimental Psychology*, 75(2), 289-296. doi:10.1177/1747021820929190 | [osf.io/y6qbv](https://osf.io/y6qbv/) | Academic Free License (AFL) 3.0 |
| `kristjansson2014plos` | Kristjánsson, Jóhannesson & Thornton (2014). Common attentional constraints in visual foraging. *PLoS ONE*, 9(6), e100752. doi:10.1371/journal.pone.0100752 | the paper's supporting information (S1) | CC BY (the article's Creative Commons Attribution licence, which covers its supporting information) |
| `tagu2022cog`, `tagu2022cog2items` | Tagu & Kristjánsson (2022). The selection balance: Contrasting value, proximity and priming in a multitarget foraging task. *Cognition*, 218, 104935. doi:10.1016/j.cognition.2021.104935 | [osf.io/e7cnd](https://osf.io/e7cnd/) | CC BY 4.0 |
| `hughes2024rsos` | Hughes, Statham & Clarke (2024). The effect of target scarcity on visual foraging. *Royal Society Open Science*, 11(12), 240060. doi:10.1098/rsos.240060 | [github.com/Riadsala/foraging_scarcity](https://github.com/Riadsala/foraging_scarcity) | CC BY 4.0 |

## `stim`: one row per item per trial (all items displayed)

| Column | Type | Meaning |
|---|---|---|
| `person` | numeric | Participant, numbered 1..L with no gaps |
| `condition` | factor | Experimental condition; every level used (K = number of levels) |
| `trial` | numeric | Trial, numbered 1..T with no gaps, unique across people and conditions |
| `id` | numeric | Item within its trial, 1..nt |
| `x`, `y` | numeric | Item position (see coordinates below) |
| `v` | factor | Optional item class (cat models); absent = spatial model |
| `trial_p` | numeric | Optional: the trial's original within-person/condition number (kept by `fix_person_and_trial()`) |

Rules:
- Rows are sorted by `trial`, and within each trial `id` is exactly 1..nt in
  row order. `found$id` is used as a position index into its trial's rows.
- Each trial belongs to exactly one `person` and one `condition`.

## `found`: one row per selection, in selection order

| Column | Type | Meaning |
|---|---|---|
| `trial` | numeric | Same numbering as `stim$trial` |
| `id` | numeric | Which item was selected (matches `stim$id` within that trial) |
| `n` | numeric | Selection order within the trial: 1 = first selection, 2 = second, ... |

(Importers usually carry the other `stim` columns through as well; `prep_stan_list()` only reads the three above.)

Rules:
- Sorted by `trial`, then `n`; `n` is exactly 1..f within each trial.
- Every trial has at least one selection. A trial with none shifts `f`
  against every later trial.
- No item is selected twice in a trial; `found` can have fewer rows than
  `stim` (e.g. distractors, or unfinished trials).

## Coordinates

- (0, 0) is the bottom-left of the display. Flip y if the raw data are in
  screen/pixel coordinates (top-left origin).
- x is scaled to [0, 1]; y is scaled by the same factor, so it runs 0..aspect
  ratio. `prep_stan_list()` then multiplies both by `xy_scaling` (default 10).

## Dataset notes

- `kristjansson2014plos`: rotated 180 degrees on import, so that first
  selections start top-left as in the other datasets. The raw
  coordinates' orientation hasn't been confirmed, so this is an
  assumption (`import_data()` says so on every load).
- `hughes2024rsos`: `hughes2024rsos_stim.csv` and `_found.csv` combine
  the per-participant files (github.com/Riadsala/foraging_scarcity).
  `scarcity = TRUE` (the default)
  gives K = 6 conditions (task x scarcity), `FALSE` K = 2 (feature,
  conjunction).
- `tagu2022cog` has 3 item classes; `tagu2022cog2items` has 2 (the
  high-value targets, and the rest).

## Adding a dataset

1. Put the files in `data/<name>/` and add `import_<name>()` to
   `functions/import_data.R`, registering it in `import_data()`'s `switch()`.
2. Use `fix_person_and_trial()` to renumber person and trial, and follow the
   rules above.
3. Add it to `tests/test_dataset_import.R` and run
   `Rscript tests/test_dataset_import.R`.
4. Add its source to the table under Sources above.
