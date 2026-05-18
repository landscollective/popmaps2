# Empirical Tuning Notes

These notes summarize local validation runs across the empirical example data
kept outside the package repository in `../popmaps_test_data`. The raw example
files and generated tuning outputs are not committed to the package.

## Current Interpretation

The empirical examples support treating POPMAPS tuning as a species-specific
model-selection step rather than as a search for universal defaults.

In the full leave-one-site-out plus single spatial-block run from
`report-20260517-184427`, all six species had unique best parameter sets under
both validation designs. Leave-one-site-out validation was usually sharper and
lower-error. Spatial-block validation was stricter and exposed cases where the
data did not strongly identify a single parameter combination.

In the repeated spatial-block run from `report-20260517-191606`, five rotated
spatial partitions were evaluated for each species. The best parameter set was
again unique for every species. Repeat-level uncertainty showed that some
species have much more block-layout sensitivity than others.

## Repeated Spatial-Block Summary

| species | best RMSE | RMSE repeat SD | near-best fraction | support | signal | best num_sites | best num_tested | half-distance km | 10%-distance km | empirical point distance km |
| --- | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| ASLO | 0.1874 | 0.0356 | 0.1087 | moderate | moderate | 5 | 2 | 55.6 | 184.7 | 0.0 |
| CLLU | 0.2487 | 0.0354 | 0.0702 | sharp | moderate | 5 | 2 | 53.1 | 176.4 | 0.0 |
| CLSE | 0.1972 | 0.0152 | 0.0345 | sharp | moderate | 16 | 3 | 83.3 | 276.8 | 200.4 |
| HEMU | 0.2377 | 0.1025 | 0.0901 | sharp | moderate | 5 | 2 | 76.0 | 252.4 | 118.0 |
| MACA | 0.2186 | 0.0270 | 0.3161 | broad | weak | 21 | 2 | 133.7 | 444.0 | 202.2 |
| SPPA | 0.3582 | 0.1491 | 0.3274 | broad | weak | 5 | 2 | 62.4 | 207.2 | 0.0 |

## Takeaways

- Parameter tuning appears useful: the examples do not collapse onto one shared
  parameter combination.
- Single spatial-block validation can overstate uncertainty for some species.
  HEMU looked broad/weak in the single-block report, but repeated spatial blocks
  recovered sharper parameter support while still showing high repeat-level
  RMSE variability.
- MACA and SPPA remain broad/weak under repeated spatial blocks. For these
  species, the current geographic-distance model should be treated cautiously:
  multiple parameter combinations perform similarly, and the ancestry surface
  should not be interpreted as finely resolved.
- CLSE favors a broader site pool and spatial rarefaction under spatial-block
  validation, which may indicate that predictions into unsampled regions benefit
  from more spatially distributed empirical information.
- Repeated spatial-block validation should be the default empirical check before
  comparing geographic distance (`surface = "G"`) with cost/resistance surfaces
  (`surface = "C"`).

## Reproducible Commands

```sh
POPMAPS_EXAMPLE_VALIDATION=spatial_block \
POPMAPS_EXAMPLE_BLOCK_REPEATS=5 \
POPMAPS_EXAMPLE_BLOCK_SEED=42 \
Rscript tools/validate-example-tuning.R \
  ../popmaps_test_data \
  ../popmaps_test_data/tuning_outputs

Rscript tools/summarize-example-tuning.R \
  ../popmaps_test_data/tuning_outputs
```

