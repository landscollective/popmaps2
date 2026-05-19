# Geographic vs. Suitability-Weighted Surfaces

This note records the intended scientific meaning of `surface = "G"` and
`surface = "C"` before the `surface = "C"` engine is modernized. The goal is to
preserve the flexibility of POPMAPS 1.03 while making the model-choice workflow
more explicit, faster, and easier to validate.

## Scientific Purpose

The original purpose of offering both `G` and `C` was to let users choose the
distance surface that best represents the biology of a focal species.

- `surface = "G"` represents isolation by distance. It asks whether ancestry
  patterns are best interpolated over geographic space alone.
- `surface = "C"` represents a biologically informed landscape surface. In the
  original examples, raster values are MaxEnt species distribution model
  logistic values: probabilities or suitability scores describing whether the
  focal species can occupy each cell. This asks whether ancestry patterns are
  better explained by movement or gene flow through suitable habitat than by
  straight-line distance alone.

The `C` option should therefore not be treated as only an implementation detail
or performance problem. It is central to testing whether genetic patterns are
shaped by isolation by distance, isolation by environment, habitat-mediated gene
flow, or similar landscape processes.

## Current POPMAPS 1.03 Behavior

The legacy implementation builds `surface = "C"` distances with:

```r
gdistance::transition(input_raster, transitionFunction = mean, directions = 8)
gdistance::geoCorrection(..., type = "c", scl = TRUE)
gdistance::costDistance(...)
```

Because the transition layer is built from the mean of neighboring raster cell
values, higher SDM logistic values behave like higher conductance or easier
movement. Lower values increase the effective least-cost distance. This matches
the original Hilaria jamesii example, where higher habitat suitability was
hypothesized to be positively related to gene flow.

Important legacy details to preserve or make explicit:

- Under `surface = "G"`, raster values do not influence distances or ancestry
  weights. The raster provides the grid, extent, resolution, and missing-value
  mask.
- Under `surface = "C"`, raster values influence the distance `D(i)` used by the
  ancestry interpolation equation.
- `threshold` skips ancestry coefficient estimation for focal cells below the
  suitability cutoff, but the legacy code does not turn sub-threshold cells into
  hard barriers for least-cost paths. If `popmaps2` adds barrier behavior later,
  it should be controlled by a separate explicit option.
- Missing raster values remain non-estimable cells and should not be silently
  converted to traversable habitat.

## Shared Ancestry Equation

Both `G` and `C` use the same ancestry interpolation equation. The surface only
changes the distance term. For each focal cell `c`, population or ancestry-axis
coefficient `k`, and selected empirical site `i`:

```text
q_k(c) = sum(q_k(i) * exp(popmod * D(i))) / num_tested
```

For `surface = "G"`, `D(i)` is geographic distance between the cell and the
empirical site. For `surface = "C"`, `D(i)` is suitability-weighted least-cost
distance across the raster surface.

The ancestry probability score is then calculated from the dominant estimated
coefficient after rescaling relative dominance above 0.5 to the 0-1 probability
range. This means surface choice affects both estimated ancestry coefficients
and downstream confidence in the hard-boundary assignment.

## Model-Selection Goal

The goal is not to force all species into one default parameter combination or
one default surface. The goal is to ask which biologically plausible surface and
parameter combination predicts withheld empirical ancestry estimates best while
remaining honest about uncertainty.

For each species, `popmaps2` should be able to compare:

- a geographic-distance model (`G`);
- an SDM suitability/conductance model (`C`);
- eventually, multiple candidate `C` surfaces or transformations when users have
  competing hypotheses about habitat, environment, dispersal, or resistance.

The comparison should use matched validation folds and the same primary metric.
Distance-scale parameters should be suggested from each candidate surface unless
the user deliberately supplies a shared grid. Geographic kilometers and
least-cost distances do not have the same units, so forcing both through one
`popmod` or `empirical_pt_dist` grid can make surface comparison look more
precise than it really is. Repeated spatial-block validation is especially
important because it asks whether the selected surface predicts across space,
not only whether it interpolates well at nearby leave-one-out sites.

## Interpreting G vs. C

Suggested decision language:

- **C supported**: the suitability-weighted surface has meaningfully lower
  validation error than `G`, the improvement is stable across repeated folds,
  and near-best support is not so broad that many incompatible models are
  effectively tied.
- **G supported**: geographic distance is as good as or better than `C`, implying
  the supplied landscape surface does not improve ancestry prediction for the
  empirical data at hand.
- **Indistinguishable**: differences between `G` and `C` are small relative to
  repeat-level uncertainty.
- **Unstable or insufficient data**: the selected surface changes across folds,
  error is high, or support is broad enough that the analysis should be treated
  as exploratory.

This is intentionally conservative. Because ancestry probability surfaces can
guide management, a complex surface should only be preferred when it provides a
repeatable predictive gain.

## Modernization Requirements

The modern `surface = "C"` implementation should:

- preserve legacy suitability-as-conductance behavior as the default;
- precompute least-cost distances from empirical sites to raster cells once per
  surface, then reuse those distances across tuning parameters and validation
  folds;
- precompute empirical-site least-cost distances for enforcing
  `empirical_pt_dist`;
- validate that conductance values are finite and non-negative, with a clear
  handling rule for zeros and missing cells;
- compare a small modern least-cost distance matrix against the legacy
  `gdistance` result before replacing the old engine;
- keep `threshold` as an output/prediction mask unless the user explicitly asks
  for thresholded cells to become movement barriers;
- record the surface type and any transformation in tuning outputs so reports
  are scientifically interpretable.

Future flexibility can include optional transformations, but the default should
remain faithful to the original SDM logistic-value workflow:

- suitability as conductance;
- resistance as cost;
- transformed suitability, such as `1 - suitability` or inverse suitability,
  only when the user requests that biological hypothesis explicitly.

## Accepted Surface Input Scope

`popmaps2` should be good at using outputs from other tools, not at replacing
those tools. Users should be able to supply candidate surfaces from:

- MaxEnt, ENMeval, biomod2, or other SDM workflows;
- Circuitscape, Omniscape, ResistanceGA, or expert-built resistance and
  conductance rasters;
- EEMS, FEEMS, or reems-derived effective migration surfaces, after users export
  them to a raster or other supported distance representation;
- categorical land-cover rasters, once the user supplies an explicit
  reclassification table;
- precomputed distance matrices, in a later advanced interface.

The package should not run SDMs, EEMS, FEEMS, Circuitscape, or ResistanceGA as
part of the core workflow. Those tools require their own modeling assumptions,
inputs, dependencies, and validation. `popmaps2` should instead document how to
bring their outputs into a POPMAPS ancestry interpolation and validation
workflow.

For `surface = "C"`, candidate rasters should declare what their values mean:

- `surface_values = "suitability"` for SDM logistic or habitat suitability
  values. These are conductance-like: higher values imply easier movement or
  stronger habitat-mediated gene flow.
- `surface_values = "conductance"` for surfaces already expressed as movement,
  migration, or connectivity. EEMS/FEEMS-derived rasters usually belong here,
  with a warning about circular validation when the same genetic data generated
  the surface.
- `surface_values = "resistance"` for surfaces where higher values mean harder
  movement. These can be converted internally to conductance with an inverse
  transform before least-cost distances are calculated.

Prediction masks and movement barriers should remain separate inputs. A mask
answers "where should POPMAPS estimate ancestry?" A barrier answers "where
should movement be disallowed or strongly constrained?" Keeping them separate
avoids silently turning a suitability threshold into a biological wall.

## Validation Plan

1. Add tiny test rasters where expected `C` distances can be verified against
   `gdistance`. `(initial helper tests added)`
2. Add a modern internal least-cost distance helper that returns cell-site and
   site-site distance matrices. `(initial helper added)`
3. Extend `popmaps(surface = "C")` to use the modern helper while preserving
   POPMAPS 1.03 behavior.
4. Extend `tune_popmaps(surface = "C")` so tuning can compare `G` and `C` with
   matched folds and repeated spatial blocks. `(initial support added)`
5. Add a candidate-surface comparison wrapper for user-supplied surfaces.
   `(surface-specific grid support added)`
6. Run the empirical examples with `G` and each available SDM-based `C` surface.
   `(initial local report added)`
7. Summarize whether each species supports `G`, `C`, indistinguishable models,
   or unstable support. `(initial local report added)`

## Open Design Questions

- Should `threshold` ever act as a least-cost barrier, or should barrier masks be
  a separate argument such as `barrier_threshold`?
- Should the public API name the default `C` behavior as `conductance` even
  though the legacy argument is `surface = "C"`?
- Which candidate surfaces should be included in public examples, and which
  should remain private empirical-validation data?
- How large does a predictive improvement need to be, relative to repeated-fold
  uncertainty, before the package labels `C` as supported?
