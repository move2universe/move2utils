# move2utils 0.4.5 — cleaned-track return no longer carries flag columns

`mt_clean_track(x)` with the default `remove = TRUE` now returns **only
the caller's original columns** (minus the flagged rows). The cascade's
annotation columns (`is_outlier`, `flagged_by_*`, `loglr_*`,
`combined_evidence`, `block_id`, `error_class`, `flag_iteration`) are no
longer attached to the removed-rows return.

**Why.** On a `remove = TRUE` return the object is *shrunk*, so those
columns no longer align to the original row indices. Indexing them by a
pre-removal position (e.g. a row number computed on the full input)
silently reads a *different, shifted fix* — which reads as "the cascade
missed an outlier it actually removed". Dropping the columns makes that
misuse impossible instead of merely discouraged.

**No behaviour change to detection.** Flags, consensus, and which rows
are removed are all unchanged. To inspect flags, call with
`remove = FALSE` (as `mt_diagnose_flags()` / `mt_diagnose_clean_track()`
already require).

## Post-release documentation, folded into 0.4.5 (2026-08-17 → 2026-09-18)

The tag `v0.4.5` was moved on 2026-09-18 to include the documentation line written
after the July release, so that `main`, the tag, the public mirror and the installed
package all describe one state. **No detection path changed:** flags, consensus and
removed rows are identical to the July build. Two user-visible *message* strings
changed, one function returns two extra informational fields, and the manual pages
were regenerated. In order:

- **2026-08-17.** The refuted "empirically validated default / best F1" claim removed
  from `?mt_clean_track`, `man/` and this file; the replacement states that there is
  *no evidence* `evidence_corroborated` detects better than `class_aware`, not that the
  two are equivalent. `HEURISTICS.md` gains the default rule's constants
  (`evidence_threshold`, `evidence_C`, `solo_cols`, the saturation multiple).
  `?mt_clean_track` gains `@section Runtime`. `.gap_aware_autodiff()` now reports
  whether gap-aware scaling is in force (fields `gap_dependent`, `gap_reason`; a note in
  the narration when it is not) instead of falling back to a constant scale silently.
- **2026-08-18.** The 55 m/s ceiling's provenance corrected — below.
- **2026-09-08.** `inst/CITATION`, `README.md` and four `@references` blocks updated:
  both manuscripts are submitted to *Methods in Ecology and Evolution* and on bioRxiv
  (Application Note 10.64898/2026.07.07.736908; outlier paper 10.64898/2026.07.11.737894).
  They had said "in preparation" for fifteen days after the outlier preprint appeared.

### The 55 m/s sanity ceiling is no longer attributed to a Hirt et al. (2017) confidence bound

The auto-cap's biological-sanity ceiling was documented in four places — the roxygen
of `mt_flag_speed_cap()` and `mt_suggest_speed_cap()`, both functions' runtime
messages, and `HEURISTICS.md` — as sitting "above the Hirt et al. (2017) 95% upper CI
of the maximum biological speed across all masses and modes (~52.6 m/s at 1.35 kg)".

**That attribution was wrong in two ways.**

- **The interval is not Hirt's.** Hirt et al. publish fitted parameters with marginal
  standard errors, not an interval on maximum speed. The ~52.6 m/s figure is this
  package's own delta-method propagation of those SEs, reproducible as the upper bound
  of `v_phys_estimate(mass = 1.35, mode = "flying")` (52.56 m/s). Citing it to the
  article attributed to it an interval it never contained.
- **"Across all masses and modes" was false.** That upper CI **exceeds 55 m/s for
  flying between 1.69 kg and 15.65 kg, peaking at 60.88 m/s at 4.89 kg**, because the
  interval widens with mass while the central prediction falls away from its peak near
  1.38 kg. The bound was cleared only at the single mass quoted.

The interval is in any case too crude to carry the claim: it uses marginal SEs with no
parameter covariance (a known limitation), giving [1.23, 60.88] m/s at 4.89 kg and
[0.00, 58.52] at 10 kg.

**What the ceiling is now documented as:** an impossibility bound on speed *sustained*
over the interval separating two fixes. Nothing living sustains 55 m/s over such an
interval, so a step implying one is positional error or a tag no longer on the animal
rather than movement. This needs no interval and no allometry, and it is what the cap
actually tests — a stooping falcon exceeding 55 m/s for seconds is not a counterexample.

The runtime message now reads `(universal sustained-speed bound, not species-specific)`
where it previously named the Hirt upper CI. **The ceiling value is unchanged at 55 m/s
and remains warning-only**; it never alters the cap.

In `HEURISTICS.md` the row moves from the **DERIVED** table to **HEURISTIC**: a round
number above a stated bound is a documented choice, not a published constant. The
per-species allometric `v_max` row is unaffected and stays DERIVED — `v_phys_estimate()`
traces cleanly to the parameters of Hirt et al.'s Supplementary Table 4, which are
printed in `?v_phys_estimate`.

Found by the companion outlier-paper project, verified here against source and arithmetic.

### Manual pages regenerated

`man/` had drifted from the roxygen it is generated from: the corrected scope caveat on
the detour leg-gate argument in `?mt_clean_track` was present in `R/` but missing from
the shipped manual page. Five pages are regenerated, so `R/` and `man/` agree again.

# move2utils 0.4.4 — block-expansion monotonicity fix

`mt_clean_track()` block expansion is now consistent across the
auto-cap and physiological-cap paths: **supplying a physiological cap
(`v_max`, or `mass` + `mode`) no longer reduces recovery of a coherent
boundary block.**

**Bug.** A sustained, boundary-anchored spoof block (a run of
internally-consistent erroneous fixes pinned off-route, so only the
seam steps look fast) was recovered in full on the auto path but
largely missed when a physiological cap was supplied. Two causes
compounded: (1) supplying a cap disabled block expansion outright, and
(2) the iterative speed-peel removed the block's seam fixes, which
widened the time gap across the seam so the implied across-seam speed
fell below the cap — the block then re-joined the main trajectory in
the connectivity graph and no block remained to flag.

**Fix.**

- The component partition underlying block expansion
  (`.compute_component_partition`) now uses original, un-diluted
  timing: across a gap created by removed fixes, the connectivity lag
  is capped at the track's median sampling interval, so a large
  displacement across a removed span stays severable. Genuine
  single-step gaps (including real missing-data gaps) are unaffected,
  and isolated spikes still correctly re-join (their kept neighbours
  are spatially close).
- Block expansion is no longer skipped when a physiological cap was
  supplied. Recovery is now monotone in user information: a supplied
  cap — the connectivity ceiling block expansion needs — can only help.
- A latent off-by-one in the partition (a phantom unit-size component)
  is corrected.

Behaviour is unchanged on clean tracks, on mid-track blocks (recovered
by the per-fix layer, with the gate correctly declining), and on the
golden-eagle K02 case (byte-identical). Regression-tested with the new
`inst/extdata/make_boundary_spoof_demo.R` fixture.


# move2utils 0.4.3 — mt_thin_distance overhaul

`mt_thin_distance()` was reworked after a check found its behaviour
matched neither its own documentation nor the legacy
`move::thinDistanceAlongTrack()` it claimed to replicate.

**Behaviour change.**

- `mt_thin_distance()` now offers two thinning rules via a new `method`
  argument. The default `method = "step"` is a greedy minimum-spacing
  walk: every retained pair is at least `distance` of along-track travel
  apart — the intuitive "one fix per `distance` of travel". The previous
  implementation used a bucket-from-origin rule that could retain fixes
  much closer than `distance`; its output is not reproduced.
- New `method = "interval"` is the distance-space twin of
  `mt_thin_time()`: it selects the longest subset whose successive
  retained pairs lie within `[distance - tolerance, distance + tolerance]`
  of travel, using the same dynamic-programming sweep. The new
  `tolerance` and `criterion` arguments serve this mode (they warn / are
  ignored under `method = "step"`).

**Documentation.**

- The `mt_thin_distance()` help and the "Interpolating and thinning"
  vignette were corrected to describe the actual behaviour; the old text
  described a since-last rule the code never implemented.

# move2utils 0.4.2 — conceptual-integrity fixes

A batch of cross-cutting consistency, correctness, and robustness fixes
surfaced by a package-wide conceptual-integrity audit (auditing the
cross-cutting *contracts* rather than individual functions).

**Breaking change.**

- `mt_dbgb_variance()` argument order changed from
  `(object, location_error, margin, window_size)` to
  `(object, location_error, window_size, margin)` so it matches
  `mt_dbbmm_variance()`. Calls using named arguments are unaffected;
  positional calls that passed `margin` before `window_size` must be
  updated.

**Behaviour changes.**

- The outlier detectors now **reject empty or non-finite ("NA") geometry**
  with a classed error (`move2utils_input_has_empty_geometries` /
  `move2utils_input_nonfinite_coords`) instead of silently dropping those
  fixes. There is no outlier to identify in a fix with no location, so it
  must be removed upstream (e.g. `x <- x[!sf::st_is_empty(x), ]` or
  `mt_filter_gps_quality(x)`). Non-finite *timestamps* are still excluded
  with a note rather than rejected.

**Fixes and consistency.**

- `ud_outer_probability()` now reprojects `sf`/`move2` query points to the
  UD's CRS automatically; previously a CRS mismatch returned silently-`NA`
  or wrong probabilities. (Same projection-independence class as v0.4.1.)
- `mt_thin_distance()` now measures distances in metres
  (`mt_distance(units = "m")`); previously the CRS's linear unit could leak
  into the threshold comparison on non-metre projected CRSs.
- `mt_clean_track()` now retains the per-detector `loglr_*` columns on its
  output, so `mt_flag_consensus()` can re-score a cleaned track under the
  default `evidence_corroborated` (or `weighted_evidence`) mode instead of
  erroring; and the standalone `mt_flag_outliers()`,
  `mt_flag_outliers_bridge()`, and `mt_flag_speed_cap()` now each emit a
  `flagged_by_<detector>` column (like `mt_flag_outliers_detour()` already
  did), so independently-run primitives are votable by `mt_flag_consensus()`.
- `mt_flag_consensus(pool_by = )` is now resolved through the same helper
  the detector primitives use: it looks the column up in
  `mt_track_data(x)`, accepts the nested `c(outer, inner)` form (calibrating
  on the outer, distribution level), and errors on a missing column instead
  of silently grouping per track.
- `mt_flag_outliers_dbgb()` / `mt_flag_outliers_dbbmm()` now accept
  `location_error` in the same forms as the other primitives (`NULL` /
  scalar / per-fix vector / column name / `"auto"`); previously a column
  name or `"auto"` reached the variance denominator unresolved and errored,
  and a per-fix vector could misalign after the optional pre-peel. The
  default changed from `0` to `NULL`.
- Duplicate- and unsorted-timestamp errors now carry the shared
  `move2utils_input_duplicate_timestamps` /
  `move2utils_input_unsorted_timestamps` condition classes (in addition to
  their previous function-specific class), so `tryCatch()` handlers keyed on
  the shared class work uniformly across functions.

# move2utils 0.4.1 — projection-independent geometry

**Coordinate reference systems no longer change the result.** Every
geometric computation is stratified onto a consistent local metric
space, so a track's analysis no longer depends on the projection it
happens to arrive in. This is a robustness fix with no API changes;
existing projected workflows are unaffected.

- **Outlier detection is projection-invariant.** `mt_clean_track()` and
  the standalone detectors (`mt_flag_outliers()`,
  `mt_flag_outliers_bridge()`, `mt_sequential_outliers()`) now compute
  all geometry in a per-track local azimuthal-equidistant (AEQD)
  projection derived from the data, regardless of the input CRS. The
  same track yields the same flags whether supplied in
  longitude/latitude, UTM, or any other system, and the flags and
  log-likelihood-ratios are returned on the original object in its
  original CRS. (Previously a caller-supplied projection could leak into
  the result.) `mt_sequential_outliers()` also gains projected-CRS
  support — it no longer relies on the longitude/latitude-only
  `mt_turnangle()`.

- **UD and variance functions accept any CRS and return it in kind.**
  `mt_dbbmm_variance()`, `mt_dbgb_variance()`, `mt_dbbmm_ud()` and
  `mt_dbgb_ud()` no longer require projected input. Longitude/latitude is
  auto-projected to a local AEQD for the bridge math and the UD is
  reprojected back to the caller's CRS and renormalised; a projected
  metric CRS is used as supplied (longitude/latitude in →
  longitude/latitude out, UTM in → UTM out). Passing an environmental
  `SpatRaster` as the output grid returns the UD on that exact grid,
  ready to stack against environmental layers. (Previously these
  functions errored on longitude/latitude input.)

# move2utils 0.4.0 — evidence-corroborated cascade default, classed conditions, projection-agnostic `mt_corridor`

**New default outlier-decision rule (`evidence_corroborated`).** The
`mt_clean_track()` cascade — and the standalone `mt_flag_consensus()` —
now turn the four detectors' flags into outlier calls with an
evidence-accumulation rule instead of the previous Boolean class-rule
(`class_aware`). Each detector emits a signed log-likelihood-ratio
(surprisal minus its own data-driven flag boundary, exposed as the
`loglr_*` columns); these are made commensurable
(`C·tanh(loglr / MAD)`), summed into a single `combined_evidence` score,
and a fix is flagged when that evidence is positive **and** either at
least two detectors corroborate, or the high-specificity detour detector
is overwhelming (saturated). It is data-driven, principled, and
unsupervised — no user input or labels. On current evidence it is not more
accurate than the previous `class_aware` rule — the two differ by less than
the realisation noise, which is an absence of evidence for a difference
rather than a demonstration of equivalence — so the default changed on
principled grounds rather than on measured skill. Coherent
spoofing/GPS-jamming blocks are still removed at full recall; which layer
removes them depends on where the block sits — the per-fix detectors and the
consensus rule for a mid-track block, graph-based block expansion for a
boundary-anchored one. The previous behaviour is available with
`consensus = "class_aware"`. The smooth `combined_evidence` column
doubles as an optional one-lever sensitivity control. See
`?mt_clean_track`, `?mt_flag_consensus`, and
`DESIGN_evidence_accumulation.md`.

The remaining changes below also ship in this release:

**Classed errors and warnings throughout.** Every `stop()` and
`warning()` in `R/` now uses `rlang::abort()` / `rlang::warn()`
with a stable `class = "move2utils_<context>_<reason>"`.  Downstream
code can now key off the class in `tryCatch()` /
`withCallingHandlers()` rather than regex-matching the message text.
Common shared classes:

- `move2utils_input_not_move2` — the input isn't a `move2` object.
- `move2utils_input_longlat_required_projected` — a function that
  needs a projected CRS got longlat input.
- `move2utils_input_nonfinite_coords`,
  `move2utils_input_unsorted_timestamps`,
  `move2utils_input_duplicate_timestamps`,
  `move2utils_input_has_empty_geometries` — track-hygiene
  preconditions raised by `.validate_move2()`.
- `move2utils_internal_*` — invariant violations / "this shouldn't
  happen" paths.

Function-specific classes follow `move2utils_<function>_<reason>`
(e.g. `move2utils_mt_clean_track_bad_v_max`,
`move2utils_mt_flag_outliers_bad_prob_type`).  Message text is
unchanged, so existing `expect_error(..., regexp = "...")` tests
continue to match.  When errors fire and what they mean is
unchanged — only how they are caught.

**`mt_corridor()` rewrite** — projection-agnostic (handles both
longlat and projected input), track-aware (via `mt_segments()`, so
multi-track inputs never produce cross-track segments), handles
empty geometries with `NA`-padding, and uses `mt_aeqd_crs()`
internally for the longlat-to-metric transform.  Thresholds are
user-supplied numbers (the legacy `speed_quantile = 0.75` /
`circ_quantile = 0.25` quantile fallback is preserved but now emits
a warning naming the resolved numeric).  `corridor_azimuth` is now
in degrees (matching the docstring; the previous implementation
stored radians).  Drops the `circular` and `geosphere` package
dependencies in this code path.  See `?mt_corridor` and the
rewritten `vignette("corridor")` for the new idioms.

**`mt_clean_track(state = ...)` no longer crashes on degenerate
state subsets.** When cleaning is dispatched by behavioural state,
each state value is cleaned on its own subset.  Small or stationary
subsets exposed two crashes in the probability detector, both
surfaced (not caused) by the 0.3.3 propagation fix that made the
gap-based threshold the effective default prob path:

- A subset with too few valid probabilities falls back from
  gap/entropy/significance to the percentile method.  The fallback
  carried the originating method's threshold — which may be `NULL`
  (deferred to a leaf formal) or a gap-spacing multiple, neither a
  percentile fraction — into `100 - threshold * 100`, producing a
  zero-length flag vector and then `replacement has length zero`.
  The percentile branch now coerces any non-`(0, 1]` value to the
  percentile default and guards against an empty probability set.
- An all-stationary subset (repeated positions, every step length
  `0`) drove the internal `(turn, step)` histogram to a terra raster
  with a zero-height `ymin == ymax` extent (`[rast,missing] invalid
  extent`).  `.turn_step_hist()` now returns `NULL` on degenerate
  input (no positive step length or no valid turn angle) and the
  caller treats that as "no kinematic signal" (those fixes are kept).

Single-state cleaning was unaffected (a full track almost never has
fewer than ten valid probabilities).  Regression tests cover both
the internal helpers and the end-to-end multistate path.

# move2utils 0.3.3 — threshold-parameter propagation fix (2026-05-25)

Fixes a layered shadowing bug whereby user-supplied entropy/gap
thresholds were silently dropped before reaching the leaf detectors
(`.entropy_threshold_lower`, `.gap_threshold_lower`).  Pre-0.3.3
roughly twelve sites in the package re-declared the
`(entropy = 0.3, gap = 3)` mapping and forwarded the literal as an
explicit argument to the leaf, shadowing both any namespace-level
override and any user setting passed at the cascade entry point.
`mt_clean_track()` additionally bypassed the user-facing primitive
wrappers entirely (calling internal `.fn_core` siblings directly),
so even the existing standalone-primitive surfaces (`threshold = X`)
could not propagate through the cascade.

Pattern A applied consistently from leaf to user-facing layer: every
wrapper defaults the relevant threshold to `NULL` and forwards only
when non-NULL; the leaf formal is the single source of truth for the
package-wide entropy / gap defaults.

Three new knobs on `mt_clean_track()`:

- `entropy_threshold = NULL` — applies wherever the cascade runs an
  entropy-valley detector (bridge primitive when
  `bridge_threshold_type = "entropy"`, prob primitive when
  `prob_threshold_type = "entropy"`, speed-cap auto path's entropy
  arm).
- `gap_threshold = NULL` — applies wherever the cascade runs a
  broken-stick gap detector (bridge / prob / speed-cap auto path's
  gap fallback / class-aware persistence post-filter).
- `persistence_filter_threshold = NULL` — per-scale gap-threshold
  passed to the class-aware persistence post-filter.

`NULL` (the default everywhere) defers to the leaf's formal default
(0.3 entropy / 3 gap) — byte-identical to pre-0.3.3 behaviour at the
defaults.  Non-NULL scalar overrides.

Three new knobs on `mt_flag_outliers_dbgb()` for the absolute-
residual axis:

- `residual_entropy_threshold = NULL`
- `residual_gap_threshold = NULL`
- `residual_dip_alpha = 0.05`

These control the `.auto_residual_max` gate that runs when
`residual_max = NULL`.

`mt_persistence_score()`'s `threshold` argument changes default from
`3` to `NULL` (defers to the leaf).  Pre-0.3.3 `mt_persistence_score(x)`
ran with hardcoded `3`; post-0.3.3 it ran with leaf-formal `3` — the
same numeric value, but now the leaf is the source of truth.

`mt_flag_outliers_dbgb()`'s `z_gap_threshold` default changes from
`3` to `NULL` (defers to `.z_gap_threshold`'s leaf formal).

Sweep harness (`sweep/2026-05-02-heuristic-sensitivity/`) consumers:
the `assignInNamespace` patches on the two leaf formals now reach
through every cascade path; the `bridge_wrapper` injection in
`01_run_one_combo.R` becomes redundant and is slated for removal in
the next sweep refresh.

Regression tests in `tests/testthat/test-propagation.R` lock the
chain via direct spies on the leaves: each user-facing surface is
exercised with a non-default threshold and the leaf's recorded call-
time value is asserted to match.  Pre-fix these tests fail; post-fix
they pass and pin the propagation chain against future drift.

See `audits/2026-05-25-parameter-propagation/findings.md` for the
full audit (propagation + adjacent code-health / hygiene dimensions).


# move2utils 0.3.2 — multi-track `location_error` slicing fix (2026-05-14)

Fixes the move→move2 porting bug whereby
`mt_dbbmm_variance()` / `mt_dbgb_variance()` rejected a per-fix
`location_error` vector on multi-track input with *"location_error
must be length 1 or equal to the number of locations"* — the
inherited single-track helper received the full object-length vector
unchanged from `.dbbmm_variance_multi` / `.dbgb_variance_multi`.
Same failure propagated into `mt_dbbmm_ud.list()` /
`mt_dbgb_ud.list()`.

`location_error` is now treated as a unified per-fix property of the
move2 object, in metres, with the same five forms accepted by the
bridge primitive:

- `NULL` (no per-fix error contribution),
- numeric scalar (uniform sigma),
- numeric vector of length `nrow(object)` (per-fix sigma; sliced
  per-track internally for multi-track input),
- column name in `object` (column values used as metres),
- the literal `"auto"` (probes `eobs_horizontal_accuracy_estimate`
  then `argos_lc`; HDOP-style columns are deliberately not
  auto-resolved — conversion to metres needs a sensor-specific
  multiplier the package cannot infer).

Two new things:

1. `mt_dbbmm_variance()` and `mt_dbgb_variance()` gain a
   `location_error_na` argument (`"median"` / `"mean"` / `"zero"` /
   `"approx"`, default `"median"`) for imputing NAs in a per-fix
   vector.  Per-track median fill is the default because zero-fill
   fakes certainty: in the variance/UD layer a zero anchor error
   tightens the bridge variance and produces an artificially
   confident UD precisely where the data are least informative.
2. The variance object now carries the resolved per-fix
   `location_error` on `track_data$location_error`, so
   `mt_dbbmm_ud()` / `mt_dbgb_ud()` re-use it automatically when
   called on a variance object or a named list of them.  The
   `location_error` argument on those methods defaults to `NULL`
   (use stored); explicit overrides still work but per-row vectors
   are rejected at the list-dispatch step because the variance
   object carries no row-mapping back to a flat user vector.

The bridge primitive (`mt_flag_outliers_bridge()`) and the
cascade entry-point (`mt_clean_track()`) are unchanged: their
NA-anchor-fallback-to-S_hat semantics are preserved, so cohort A/B
remains byte-identical.

Twelve new regression tests in
`tests/testthat/test-multi-track-location-error.R` lock the
contract: per-fix vector on multi-track now matches per-track fits
byte-identically; column-name and scalar inputs still work;
wrong-length vectors error clearly; NA imputation defaults to
per-track median; UD list dispatch flows stored vectors through.

Full test suite after the change: 1011 PASS / 0 FAIL / 1 SKIP / 0 WARN.

# move2utils 0.3.1 — beta-release hygiene, consensus refactor, vignette pedagogy (2026-05-13)

A beta-release polish pass.  No code regressions; all 914 tests
remain green.  Most changes are user-facing improvements: a cleaner
public API for the consensus rule, friendlier first-run messages,
a vignette rewrite under a "speak as people speak" voice directive,
a nested-pool extension to `pool_by`, and a new post-run
diagnostic helper.

## New: `mt_diagnose_flags()` for auditing a cleaned object

`mt_diagnose_flags()` reads the per-detector flag columns that
`mt_clean_track()` (or any of the four primitives) attaches and
returns:

- An **error_class breakdown** of where flags land in the
  cascade hierarchy.
- A **per-detector fire table** showing which detectors are doing
  the work on this dataset.
- A **co-fire histogram** crossing "how many detectors fired" with
  `is_outlier`.
- A **near-miss table**: for each detector, how many fixes fired on
  it but were not flagged because no other detector corroborated.
  A loud detector with many near-misses is a candidate for either a
  tighter threshold or a looser consensus rule.
- A **consensus-mode comparison**: under each built-in mode of
  `mt_flag_consensus()` (`class_aware`, `strict`, `majority`,
  `speed_trusted`, `any`), how many fixes that mode would flag and
  the delta from the current `is_outlier`.  No per-fix detectors
  are re-run; the math is post-hoc on existing flag columns.
- A **five-panel map** (single AEQD projection): four small panels
  per detector marking confirmed-by-consensus vs near-miss fires,
  plus one large panel showing the consensus decision coloured by
  `error_class`.

The function is intended as the diagnostic-then-tune loop: run the
cleaning, call `mt_diagnose_flags()`, see where the near-misses
cluster, decide whether to tighten a specific detector, loosen the
consensus rule, or call `pre_peel_aux = "primitives"` for a
geometric-spike pre-peel pass.



## Nested `pool_by`: separate fit-source and operating unit

`pool_by` on `mt_clean_track()` and the four primitives
(`mt_flag_outliers`, `mt_flag_outliers_bridge`,
`mt_flag_outliers_detour`, `mt_flag_speed_cap`) now accepts a
length-2 character vector `c(outer, inner)` in addition to the
existing length-1 form.

- `outer` names the column whose value defines the **fit set**:
  threshold-fitting primitives draw their reference distribution
  from the union of events whose tracks share this column's value.
- `inner` names the column whose value defines the **operating
  unit**: the post-cascade flag union acts within tracks sharing
  this column's value, and the orchestrator's post-cascade sweep
  iterates over these groups.

Length 2 requires strict nesting: every distinct `inner` value
must map to exactly one `outer` value. Inputs that violate this
error with the offending value named.

**Length \eqn{> 2} is rejected** with a deliberately verbose
message. `pool_by` has exactly two semantic roles; a deeper
hierarchy (e.g. species/population/individual/tag) would only earn
its keep under hierarchical / partial-pooling threshold
estimation, which the cascade does not perform. Users with many
nested levels pick the *pair* of columns that captures their
trust claim and operating unit.

Length-1 input is preserved byte-identically (treated as
`c(col, col)`). Validation error messages were rewritten to name
the 2-level cap and the rationale.

New vignette `heterogeneous_error_regimes` walks through the
companion case the buffalo dataset surfaced: mixed-sensor tracks
(GPS + Sigfox interleaved on one device) require **pre-splitting
by sensor** before any pool, since the cascade's per-call
homogeneous-error-regime contract holds. The vignette shows the
`dplyr`-based pre-split + re-merge pattern and how nested `pool_by`
combines with it on multi-individual datasets.

## Externalised consensus rule: `mt_flag_consensus()` + `consensus =`

The per-fix outlier decision (how the four detectors' flags combine
into `is_outlier`) is now a separate concern from the cascade
orchestration, exposed both as a parameter on `mt_clean_track()`
and as a stand-alone function.

**New exported function** `mt_flag_consensus(x, mode = ...)` takes a
move2 object already carrying per-detector flag columns and returns
it with `is_outlier` set according to the chosen consensus mode.
This lets users running their own cascades (or running the four
primitives standalone) apply the same consensus rules to their own
flag columns. Six modes:

- `"class_aware"` (default) — empirically validated one-fits-all
  default. A fix flags if any of the class rules fires (consensus
  ≥ 3 of 4, geometric_spike, state_anomaly, kinematic_confluence).
- `"strict"` — `(bridge | detour) & (prob | speed)`; highest
  precision.
- `"majority"` — flag iff ≥ 2 of 4 detectors agree.
- `"speed_trusted"` — `speed | ((bridge | detour) & prob)`.
- `"any"` — union of all four detectors; maximum recall.
- `"custom"` — user-supplied function of four logical vectors.

`mt_clean_track(consensus = ...)` defaults to `"class_aware"` and
delegates to the same internal core.  The `class_aware` taxonomy
used by `error_class` is independent of the chosen `consensus`
mode — even under e.g. `"majority"` voting, flagged fixes still get
a class label describing the agreement structure.

**Retired** as part of the same refactor: the previous
`flagging_strategy` and `strictness` parameters of
`mt_clean_track()`.  Their combinations are now reachable via
the `consensus =` parameter (`flagging_strategy = "legacy_strict"`
with the four `strictness` modes maps onto
`consensus = "strict"|"majority"|"speed_trusted"|"any"`).
The retired parameters are documented in
`archive/2026-05-13-legacy-strict-retired/README.md` along with how
to reproduce the v0.1.x behaviour if needed.

## First-adopter friction fixes

A first-adopter simulation on a public Movebank Hornbill track
surfaced several rough edges; this release addresses them.

- **Movebank licence acceptance** is now documented next to the
  download example in `vignette("getting_started")`.  Previously,
  the first call to a CC0 study returned an opaque wall of licence
  text and the user had to dig out the `license-md5 = "..."` hash
  themselves.
- **Auto-cap warning is friendlier.**  The previous message
  lectured about K02 spoofs, "block-shaped contamination", and
  "homing-pigeon racing flight + perch" — intimidating jargon for
  first-time users.  The new message leads with "this works well
  for most cases" and points at `?v_phys_estimate` and
  `?mt_clean_track` for detail.
- **Zero-flag console hint.**  When `mt_clean_track()` returns no
  outliers, the console now adds a one-line follow-up:
  "This is normal for clean tracks; for tracks with multiple
  behavioural states the cascade can come up empty until you
  switch to state-conditional cleaning. Run
  `mt_diagnose_clean_track()` to confirm."  Closes a "did it
  work?" interpretation gap.
- **Per-segment dispatch warnings are aggregated.**  State-
  conditional cleaning on a track with many short segments
  used to emit 50+ "Too few locations" warnings.  Now a single
  summary: "State-dispatch summary: N of M segments were too
  short for the cascade and left unflagged (X with <3 fixes, Y
  with 3-9 fixes)."
- **Auto-cap warning no longer repeats in recursive calls.**
  Previously the message fired once per recursive
  `mt_clean_track()` call (per state segment or per individual);
  now it fires at most once at the outer level.
- **Internal function reference removed from a vignette.**
  `vignette("state_conditional")` used to suggest
  `move2utils:::.step_lengths_fast(x)` in user-facing example
  code; replaced with the exported `mt_speed(x, units = "m/s")`.

## Vignette voice rewrite

Six vignettes rewritten under a pedagogic / conversational /
verbose-over-tight voice directive: open with the user's
situation, lead with intuition before mechanism, demote math
notation to "for the curious" sub-sections, write as people
speak:

- `outlier_bridge.Rmd`: new opening + "How it works — intuition
  first" lead, formal math demoted to a sub-section.
- `persistence_score.Rmd`: reordered so the use case ("how
  confident should I be in this flag?") precedes the formal
  mechanism.
- `dbgb_ud.Rmd`: open with the use case (migration leg, daily
  commute) instead of "isotropic Brownian diffusion".
- `bursted_uds.Rmd`: retitled "Producing one UD per behavioural
  state (day/night, season, etc.)".  Legacy `MoveBurst` note now
  parenthetical, not the lead.
- `gap_aware_ud.Rmd`: open with the use case (long gap → smeared
  UD) and a clear two-option overview.
- `getting_started.Rmd`: the ~90-line Movebank download section
  has been moved to an Appendix at the end of the vignette.  The
  main path now has a single cross-reference, so first-time
  readers reach the cleaning workflow immediately.

## Documentation and navigation hygiene

- `_pkgdown.yml` articles regrouped by audience (Outlier cleaning
  → Worked examples → Utilisation distributions → Trajectory
  utilities → Reference notes); explicit "read in this order"
  labels on the cleaning vignettes.  `mt_flag_consensus` added
  to the reference index.
- `README.md`: function-name catch-up (`mt_ud_volume` →
  `ud_volume`, `mt_outer_probability` → `ud_outer_probability`,
  `mt_v_phys_estimate` → `v_phys_estimate`); citation bumped to
  v0.3.1.
- `getting_started.Rmd` "Further reading" annotated with a
  recommended reading order.

## Minimum R version honesty

`DESCRIPTION` now declares `R (>= 4.5.0)` instead of `R (>= 4.1.0)`.
The package has only ever been tested on R 4.5.x (Raven container)
and R 4.6.x (GitLab CI runner, local development).  Claiming 4.1+
without actually testing on it was an overclaim from an earlier
phase.  R 4.5 was released April 2025; users on older R should
upgrade.

## Repository hygiene

- Five top-level dev-history docs (`AUDIT_*.md`,
  `CASCADE_AUDIT_*.md`, `FUNCTION_REVIEW_*.md`,
  `HYGIENE_AUDIT_*.md`, `INVENTORY.md`), the superseded
  `audits/2026-04-29-stratified-movebank/` directory, and the
  contents of `scratch/` moved to
  `archive/2026-05-13-hygiene-archive/`.
- The `flagging_strategy`/`strictness` comparison test moved to
  `archive/2026-05-13-legacy-strict-retired/`.
- Two READMEs document everything in the archive, force-tracked
  in git despite `archive/` being gitignored, so the institutional
  memory survives across machines.

# move2utils 0.3.0 — pool_by, primitive-knob exposure, sequential bug fix, naming hygiene

## `pool_by` cross-track threshold pooling (2026-05-12)

All four outlier-detection primitives (`mt_flag_outliers_bridge`,
`mt_flag_outliers_detour`, `mt_flag_outliers`, `mt_flag_speed_cap`)
and the orchestrator `mt_clean_track` gain a new optional `pool_by`
parameter: a character scalar naming a column in
`mt_track_data(x)`. When supplied, tracks sharing the same value in
that column are pooled for threshold fitting — one threshold is
fit per group from the union of the group's data, then evaluated
per track. Typical scopes: `"individual_id"` (pool deployments of
one animal — the Rhino dual-collar motivating case),
`"study_id"` (population-level), `"species"` (across-species for
sparse-per-individual datasets).

Semantics: pool flags are ADDITIVE only — never un-flag what
per-track threshold caught. `NULL` (default) preserves per-track
behaviour byte-identically. Pool-added flags in
`mt_clean_track` output are tagged `error_class = "pool"`.

**What pool_by helps with** — clarified after empirical scrutiny:
distribution insufficiency for thresholding (not "short tracks" by
fix count). Reasons include too few fixes, clean track with no
extreme tail (Rhino track 2: 1206 fixes, max speed 1.33 m/s — no
anomaly to break from but inherits a useful cap from the sibling
deployment), contamination ratio above the 5% gate, or
state-imbalanced deployment.

`reference =` (single-external-clean-track override on
`mt_flag_outliers`) and `pool_by =` are mutually exclusive.

## Primitive-knob overrides at the cascade level (2026-05-12)

`mt_clean_track` gains five `NULL`-default parameters exposing
methodological choices previously hardcoded inside the cascade:
`bridge_method`, `bridge_threshold_type`, `bridge_iterations`,
`prob_threshold_type`, `detour_threshold_type`. `NULL` preserves the
cascade's empirically-tuned hardcoded value byte-identically.

Two architectural knobs stay hardcoded with explicit rationale
(documented in the new `@section Primitive-knob overrides` of
`?mt_clean_track`): detour `min_leg = 0` (gating is owned by the
conjunction rule) and speed_cap `threshold_type = "auto"` (cascade
decouples the relative-speed conjunction flag from the absolute
`v_max` block-expansion cap).

Empirical demo: `bridge_method = "directional"` on CPF_A yields
F1 = 1.000 (vs 0.958 with the cascade's `"combined"` default) —
eliminates 2 FPs without sacrificing TPs.

## `mt_sequential_outliers` two-threshold calibration fix (2026-05-12)

Pre-fix, the sequential scan calibrated a single threshold from
the full reference formula (`stp * (dsp * dtp)^autodiff_alpha`) and
applied it to both scan-time regimes. On multi-state tracks where
the gap-aware autodifference KDE has very tall peaks on the
bimodal Delta-step distribution (e.g.\ WH17: a Bonelli's eagle's
rest + flight modes produced autodifference densities up to
8 × 10⁷), the single threshold sat orders of magnitude above the
first-step `prob = stp` values that the scan returns when
`prev_anchor` is `NA` — so every first-step was flagged. The
anchor never advanced, `max_skip` cycled, and 99.99% of WH17's
fixes were flagged.

The fix calibrates two thresholds, each as the 0.5%-percentile of
its matching reference distribution: `threshold_partial` for the
partial-formula regime (`prob = stp` only) and `threshold_full`
for the full-formula regime. The scan compares each step's score
to its regime-matched threshold; final flagging is
`prob / threshold-for-this-regime < 1`. User-supplied scalar
`threshold = X` continues to set both thresholds to `X`
(back-compat).

Empirical: WH17 99.99% → 1.42%. CPF synthetic byte-identical (the
single-state autodifference KDE returns ~1 at the bulk, both
thresholds agree numerically). Paper claims (§5.2, §6.4, §7.4)
preserved — only the regime-matching is corrected.

## Naming hygiene: `mt_` prefix is now reserved for move2-touching exports (2026-05-12)

Three exports renamed to drop the `mt_` prefix (the `emd` export
was the existing precedent):

- `mt_v_phys_estimate()` → `v_phys_estimate()` — takes
  `(mass, mode)` scalars, returns a numeric.
- `mt_ud_volume()` → `ud_volume()` — takes a `terra::SpatRaster`
  UD.
- `mt_outer_probability()` → `ud_outer_probability()` — takes a
  `terra::SpatRaster` UD plus query points.

No deprecation shims. Update calls to the new names directly.

## Round-4 micro-audit cycle complete (2026-05-12)

Nine micro-audits closed (A B C D F H I E G); empirical record in
`FUTURE_IMPROVEMENTS.md` and `CASCADE_AUDIT_2026-05-11.md`. Net
code shipped: pool_by (5 commits), Item E sequential fix, Item G
knob exposure, Item F R5 magic-number ratio replacement, Item H
+ Item I closed on grounds (no code change). v0.3.0 ship gate
confirmed bug-free.

Queued for v0.3.x / v0.4 (NOT bugs — methodology improvements):

- **M-12**: per-fix p-value calibration on `mt_sequential_outliers`
  replacing the 0.5%-percentile heuristic with a parametric-null
  FWER/FDR criterion. Estimated 1–2 days + paper §7.4 update.
- **Approach (i) iteration-aware bridge pooling**: passes the
  concatenated pool object to `.bridge_fn_core` with cross-track
  boundary fixes masked. Removes the documented approach-(ii)
  approximation. Defer until empirical evidence warrants.

# move2utils 0.2.0 — class-aware flagging architecture

## `mt_persistence_score()` is now CRS-invariant (2026-05-11)

`mt_persistence_score()` previously dispatched between Haversine
(lon/lat input) and Euclidean (projected input) distance computation
via a coordinate-magnitude heuristic in `.scale_k_geometry`.  Both
methods return metres in principle, but the histogram-based gap
threshold builds bin breaks from `seq(min(log_step), max(log_step),
length.out = n_breaks + 1)`; sub-percent CRS-induced differences in
step length therefore shifted bin boundaries enough to push some
flagged fixes across the gap threshold.  Empirically: on CPF_D's 59
cascade flags, 14 fixes (24%) had different `persistence_count`
between the lon/lat and AEQD paths.

The function now auto-projects lon/lat input to AEQD internally
(matching the convention of `mt_flag_outliers_bridge` and
`mt_flag_outliers_dbgb`) and computes geometry on a single canonical
projection.  Output `persistence_count` and `persistence_at_scale_*`
columns are now byte-identical regardless of whether the caller
supplies lon/lat or projected input.  No user-facing API change.

Internal helper `.scale_k_geometry()` simplified accordingly: the
lon/lat-vs-projected dispatch is removed (always sees projected
input from the caller); `.bearing_pair()` (only used by the lon/lat
branch) is removed.

Tests: 1 new in `tests/testthat/test-mt_persistence_score.R`
asserting `persistence_count` is identical for lon/lat and AEQD
inputs on CPF_D.

## Asymmetric pre-peel as opt-in on `mt_clean_track()` (2026-05-11)

**New parameter `pre_peel_aux = c("none", "primitives")`** on
`mt_clean_track()`.  Default `"none"` preserves the symmetric
pre-peel behaviour exactly.  Setting `pre_peel_aux = "primitives"`
runs `mt_flag_outliers_bridge()` and `mt_flag_outliers_detour()`
once on the raw input, builds a per-fix auxiliary score by
rank-normalising and summing their magnitudes, and passes it to
`mt_peel_speed()`'s `aux_scores` argument.  The peel then flags only
the higher-scoring endpoint per offending edge instead of both
endpoints — sparing clean neighbours of 1-fix spikes.

**Empirical motivation** (`CASCADE_AUDIT_2026-05-11.md` Section 6.1):
the symmetric pre-peel removes the spike *and* its two clean
neighbours on every 1-fix-spike pattern; round-3 #1 demonstrated that
asymmetric peel rescues those neighbours when an aux score
discriminates them from the spike.  Wiring the aux into
`mt_clean_track()`'s pre-peel closes the structural gap that round-3
#1 had only addressed at the primitive level.

**CPF benchmark with mass-mode = (1.0 kg, "flying")**:

| Track | symmetric F1 | asymmetric F1 | ΔF1 |
|---|---|---|---|
| CPF_A | 0.885 | **0.958** | **+0.074** |
| CPF_C | 1.000 | 1.000 | 0 |
| CPF_D | 0.674 | **0.789** | **+0.115** |
| CPF_E | 0.994 | 0.994 | 0 |
| CPF_F | 0.533 | 0.533 | 0 |

Mean ΔF1 +0.036.  Two clear wins (CPF_A halo, CPF_D block boundary
cleanup), no losses on synthetic.

**Cluster-outlier caveat — do NOT enable on K02-style sustained
spoofs.**  On a coherent multi-fix cluster the auxiliary score may
not discriminate between spoof-interior and clean-track-interior
boundary fixes; the asymmetric peel can walk inward from one side
only or oscillate, producing FPs at the cluster boundary.  K02
benchmark: F1 0.776 (symmetric) → 0.572 (asymmetric).  Default
`"none"` is the right choice on data with known cluster-shape
contamination (Argos PTT spoofs, deployment confusion blocks).
Documented in the `pre_peel_aux` @param section of
`?mt_clean_track`.

Tests: 4 new in `tests/testthat/test-mt_clean_track.R` covering
default = "none" backward-compat, CPF_A win, CPF_D win, and
no-op-without-v_max.

## Multi-scale persistence annotation (2026-05-09 evening)

**New function `mt_persistence_score()`** is a detector-agnostic
annotation helper.  It takes the output of any flagger that emits an
`is_outlier` column (`mt_clean_track`, `mt_flag_outliers_bridge`,
`mt_flag_outliers_detour`, `mt_flag_outliers`, `mt_flag_speed_cap`,
`mt_sequential_outliers`, `mt_combined_outliers`) and adds a
per-flag confidence score quantifying how anomalous each flag's
local geometry remains when viewed over wider temporal windows
(default `scales = c(2, 4, 8)`).  Adds columns
`persistence_count` (1 = flagged only at native resolution; up to
4 = flagged at every validation scale) and one boolean column per
scale.  `is_outlier` is never modified.

The mechanism is *windowed*, not thinned: at each scale `k`, the
function computes a scale-`k` view of *every* fix's local geometry
(step lengths and turn over `i±k`-fix windows) and asks whether
the flagged fix's scale-`k` geometry is anomalous against the
track's full scale-`k` distribution.  No thinning, no index-parity
bias.

**Empirical class-conditional finding** (synthetic CPF + cohort
analysis 2026-05-09): on cascade output, persistence cleanly
discriminates true positives from false positives within the
`state_anomaly` class (TPs persist 96% at p≥3 vs FPs at 57%, gap
+39 pp) and the `consensus` class (82% vs 44%, gap +38 pp).  The
`geometric_spike` class is empirically pure on the synthetic --
no FPs to filter.  The recommended pattern is therefore *class-
aware filtering*: gate any persistence-based filter on
`error_class %in% c("state_anomaly", "consensus")`.  The vignette
`vignette("persistence_score")` documents the analysis and the
usage pattern.

**Retirement of `mt_flag_outliers_multiscale()`.**  The previous
multi-scale function thinned the track at multiple scales and
voted across thinned-detector outputs, which has a structural
index-parity bias (a fix at original index 7 with
`scales = c(1, 2, 4, 8)` is in the thinned subset only at scale 1
and can vote at most once -- so even-indexed fixes could never
reach the default `min_votes = 2`).  The empirical work that
motivated the redesign is recorded in `FUTURE_IMPROVEMENTS.md`
"Round-3 follow-up -- Multi-scale persistence redesign".  The
retired source is preserved at
`archive/2026-05-09-multiscale-retired/` for historical reference;
it is no longer exported, no longer documented, and no replacement
shim is provided -- users who need the old behaviour should pin
to v0.1.x.

## Unified state-aware bridge primitive (2026-05-07)

**New function `mt_flag_outliers_dbgb()`** is now the canonical state-
aware bridge primitive.  It estimates per-axis (parallel /
orthogonal) motion variance via `mt_dbgb_variance()` and applies an
**envelope flag rule** across three Z channels:

- `bridge_z_para` -- residual component along the local travel
  direction normalised by the parallel-axis sigma; |Normal(0,1)|
  under the bridge null.
- `bridge_z_orth` -- across-axis component normalised by the
  orthogonal-axis sigma.
- `bridge_z_chisq = sqrt(Z_para^2 + Z_orth^2)` -- the diagonal-Sigma
  Mahalanobis magnitude; Rayleigh(1) under the null.  Reduces to the
  dBBMM scalar Z when sigma_para = sigma_orth.

A fix is flagged when **any** of the three channels exceeds its own
threshold, plus a `bridge_z_class` diagnostic label tagging each
flagged fix with which channel(s) fired -- including chisq-only
sub-classification by axis dominance (`isotropic_para_dom`,
`isotropic_orth_dom`, `isotropic`) so directional information is
recoverable even when only the joint Z fires.

**Threshold methods.** `z_threshold_method = "bonferroni"` (default)
applies per-channel parametric thresholds with joint FWER <= 0.05
(alpha/2 to chisq, alpha/4 to each per-axis channel).  `"bh_fdr"`
applies the same alpha allocation under BH-FDR.  `"gap"` runs the
package's broken-stick + tail-decay inflection detector per channel
(self-thresholding; useful for diagnostic exploration but mis-
calibrates on raw Z whose null is half-normal / Rayleigh with thin
natural tails -- see `?mt_flag_outliers_dbgb` Details).

**Performance.** Underlying dBGB variance estimation runs in native
C (`src/bgb_var_window_c.c`), bringing it to ~1.7x the cost of dBBMM
(was ~200x slower under the previous all-R implementation with one
`optim()` call per breakpoint candidate per window).  Numerical
equivalence with the original R reference is asserted in
`tests/testthat/test-bgb_var_break_c_vs_r.R`.

**`mt_flag_outliers_dbbmm()`** is retained as a thin wrapper around
`mt_flag_outliers_dbgb(variance = "dbbmm", ...)` -- the
sigma_para = sigma_orth submodel, faster (one 1-D Brent fit per
window rather than two) and emitting a simpler scalar-Z output
schema.  Equivalent to calling the dbgb function with the variance
constraint; recommended only when the per-axis decomposition isn't
needed.

**`mt_flag_outliers_bridge()` (legacy, leverage-immune) is
unchanged** -- it remains the static bridge primitive that doesn't
estimate motion variance and is appropriate when leverage protection
is needed.  `mt_clean_track()`'s cascade still uses it for the
bridge-residual stage.

## Naming consistency: obs_error -> location_error (2026-05-07)

**`obs_error` parameter renamed to `location_error`** in
`mt_flag_outliers_bridge()` and `mt_clean_track()`, plus the
internal helper file `R/obs_error.R` -> `R/location_error.R` with
matching function renames (`.resolve_obs_error` ->
`.resolve_location_error`, etc.).  Brings the package onto a single
canonical name; previously the variance estimators
(`mt_dbbmm_variance()`, `mt_dbgb_variance()`) and the UD
constructors used `location_error` while the bridge primitive used
`obs_error` -- an internal inconsistency that confused users mixing
primitives.  Tests updated.  No deprecation shim provided.

## Bug fixes (2026-05-07)

**`mt_flag_outliers(method = "copula")` recall regression** -- the
von Mises fit was being applied to angular *velocity* (rad/s)
instead of raw turn angles.  Under `time_normalize = TRUE` with the
C1 unit fix giving proper seconds, angular-velocity values compress
proportionally and `dvonmises()` numerically rounds to zero
everywhere off the mean -- driving every fix's joint probability to
zero and disabling outlier detection.  Fixed by fitting the von
Mises to the raw turn angle (circular, scale-invariant) while
keeping the Weibull on speed.  Recall on the synthetic benchmarks
restored.

**`time_normalize = FALSE` "backward compat" test removed.**  The
test asserted that the default call equals an explicit
`time_normalize = FALSE`; with the audit pass's C1 unit fix this
correctly fails (the default is `time_normalize = TRUE`, so the two
paths must differ by definition).  The test only ever passed
because the C1 unit bug had silently zero-broadcast the time
normalisation.

## Breaking changes

**`mt_emd()` renamed to `emd()`.**  The `mt_` prefix was dropped
because the function operates on UD rasters, not on `move2`
trajectory objects, so the prefix conveyed no useful dispatch
signal.  Update call sites; no deprecation shim is provided.

**`mt_clean_track()` now uses class-aware flagging by default**
(`flagging_strategy = "class_aware"`).  Set
`flagging_strategy = "legacy_strict"` to reproduce the v0.1.x behaviour
exactly.

The previous strict conjunction rule
`(bridge OR detour) AND (prob OR speed)` rested on the assumption
that a geometric flag must be confirmed by a kinematic detector.
That assumption fails on **symmetric out-and-back GPS spikes at
sparse sampling**: at 1-h GPS the spike's implied step speed is well
below physiological caps (so speed cannot agree), and the autodiff
framework is structurally blind to symmetric returns (so prob cannot
agree either).  The geometric primitives fire correctly, but the
conjunction rule discards them — exactly the case the package should
catch.

The new class-aware strategy decomposes flagging into independent
class rules, each using only the tools that are *in scope* for that
class:

  - `consensus`              : ≥ 3 of {bridge, detour, prob, speed}
  - `geometric_spike`        : bridge AND detour (no kinematic needed
                                — geometric impossibility is sufficient)
  - `state_anomaly`          : (bridge|detour) AND speed_state
  - `kinematic_confluence`   : (bridge|detour) AND prob

A fix is flagged if ANY class fires.  The output `error_class` names
the highest-confidence class that fired.  Single-detector fires never
flag (preserving precision).  This avoids one tool's *weakness*
disqualifying another tool's *strength*.

**Verification on synthetic ground truth**: CPF_A 23/3/0
(TP/FP/FN; one extra FP vs v0.1.x's 23/2/0), CPF_B 0/0/0 (no false
positive on clean data preserved), CPF_C 4/0/0 (unchanged), K02
spoof: 1162 flagged (unchanged).

### Class-aware vs strict head-to-head on every CPF track (CASCADE_AUDIT_2026-05-11 §5)

The 2026-05-11 cascade orchestrator audit measured F1 on each
CPF track under the new `class_aware` default and the four
`legacy_strict` sub-modes.  The headline finding is that
`class_aware` is **structurally necessary** on at least one CPF
pattern, not stylistic:

| Strategy | CPF_A | CPF_C | CPF_D | CPF_E | CPF_F |
|---|---|---|---|---|---|
| **class_aware** (default) | **0.958** | **1.000** | **0.674** | **0.994** | **0.533** |
| legacy_strict (0.1.x) | 0.958 | 1.000 | 0.674 | **0.140** ⚠ | 0.533 |
| legacy_strict, majority | 0.920 | 1.000 | 0.674 | 0.976 | 0.035 |
| legacy_strict, speed_trusted | 0.807 | 0.727 | 0.667 | 0.140 | 0.015 |
| legacy_strict, any | 0.123 | 0.190 | 0.161 | 0.429 | 0.014 |

**CPF_E (halo-spike pattern) is where class_aware earns its keep**:
F1 = 0.994 under class_aware vs **0.140 under strict** — a
catastrophic strict under-flag of single-fix spikes whose kinematics
are too subtle for prob/speed confirmation under the strict
conjunction's required-confirmer rule.  The geometric_spike class
(bridge AND detour, no kinematic confirmer needed) catches these
correctly.

The 2026-05-06 architecture switch was not a refactor; it was a
substantive precision-recall correction.  Backward-compat is
preserved via `flagging_strategy = "legacy_strict"` for any caller
that needs to reproduce 0.1.x exactly.

### Cohort verification: synthetic-truth proxy on Saline (CASCADE_AUDIT_2026-05-11 §3.2.1)

On real cohort tracks ground truth doesn't exist, so a flag-fraction
shift from 0.1.x to 0.2.0 is interpretive on its own.  The audit
removed the ambiguity for the Saline (homing pigeon) track by
injecting 25 controlled synthetic spikes and running both flagging
strategies under the SAME working-tree code (isolating the strategy
change from any other 0.1.1 → 0.2.0 drift):

| Strategy | n_flag | TP | FP | FN | F1 |
|---|---|---|---|---|---|
| class_aware / auto-cap | 47 | 25 | 22 | 0 | **0.694** |
| legacy_strict + strict / auto-cap (0.1.x default) | 98 | 25 | 73 | 0 | 0.407 |

Both achieve recall 1.0 (all 25 truth caught); class_aware adds 51
fewer FPs.  Genuinely measurable precision improvement on this
cohort track + spike pattern; not silent regression.

### Why this is the right architecture (paper-relevant)

Each tool measures a different facet of "impossible" with a precise
scope of validity.  Bridge measures perpendicular geometric deviation
in metres; detour measures dimensionless path-vs-displacement ratio
(scale-invariant); speed measures step velocity (scale-dependent);
autodiff measures kinematic *change* between consecutive steps (blind
to symmetric returns).  A principled cleaner runs each tool where it
is in scope and combines evidence *within each class*, not across
classes — preventing a tool's structural blind spot from gating out
another tool's signal.

## Combine-with-global pass for state-conditional dispatch

When `state =` is supplied, `mt_clean_track()` now ALSO runs a
single-state baseline on the full input and unions the
per-detector fingerprints with the per-segment results.  Geometric
impossibility (bridge AND detour both firing) is state-
INDEPENDENT, so the geometric primitives' data-driven thresholds
should be calibrated to the full-track residual distribution
rather than per-segment.  Inside a migrate segment the bridge
threshold is inflated by legitimate flight residuals, so colony
spikes embedded in migrate sub-trains slip through; the global
baseline rescues them.  Per-state thresholds remain in force for
kinematic primitives (speed, prob) where state-relative tightness
is the right scope.

Empirical effect on Víctor's gull track: 252 `geometric_spike`
flags (vs 80 without the global pass) — recovers all 190 colony
halo fixes that single-state class_aware catches AND keeps the 214
wintering-area state-anomaly fixes that per-state alone catches.

## New parameters on `mt_clean_track()`

* `flagging_strategy` — `"class_aware"` (default) or `"legacy_strict"`.
* `transition_buffer` — non-negative integer (default `1L`).  Width of
  the buffer zone around state transitions inside which state-
  dependent flag classes (`state_anomaly`, `kinematic_confluence`)
  are demoted to `"state_transition_buffered"` (kept).  Geometric
  consensus and block expansion still flag at transitions because
  geometric impossibility doesn't depend on which state assignment
  we trust.

## New `error_class` taxonomy (class_aware mode)

* `physiological`            — pre-peeled (flag_iteration == 0)
* `block`                    — block expansion (block_id non-NA)
* `consensus`                — ≥ 3 of 4 detectors fired
* `geometric_spike`          — bridge AND detour, < 3 detectors
* `state_anomaly`            — (bridge|detour) AND speed, no prob
* `kinematic_confluence`     — (bridge|detour) AND prob, no speed
* `state_transition_buffered`— kept (not flagged); the fix would have
                                fired a state-dependent class but is
                                in the transition buffer

# move2utils 0.1.1

## New

* `mt_flag_outliers_detour()` — fourth point-level outlier-detection
  primitive. Flags fixes with anomalous path-vs-displacement ratio
  over a window of k fixes:
  `ratio_k(N) = path_length(N-k..N+k) / displacement(N-k, N+k)`.
  Time-insensitive (uses only successive coordinates) and scale-
  invariant (dimensionless ratio), complementing the Brownian-bridge
  perpendicular residual where bridge sigma-scaling loses sensitivity
  at sparse sampling rates. Catches out-and-back GPS spikes whose
  implied step speed stays below physiological caps -- the dominant
  failure mode at the colony in 1-h sampled tracks.

* `mt_clean_track()` — detour primitive integrated as a fourth peer
  detector via three new parameters: `use_detour = TRUE` (default),
  `detour_k = 1L`, `detour_threshold = 8`. Strictness rules
  generalised to four detectors:
  - `strict`: `(bridge OR detour) AND (prob OR speed)`
  - `majority`: ≥ 2 of `(bridge, detour, prob, speed)`
  - `speed_trusted`: `speed OR ((bridge OR detour) AND prob)`
  - `any`: union
  When `use_detour = FALSE` the rules collapse to the earlier
  three-detector form; existing user code is unaffected. The default
  `detour_threshold = 8` is the minimum value at which the synthetic
  ground truth (CPF_A/B/C) is preserved exactly under the new
  rule (no false positive on CPF_B, no recall loss on CPF_A or
  CPF_C). New `error_class` value `"detour"` for fixes flagged only
  by the detour detector.

* `mt_diagnose_clean_track()` — post-run diagnostic suite for
  `mt_clean_track()` results, analogous to `plot.lm()` for a GLM.
  Six panels: log-speed density with detected modes and the v_max
  used; rolling flag rate over time; per-detector activity counts;
  cumulative flagging by iteration; consecutive-flag run-length
  histogram; per-individual flag-rate dotplot for multi-track
  results. Auto-prints concerns it identifies (multi-mode behaviour,
  sustained migration band, single-detector dominance, non-converging
  iteration, heavy run-length tail) with actionable verdicts pointing
  the user at remedies. Companion vignette
  `vignette("diagnose_clean_track")` walks through the panels on a
  contrasted PASS / OVER_FLAG pair.

* `mt_flag_outliers_dbbmm()` — new (work-in-progress) state-aware
  companion to `mt_flag_outliers_bridge()`. Normalises the
  Brownian-bridge residual by the local Brownian motion variance from
  `mt_dbbmm_variance()` / `mt_dbgb_variance()`. The resulting Z-score
  follows Rayleigh(1) under the bridge null, so the threshold is
  principled (Bonferroni FWER = 0.05 by default). Targets the
  documented multi-state failure mode the legacy bridge primitive
  exhibits on tracks with mixed behavioural states. Trade-off: gains
  state-awareness, exposes the leverage problem (extreme outliers
  inflate sigma^2 then bound their own Z); mitigated via
  `pre_peel_v_max` parameter. Not yet wired into `mt_clean_track()`.

* `mt_suggest_dbbmm_window()` — new exported diagnostic helper that
  inspects a track's temporal sampling and suggests `window_size` and
  `margin` for `mt_dbbmm_variance()` / `mt_dbgb_variance()`. Picks
  `window_size` as the nearest odd integer covering a user-specified
  target time per window (default 4 h), clamped to the literature
  floor of 11 and to `floor(n/4)` so several non-overlapping windows
  fit in the track. Multi-track aware; prints summary and optionally
  plots the time-lag distribution. Returns an
  `mt_dbbmm_window_suggestion` S3 with `print()` method.

* `mt_flag_outliers_bridge(obs_error = ...)` — new parameter ingesting
  per-fix observation-error priors at the bridge anchors. Accepts a
  numeric scalar (uniform sigma in m), a per-fix vector, a column name
  (e.g. `"eobs_horizontal_accuracy_estimate"`), or `"auto"` to probe
  Movebank quality columns and the Argos location-class table
  (Vincent et al. 2002 ladder). The bridge denominator absorbs anchor
  variance via `w_eff^2 = w^2 + (w_p^2 sigma_prev^2 + w_n^2 sigma_next^2)
  / S_hat`, where `S_hat = median(r^2/w^2)` is a distribution-free
  empirical residual scale estimated on the active mask in iteration 1.
  The target fix's own sigma deliberately does not enter, preserving
  leverage immunity. Most useful on Argos / mixed-mode tracks where
  per-fix accuracy varies by orders of magnitude. New diagnostic
  column `bridge_obs_inflation = w_eff/w` (== 1 when injection is off
  or all anchors have unknown sigma). `mt_clean_track(obs_error = ...)`
  forwards through to the bridge primitive. Default `NULL` preserves
  existing behaviour.

* `mt_peel_speed()` — new exported primitive for iterative speed peel
  at a user-supplied physiological cap. Removes fixes with
  `max(step_in, step_out) > v_max`, recomputes step speeds on the
  survivor set, and repeats until convergence. Catches coherent
  multi-fix error clusters (e.g. GPS spoofs with internally-consistent
  fake trajectories) that per-fix detectors cannot resolve, because
  peeling a boundary fix exposes the next interior fix as a new
  boundary with an impossible step speed. Fully vectorised; multi-track
  dispatch; auto-projection for longitude/latitude input. Adds
  `is_outlier`, `peel_iteration`, `step_speed` columns and attributes
  `v_max_used`, `n_peel_iterations`, `converged`.

* `v_phys_estimate(mass, mode)` — new exported helper. Returns a
  body-mass / locomotor-mode estimate of the species' physiological
  maximum burst speed using the Hirt et al. (2017, *Nat. Ecol. Evol.*,
  doi:10.1038/s41559-017-0241-4) general scaling law. Mode is one of
  `"flying"`, `"running"`, `"swimming"`. Output is a printable scalar
  in m/s with a parameter-uncertainty confidence interval propagated
  from the published parameter standard errors via the delta method.
  Designed to be passed directly into `mt_clean_track()` as the
  principled `v_max` default when no species-specific published
  maximum is at hand.

* `mt_clean_track()` — new `state =` parameter for state-conditional
  cleaning. Accepts `NULL` (default; current behaviour), a column name
  on the move2 object, or a per-fix vector. Partitions each track into
  contiguous runs of constant state and runs the full pipeline
  independently on each segment, so the per-fix detectors see a single
  kinematic distribution per segment (avoids the "cut between modes
  vs. cut beyond modes" tradeoff on multi-state tracks). Segments
  shorter than 3 fixes pass through unflagged; `NA` is treated as its
  own state value. Segment-local `block_id`s are renumbered globally
  on recombination. The package's contract is to *respect*
  user-supplied state labels — segmentation itself (speed-threshold,
  HMM, BCPA, manual annotation) is the user's responsibility, with
  worked examples in the rewritten `state_conditional` vignette.

## Vignettes

* `vignette("diagnose_clean_track")` — six-panel post-run health check
  walked through on a contrasted PASS (CPF_B synthetic) / OVER_FLAG
  (Pettstadt1 white stork) pair, with a decision-tree table mapping
  each diagnostic signature to its remedy.
* `vignette("state_conditional")` — rewritten around the new `state =`
  parameter on `mt_clean_track()`. Three usage patterns (column on
  the move2 object / speed-threshold-derived vector / HMM-decoded
  state) plus a worked synthetic example and limitations
  (segmentation responsibility, endpoint effects, boundary outliers).

## Empirical evidence

* `audits/2026-04-29-stratified-movebank/` — 65-individual stratified
  Movebank audit across 4 movement classes (storks, marine birds,
  raptors, terrestrial mammals) and 22 studies. With naive defaults
  and *no* per-track tuning, 75% PASS at <0.5% flag rate, 95%
  defensible, 0 pipeline crashes. The 5% OVER_FLAG cases all hit the
  documented multi-state-with-thin-activity-mode failure and were
  identified correctly by `mt_diagnose_clean_track()`. A parameter-
  ladder follow-up (tier 1 mass+mode, tier 3 obs_error, tier 6
  no-block-expansion) shows each adjustable targets a specific
  empirically-distinguishable failure; mass+mode rescues thin-activity
  cases but can regress multi-state ones, so the diagnostic-driven
  user workflow is the supported path.

## Changed

* `mt_clean_track()` — substantially faster on lon/lat tracks
  (~2.8x speedup on a 855k-fix track in profiling). Three changes
  combined:
  (1) the per-iteration `mt_distance()` calls in
  `mt_flag_outliers()`, `mt_flag_speed_cap()`, and the internal
  block-expansion step now use a vectorised Haversine
  (`.step_lengths_fast()`), about 60x faster than s2 on lon/lat with
  < 2 m worst-case error vs. the WGS84 ellipsoidal reference;
  (2) the per-iteration `mt_turnangle()` call in `mt_flag_outliers()`
  now uses a direct spherical/Euclidean azimuth helper
  (`.turn_angles_fast()`), about 25x faster than the `s2` path,
  with < 0.003 rad max divergence from the ellipsoidal reference;
  (3) `mt_clean_track()` pre-projects lon/lat input to AEQD once,
  rather than letting the bridge primitive auto-project on every
  iteration. Behavioural sanity preserved: identical TP/FP/FN on
  CPF_A/B/C synthetic ground truth (CPF_A 23/23 truth + 2 FPs;
  CPF_B clean; CPF_C 4/4 truth + 0 FPs).

* `mt_flag_outliers()` — accepts projected input directly. Earlier
  the function aborted on projected input because
  `move2::mt_turnangle()` calls `mt_azimuth()`, which currently only
  supports lon/lat. With the new internal `.turn_angles_fast()`
  helper, projected input runs natively (Cartesian `atan2(dx, dy)`)
  with no transform.

* `mt_clean_track(v_max = ...)` — when the user supplies `v_max`,
  the pipeline now calls `mt_peel_speed()` as a pre-step, then runs
  the existing per-fix conjunction loop on the peeled remainder with
  internal speed-cap set to `auto` (preventing the cascading that the
  old hard-cap conjunction + block expansion produced on heavy-tail
  data: K02 with `v_max = 30` previously over-flagged 72 % of the
  track; now produces a clean 0.06 % result with 100 % spoof
  coverage). Block expansion is skipped when pre-peel ran, since the
  peel's iterative nature already handles cluster propagation. Default
  behaviour (`v_max = NULL`) is unchanged — identical to 0.1.0.

* `mt_suggest_speed_cap()` — now reports a **congruence diagnostic**:
  runs both entropy- and gap-based break detectors and reports the
  disagreement magnitude between the two candidates. Large disagreement
  (> 3×) hints at multi-scale tail structure (coherent outlier cluster
  plus extreme singletons) that warrants user inspection before setting
  `v_max`. Does not change the returned value — purely diagnostic.

* `mt_clean_track(mass, mode)` — accepts the same `(mass, mode)` pair
  as `v_phys_estimate()`; when both are supplied (and `v_max` is
  `NULL`), the function derives the cap allometrically and uses it
  for the iterative speed peel. Mutually exclusive with `v_max`.

* `mt_suggest_speed_cap(mass, mode, v_phys)` — extended with three
  optional arguments that overlay independent perspectives on the
  diagnostic plot. Supplying `(mass, mode)` adds the allometric
  prediction; supplying just `mass` shows all three mode-specific
  lines so the user can compare predictions for ambiguous cases
  (mode-switching species, uncertain locomotor categorisation);
  supplying `v_phys` adds the user's own estimate. The function
  prints a triangulation diagnostic — empirical break vs allometric
  prediction vs user estimate — that flags contamination (empirical
  > allometric) and within-distribution structure (empirical <
  allometric) explicitly.

* `mt_clean_track()` — adds an `error_class` column when
  `remove = FALSE` that maps the per-detector flag combination to a
  categorical interpretation of *why* each fix was caught:
  `"block"` (member of a topologically-isolated cluster, typical of
  GPS spoofs and timestamp glitches), `"consensus"` (≥ 2 detectors
  agree on the same fix; highest confidence), `"speed_cap"`
  (speed-only; physical impossibility on an isolated transition),
  `"jitter"` (bridge-only; geometric outlier without kinematic or
  speed violation, typical of GPS multipath), or `"kinematic"`
  (probability-only; behaviour-metric anomaly too subtle to disturb
  the bridge or speed-cap detectors). `NA` where `is_outlier` is
  `FALSE`. The categorical column complements the existing
  `flagged_by_bridge` / `flagged_by_prob` / `flagged_by_speed`
  Boolean provenance columns.

* `mt_flag_outliers_bridge()` — default diagnostic plot now leads
  with the directional decomposition: a
  `log10(eta_para + 1)` vs `log10(eta_perp + 1)` scatter with the
  diagonal `y = x` reference line that separates isotropic jitter
  from along-track and across-track errors. Flagged fixes
  highlighted in red; the existing track map remains as the right
  panel. Falls back to the old sorted-eta plot when only the
  isotropic score is available.

* `mt_flag_speed_cap()` — default diagnostic plot now leads with the
  step-speed empirical CDF on log-x with the applied `v_max` line
  and a `n above (frac %)` annotation that makes the cost of the
  chosen cap visible. The existing track map remains as the right
  panel.

## Internal

* Removed stale `.compute_bridge_residuals()` helper (the isotropic-only
  version superseded by `.compute_bridge_residuals_dbgb()`).

# move2utils 0.1.0

Initial release.

## Outlier detection

* `mt_flag_outliers()` — probability-based outlier detection using the
  empirical joint distribution of speed, angular velocity, and their
  gap-aware auto-differences. Four threshold types: `"gap"`,
  `"entropy"`, `"significance"`, `"percentile"`. Multi-individual
  dispatch, iterative refinement, reference-distribution support, and
  GPS quality weighting.
* `mt_flag_outliers_bridge()` — geometric outlier detection based on
  the time-weighted Brownian-bridge mean of each location's temporal
  neighbours. Method options `c("combined", "isotropic", "directional")`;
  the directional variant decomposes residuals into parallel and
  perpendicular components for error-morphology classification.
  Leverage-immune by construction (no locally-estimated variance as
  denominator). Auto-projects longitude/latitude input to a local AEQD
  for the Euclidean math and returns the result in the caller's
  original CRS. Optional `residual_floor` parameter gates rate-flagged
  fixes by a minimum absolute residual (metres) — prevents GPS-noise
  jitter in burst-sampled segments from being flagged as outliers.
* `mt_flag_speed_cap()` — step-level speed cap. The only primitive
  that scores steps (not fixes): catches the boundary transitions of
  coherent outlier blocks that per-fix detectors structurally cannot
  resolve. `threshold_type = "auto"` (default) uses a dip-test-validated
  entropy-or-gap break on `-log(speed)` with a 5% safety guard that
  rejects bulk-mode separators. Hard physiological cap via
  `threshold_type = "hard"` and user-supplied `v_max` in m/s.
* `mt_suggest_speed_cap()` — diagnostic that prints speed quantiles
  and suggests a data-driven `v_max` from the structural break in the
  speed distribution. Returns `NA` with an actionable message when the
  distribution has no meaningful outlier tail (pointing at the
  hard-cap route as the appropriate alternative).
* `mt_clean_track()` — unified one-call entry-point that iterates the
  three primitives with a per-fix conjunction rule (configurable via
  `strictness`) plus a topological block-expansion step that
  disentangles disconnected components of the survivor graph. Iterates
  until convergence (`iterations = "until_clean"` default, safety cap
  200) with three stop criteria (no-new-flags, flag-fraction breach,
  max-iterations). Returns the input object with per-detector
  diagnostic columns (`flagged_by_bridge`, `flagged_by_prob`,
  `flagged_by_speed`, `block_id`, `flag_iteration`).
* `mt_flag_outliers()` — probabilistic detector gains an optional
  `step_floor` parameter that complements `residual_floor` on the
  bridge side: gate rate-flagged fixes by a minimum absolute step
  length.
* `mt_filter_gps_quality()` — new `drop_empty = TRUE` default drops
  rows with empty point geometries as part of the standard hygiene
  pass.
* `mt_sequential_outliers()` — sequential scan with three strategies
  (`"forward-backward"`, `"greedy"`, `"random"`) for clusters of
  consecutive outliers.
* `mt_combined_outliers()` — majority vote across gap, entropy, and
  sequential methods.
* `mt_flag_outliers_multiscale()` — multiscale detection with voting
  across temporal resolutions.

## Utilisation distributions and variance estimation

* `mt_dbbmm_variance()` / `mt_dbbmm_ud()` — dynamic Brownian-bridge
  movement model variance and utilisation distribution on a common
  grid. Ported from the `move` package C kernels to a modern
  `move2` / `sf` / `terra` spatial stack.
* `mt_dbgb_variance()` / `mt_dbgb_ud()` — directional (bivariate
  Gaussian bridge) variance and utilisation distribution.
* `mt_motion_variance()` — accessor returning estimated variances as a
  numeric vector or data frame.
* `ud_volume()` — cumulative-volume transform of a UD raster
  (analogue of `move::getVolumeUD()`); operates on single- or
  multi-layer `terra::SpatRaster`.
* `emd()` — pairwise Earth-mover's distance between UDs (port of
  `move::emd()`, Kranstauber, Smolla & Safi 2016). Default uses
  Sinkhorn entropic regularisation (Cuturi 2013) with a volume-based
  pre-mask, giving hundreds-of-times speedup over the legacy
  simplex-LP solver at sub-percent accuracy. Exact mode via the
  optional `emdist` package retained for validation.
* `ud_outer_probability()` — cumulative-volume quantile at a set of
  query locations (analogue of `move::outerProbability()`). Thin
  wrapper over `ud_volume()` + `terra::extract()`.
* `mt_thin_distance()` — cumulative along-track distance thinning
  (analogue of `move::thinDistanceAlongTrack()`). Multi-track
  dispatch; `thin_selected` column labelling convention matches
  `mt_thin_time()`.
* `mt_corridor()` — corridor-segment detection
  (port of `move::corridor()`, modern `sf` R-tree neighbourhood
  search, LaPoint et al. 2013).
* `mt_thin_time()` — strict tolerance-constrained time-thinning via a
  linear-time DP sweep; scales to \eqn{10^5}{1e5}+ fixes.

## Design notes

* Non-parametric gap-aware auto-difference normalisation makes
  detection robust to irregular sampling without assuming a movement
  model.
* ACF-derived `alpha` is the default weighting between the step-turn
  and auto-difference components of the joint probability.
* Entropy-valley threshold is a no-op on clean, unimodal data: the
  bridge primitive returns zero flags when no valley exists.
* Hartigan's dip test (`diptest` package, Suggests) validates gap-
  threshold candidates on the speed distribution. Rejects bulk-mode
  separators on real animal tracks (where most fixes are stationary
  or slow-moving) by refusing to accept a cap that would flag more
  than 5% of the distribution as outliers.
