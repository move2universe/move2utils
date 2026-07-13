#' Unified outlier detection for a move2 track
#'
#' One-stop detector that fuses the three primitives
#' (\code{\link{mt_flag_outliers_bridge}},
#' \code{\link{mt_flag_outliers}},
#' \code{\link{mt_flag_speed_cap}}) into a single iterative pipeline
#' and adds a topological block-expansion step.  Each primitive works
#' at a different grain of analysis (point, segment, step) and catches
#' a different class of error; combining them via a conjunction rule
#' and iterating to convergence catches the long outlier trains that
#' per-fix scoring structurally cannot resolve.
#'
#' @details
#' **Per-iteration logic:**
#' \enumerate{
#'   \item Run each primitive in single-pass mode on the current
#'         active (unflagged) set:
#'     \itemize{
#'       \item \code{mt_flag_outliers_bridge()} -> point-level bridge
#'             residual. Default \code{method = "combined"},
#'             \code{threshold_type = "entropy"} (strict).
#'       \item \code{mt_flag_outliers()} -> segment/vertex-level joint
#'             probability. Default \code{threshold_type = "gap"}
#'             (sensitive).
#'       \item \code{mt_flag_speed_cap()} -> step-level implied speed.
#'             Uses \code{v_max} if supplied; otherwise runs the
#'             detector in \code{threshold_type = "auto"} mode (entropy
#'             valley with dip-test-validated broken-stick fallback).
#'     }
#'
#'   \item Apply the **conjunction rule** for per-fix flagging: fix
#'         \eqn{i} is flagged if the bridge says its point-level
#'         residual is suspect AND at least one of (probability,
#'         speed-cap) also says an incident transition is suspect.
#'         The conjunction rewards agreement across grains and
#'         prevents any single score from over-flagging.
#'
#'   \item **Block expansion.** After flagging, partition the kept
#'         fixes into connected components: two kept fixes \eqn{i}
#'         and \eqn{j} are in the same component only if the step
#'         \eqn{i \to j} has an implied speed below the block-speed
#'         threshold (either \code{v_max} or whatever the speed
#'         detector chose in \code{"auto"} mode).  Small components
#'         (size \eqn{< \text{max\_flag\_fraction} \cdot n}) that are
#'         disconnected from the main trajectory through suspect
#'         transitions are flagged as blocks.  This is the step that
#'         dissolves long trains the per-fix detectors structurally
#'         cannot catch --- the train interior is locally coherent
#'         but topologically isolated after stage 1 flags its
#'         boundaries.
#'
#'   \item **Stop criteria** checked at the end of every iteration:
#'     \itemize{
#'       \item No new fixes flagged -> converged.
#'       \item Cumulative flags exceed \code{max_flag_fraction} of
#'             the track -> abort with warning (likely a pipeline
#'             misuse).
#'       \item Iteration count reaches \code{max_iterations}.
#'     }
#' }
#'
#' **CRS handling:** like the other detectors, this function
#' auto-projects longitude/latitude input to a local AEQD for the
#' Euclidean math and returns in the original CRS.
#'
#' **Composition with the individual primitives:** the three
#' primitives remain exported for users who want to inspect which
#' score caught what.  \code{mt_clean_track()} is the convenient
#' one-call entrypoint.
#'
#' @param x A \code{move2} object.  Any CRS (auto-projected if needed).
#' @param v_max Numeric scalar or \code{NULL}.  Physiological speed
#'   cap in m/s.  If \code{NULL} (default), the speed detector runs
#'   in \code{"auto"} mode (dip-test-validated data-driven threshold)
#'   on every iteration.  Supply a positive scalar for a hard cap
#'   override.  Mutually exclusive with the \code{(mass, mode)}
#'   allometric route below.
#' @param mass,mode Optional pair.  When both are supplied (and
#'   \code{v_max} is \code{NULL}), the function derives a principled
#'   physiological cap from species body mass and locomotor mode via
#'   \code{\link{v_phys_estimate}} (Hirt et al. 2017 general
#'   scaling law) and uses the central estimate as \code{v_max}.
#'   This is the recommended route for users without a species-specific
#'   published maximum.  Specialist sprinters (cheetah, pronghorn) and
#'   3D-tracked diving / stooping data should override with a published
#'   value instead.
#'
#'   Both can be passed as either a single value applied to every
#'   individual, or a \emph{per-track} named vector whose names match
#'   the track ids.  Per-track \code{mass} is the right choice for
#'   multi-individual studies where individuals differ in body mass
#'   (juveniles vs adults, sexual dimorphism); per-track \code{mode}
#'   is rarely needed within a single-species study but supports
#'   mixed-locomotion datasets.
#'
#'   \code{mass} is in kg.  If \code{mass} carries a \code{units}
#'   attribute (e.g. as returned by Movebank's \code{animal_mass} in
#'   grams) it is auto-converted to kg.  A bare scalar \code{mass > 100}
#'   triggers a warning that the value looks like grams.  \code{mode}
#'   is one of \code{"flying"}, \code{"running"}, \code{"swimming"}.
#' @param state Optional behavioural-state assignment.  Default
#'   \code{NULL} (no segmentation; the cleaner sees one global
#'   distribution).  Accepted forms:
#'   \itemize{
#'     \item character of length 1: name of a column on \code{x}
#'           holding per-fix state values.
#'     \item vector of length \code{nrow(x)}: per-fix state values
#'           directly.  Numeric, character, factor, or logical are
#'           all accepted; \code{NA} is treated as its own state.
#'   }
#'   When supplied, the cleaner partitions each track into
#'   contiguous runs of constant state and runs the full pipeline
#'   independently on each segment.  This is the right behaviour for
#'   tracks with kinematically distinct states (rest + flight, perch +
#'   migration glide), where pooling the speed / residual
#'   distributions across states forces the threshold detectors to
#'   choose between over- or under-flagging the smaller mode.  The
#'   package's contract is to \emph{respect} user-supplied state
#'   labels --- segmentation itself (speed-threshold, HMM, manual) is
#'   the user's responsibility; see
#'   \code{vignette("OUTLIER_3_state_conditional", package = "move2utils")}.
#'
#'   Segments shorter than 3 fixes pass through unflagged (the per-fix
#'   detectors require at least 3 points for a residual / auto-
#'   difference); a one-line note is emitted unless \code{silent}.
#' @param consensus Character.  Which rule to use to combine the
#'   detector outputs into a per-fix outlier flag.  Defaults to
#'   \code{"evidence_corroborated"}, the likelihood-ratio rule: each
#'   detector's per-fix score is converted to a signed log-LR (surprisal
#'   minus its own data-driven flag boundary), the detectors are made
#'   commensurable and summed into one calibrated evidence score, and a
#'   fix is flagged when that evidence is positive AND it is either
#'   corroborated (at least two detectors agree) or one detector is
#'   overwhelming (its evidence is saturated -- far beyond its own
#'   distribution).  This is the empirically validated default
#'   (benchmark 2026-06-07: best F1 at canonical false-positive level,
#'   keeps the conspicuous-excursion catches on slow-species data, and
#'   recall 1.0 on spoofing/jamming).  Alternatives:
#'   \code{"class_aware"} (the previous Boolean class-rule default, still
#'   available), \code{"weighted_evidence"} (net evidence \eqn{> 0}
#'   without the corroboration safeguard), \code{"strict"},
#'   \code{"majority"}, \code{"speed_trusted"}, \code{"any"},
#'   \code{"custom"}.  See \code{\link{mt_flag_consensus}} for each mode
#'   and \code{DESIGN_evidence_accumulation.md} for the evidence
#'   framework.  The block-expansion step
#'   (when \code{expand_blocks = TRUE}) is independent of
#'   \code{consensus} -- it identifies disconnected components
#'   directly from the speed distribution, so coherent outlier blocks
#'   are caught regardless of the per-fix rule.
#' @param consensus_custom Function used when \code{consensus =
#'   "custom"}.  Receives four logical vectors of equal length
#'   (\code{by_bridge, by_prob, by_speed, by_detour}) and must return
#'   a single logical vector of the same length giving the per-fix
#'   outlier decision.  Ignored unless \code{consensus = "custom"}.
#'   See \code{\link{mt_flag_consensus}} for examples.
#' @param transition_buffer Non-negative integer.  Default \code{1L}.
#'   Width (in fixes) of the buffer zone around state transitions
#'   inside which state-dependent flag classes (\code{state_anomaly},
#'   \code{kinematic_confluence}) are demoted to
#'   \code{state_transition_buffered} (kept).  Set to \code{0L} to
#'   disable.  Geometric consensus (\code{geometric_spike},
#'   \code{consensus}) and \code{block} expansion still flag at
#'   transitions.  Only takes effect when \code{state =} is supplied.
#' @param use_detour Logical.  If \code{TRUE} (default), include the
#'   path-vs-displacement detour ratio (see
#'   \code{\link{mt_flag_outliers_detour}}) as a fourth point-level
#'   detector in the conjunction rule.  Detour is time-insensitive
#'   and scale-invariant, complementing the bridge primitive at
#'   sparse sampling rates where bridge \eqn{\sigma}-scaling loses
#'   sensitivity to single-fix spikes whose implied step speed stays
#'   below physiological caps.  Set to \code{FALSE} to reproduce the
#'   strict three-detector behaviour from earlier package versions.
#' @param detour_k Integer.  Window radius for the detour primitive.
#'   Default \code{1L} (single-fix spikes).  See
#'   \code{\link{mt_flag_outliers_detour}} for the multi-k diagnostic
#'   form.
#' @param detour_threshold Numeric, > 1.  Detour ratio threshold.
#'   Default \code{8}.  EMPIRICALLY TUNED -- the minimum value at
#'   which the synthetic ground-truth set (CPF_A/B/C) is preserved
#'   exactly under \code{use_detour = TRUE} (no false positive on
#'   CPF_B, no recall loss on CPF_A or CPF_C).  Plausible range:
#'   5--15; lower is more aggressive at sparse sampling.
#' @param iterations Either the string \code{"until_clean"} (default),
#'   \code{Inf}, or a positive integer.  \code{"until_clean"} and
#'   \code{Inf} are equivalent and mean "iterate until no new flags
#'   are produced or \code{max_iterations} is reached".  A positive
#'   integer caps iteration at that number of passes.
#' @param max_iterations Integer.  Hard safety cap on iteration count
#'   even with \code{iterations = "until_clean"}.  Default 100 is
#'   EMPIRICALLY TUNED (CASCADE_AUDIT_2026-05-11.md Section 3.4 +
#'   8): the iteration-count distribution across all audit tracks
#'   tops out at 22 iter for legitimate convergence on CPF / cohort;
#'   K02 (block contamination, auto-cap path) is the longest
#'   documented legitimate convergence at 61 iter.  100 covers K02
#'   with ~39-iter margin AND caps multi-state failure modes
#'   (WH17-class tracks that never converge) at half the previous
#'   wallclock.  Tracks that hit the cap return with
#'   \code{convergence = "max_iterations"} -- a visible signal the
#'   user can opt up via \code{max_iterations = 200L} (the previous
#'   default) or higher.  Plausible range: 50--500.
#' @param max_flag_fraction Numeric in (0, 1].  If cumulative flags
#'   exceed this fraction of the track, abort.  Default 0.2 is
#'   HEURISTIC -- the convention "if cleaning >20\% of the track,
#'   you are misusing the tool" rather than a derived bound.
#'   Plausible range: 0.1--0.3.  Lower values abort earlier on
#'   pathological tracks; higher values give the iteration loop more
#'   room.
#' @param expand_blocks Logical.  If \code{TRUE} (default), run the
#'   topological block-expansion step.
#' @param pre_peel_aux Character, one of \code{"none"} (default) or
#'   \code{"primitives"}.  Controls the pre-peel mode when
#'   \code{v_max} (or \code{(mass, mode)}) is supplied.  With
#'   \code{"none"}, the pre-peel is symmetric: every fix on either
#'   end of an offending edge (\code{step_speed > v_max}) is removed.
#'   With \code{"primitives"}, the function first runs
#'   \code{\link{mt_flag_outliers_bridge}} and
#'   \code{\link{mt_flag_outliers_detour}} once on the raw input,
#'   builds a per-fix auxiliary score by rank-normalising and summing
#'   their magnitudes, and passes it to \code{\link{mt_peel_speed}}'s
#'   \code{aux_scores} argument; the peel then flags only the
#'   higher-scoring endpoint per offending edge.
#'
#'   The asymmetric mode is designed for tracks whose contamination
#'   is dominated by 1-fix spikes whose clean neighbours the
#'   symmetric default would peel along with the spike.  Empirical
#'   benchmark (synthetic CPF, audit 2026-05-11): mean F1 +0.036,
#'   CPF_A F1 0.885 -> 0.958, CPF_D F1 0.674 -> 0.779, no losses.
#'
#'   \strong{Cluster-outlier caveat}.  On coherent multi-fix
#'   contamination (sustained spoofs, deployment confusion blocks
#'   like the K02 benchmark) the asymmetric peel can walk inward
#'   from one boundary only or oscillate, producing false positives
#'   at the cluster boundary.  Empirical (K02 benchmark, audit
#'   2026-05-11): F1 0.776 (symmetric) -> 0.572 (asymmetric).  Use
#'   the default \code{"none"} on tracks with known cluster-shape
#'   contamination.  See also \code{\link{mt_peel_speed}} -> Asymmetric
#'   peel.
#' @param persistence_filter Character, one of \code{"none"} (default)
#'   or \code{"class_aware"}.  When \code{"class_aware"}, runs
#'   \code{\link{mt_persistence_score}} on the cascade output and
#'   demotes flagged fixes in the \code{state_anomaly} +
#'   \code{consensus} error classes whose \code{persistence_count < 3}
#'   to \code{is_outlier = FALSE}.
#'
#'   \strong{Empirical motivation} (the 2026-05-09 class-conditional
#'   analysis): on cascade output, persistence cleanly separates TPs
#'   from FPs only on the \code{state_anomaly} (+39.3 pp gap at
#'   \eqn{p \geq 3}) and \code{consensus} (+37.9 pp) classes.  On
#'   \code{kinematic_confluence} the gap reverses (-14.3 pp);
#'   \code{geometric_spike} was class-pure on synthetic.  Filter is
#'   therefore class-aware by design -- universal application would
#'   actively regress halo-spike detection.
#'
#'   \strong{Per-track CPF effect under the CRS-invariant
#'   \code{mt_persistence_score} (post 2026-05-11 fix)}: CPF_A 0.958
#'   (unchanged); CPF_C 1.000 (unchanged); CPF_D 0.674 -> 0.675
#'   (+0.001); CPF_E 0.994 -> 0.975 (-0.019); CPF_F 0.533
#'   (unchanged).  Mean DF1 = -0.003.
#'
#'   Net-near-neutral on synthetic CPF: the class-conditional
#'   aggregate finding holds across many fixes but doesn't translate
#'   to per-track wins on the small validation set.  The filter
#'   ships as strict opt-in (round-3 mixed-CPF rule); the user-side
#'   case for enabling it is per-dataset.  The class taxonomy used
#'   by this filter (\code{error_class}) is independent of the
#'   \code{consensus} mode and is therefore always available.
#' @param location_error Per-fix observation-error prior (1-sigma, m).
#'   Default \code{NULL}.  Forwarded to the bridge detector; see
#'   \code{\link{mt_flag_outliers_bridge}} for accepted forms
#'   (\code{NULL}, scalar, vector, column name, or \code{"auto"} for
#'   Movebank quality columns).
#' @param residual_floor Passed to the bridge detector; see
#'   \code{\link{mt_flag_outliers_bridge}}.
#' @param step_floor Passed to the probability detector; see
#'   \code{\link{mt_flag_outliers}}.
#' @param pool_by Optional character vector of length 1 or 2 naming
#'   column(s) in \code{mt_track_data(x)}.  Pool_by has two semantic
#'   roles: a \emph{fit set} (which tracks' events contribute to the
#'   threshold-fitting distribution) and an \emph{operating unit}
#'   (within which pool-added flags are unioned and the post-cascade
#'   sweep iterates).
#'
#'   \itemize{
#'   \item Length 1 (e.g. \code{pool_by = "individual_id"}): the
#'     same column is used for both roles -- pool deployments of one
#'     animal, with the union scoped to that animal.  This is the
#'     legacy single-level behaviour, preserved byte-identically.
#'   \item Length 2 (e.g. \code{pool_by = c("study_id",
#'     "individual_id")}): the first element is the \emph{outer}
#'     column (fit source) and the second is the \emph{inner} column
#'     (operating unit).  Threshold-fitting primitives draw their
#'     reference distribution from the union of events sharing the
#'     outer value; the post-cascade flag union acts within the
#'     inner value.  This lets you, say, fit thresholds from a
#'     population-wide distribution but keep the union scoped to
#'     each animal.
#'   }
#'
#'   The length-2 form requires \emph{strict nesting}: every distinct
#'   \emph{inner} value must map to exactly one \emph{outer} value.
#'   Inputs that violate this (e.g. an individual that appears in
#'   two studies) error with the offending value named.
#'
#'   Length \eqn{> 2} is rejected with a deliberately verbose
#'   message: pool_by has exactly two semantic roles; a deeper
#'   hierarchy (e.g. species/population/individual/tag) would only
#'   earn its keep under hierarchical / partial-pooling threshold
#'   estimation, which the cascade does not perform.  Users with
#'   many nested levels should pick the \emph{pair} of columns that
#'   captures their trust claim (which level's distribution to fit
#'   from) and their operating unit.
#'
#'   The cascade itself remains per-track (preserved byte-identically
#'   when \code{pool_by = NULL}); a post-cascade pool sweep then
#'   runs each pool-aware primitive wrapper (bridge, speed-cap,
#'   detour) once on the original multi-track input with
#'   \code{pool_by} set, and unions the pool-added flags into
#'   \code{is_outlier}.  Pool-added flags are tagged
#'   \code{error_class = "pool"} when the prior class was empty.
#'
#'   \code{NULL} (default) preserves per-track behaviour
#'   byte-identically.  NA values in the named column(s) cause
#'   those tracks to fall back to per-track processing with a
#'   warning.  Errors when supplied with \code{state} if any
#'   \emph{inner} pool group's tracks have entirely disjoint state
#'   vocabularies (different labelling conventions across
#'   deployments) -- outer-group state inconsistency is not
#'   enforced since outer is a fit source, not a union target.
#'
#'   The prob primitive's contribution to pool semantics is
#'   currently integrated (not post-hoc): \code{mt_flag_outliers}
#'   uses the \emph{outer} column to fit one reference distribution
#'   per outer group and injects it into per-track dispatch.  The
#'   inner column has no role in prob's pool path -- users wanting
#'   prob pooling can call \code{\link{mt_flag_outliers}} standalone
#'   with \code{pool_by} and merge.
#'
#'   Heterogeneous error regimes (e.g. mixed GPS and Sigfox fixes
#'   in one track) violate the per-call homogeneity contract of the
#'   primitives and should be split with \code{dplyr} or
#'   \code{move2} filtering before invoking \code{mt_clean_track}.
#'   See \code{vignette("OUTLIER_heterogeneous_error_regimes", package = "move2utils")} for the
#'   pre-split + re-merge pattern.
#' @param bridge_method Optional character override for the bridge
#'   primitive's \code{method}.  \code{NULL} (default) uses the
#'   cascade's empirically-tuned value \code{"combined"}; valid
#'   alternatives \code{"isotropic"} (scalar Rayleigh residual) and
#'   \code{"directional"} (perpendicular component of the
#'   isotropic/directional decomposition, useful for error-morphology
#'   classification).  See \code{?mt_flag_outliers_bridge} for full
#'   method semantics.
#' @param bridge_threshold_type Optional override for the bridge
#'   primitive's \code{threshold_type}.  \code{NULL} (default) uses
#'   the cascade's value \code{"entropy"} (sweep-validated
#'   density-ratio break detector); alternative \code{"gap"}
#'   (broken-stick + tail-decay; more sensitive but can over-flag).
#' @param bridge_iterations Optional override for the bridge
#'   primitive's iterative refinement count.  \code{NULL} (default)
#'   uses \code{3L}; convergence is typically reached in 1--2 passes
#'   on real data so the override is rarely needed.
#' @param prob_threshold_type Optional override for the probability
#'   primitive's \code{threshold_type}.  \code{NULL} (default) uses
#'   the cascade's empirical default \code{"gap"}; alternatives
#'   \code{"entropy"}, \code{"significance"}, \code{"percentile"} per
#'   \code{?mt_flag_outliers}.
#' @param detour_threshold_type Optional override for the detour
#'   primitive's \code{threshold_type}.  \code{NULL} (default) uses
#'   \code{"fixed"} (the cascade's permissive value gated by the
#'   conjunction rule; threshold = \code{detour_threshold});
#'   alternative \code{"auto"} adapts to each track's
#'   \eqn{-\log(\text{ratio})} distribution.  Note: when
#'   \code{pool_by} is set, the post-cascade pool sweep always uses
#'   \code{threshold_type = "auto"} regardless of this argument so
#'   pool detour has effect; this argument controls only the per-track
#'   cascade detour call.
#' @param entropy_threshold Numeric in (0, 1) or \code{NULL}.
#'   Density-ratio threshold used wherever the cascade runs an
#'   entropy-valley detector (bridge primitive when
#'   \code{bridge_threshold_type = "entropy"}, prob primitive when
#'   \code{prob_threshold_type = "entropy"}, speed-cap auto path's
#'   entropy arm).  \code{NULL} (default) defers to
#'   \code{.entropy_threshold_lower}'s leaf formal -- the package-wide
#'   single source of truth (0.3, sweep-validated 2026-05-06).
#'   Plausible range 0.3--0.7.
#' @param gap_threshold Positive numeric or \code{NULL}.  Break-size
#'   multiplier used wherever the cascade runs a broken-stick gap
#'   detector (bridge / prob / speed-cap auto path's gap fallback /
#'   the class-aware persistence post-filter).  \code{NULL} (default)
#'   defers to \code{.gap_threshold_lower}'s leaf formal (3,
#'   "3-sigma" convention).  Plausible range 2--5.
#' @param persistence_filter_threshold Positive numeric or \code{NULL}.
#'   Per-scale gap-threshold passed to the class-aware persistence
#'   post-filter (active only when
#'   \code{persistence_filter = "class_aware"}).  \code{NULL}
#'   (default) defers to \code{.gap_threshold_lower}'s leaf formal.
#'   Conceptually a separate knob from \code{gap_threshold} because it
#'   operates on a per-scale persistence statistic, not the cascade's
#'   primary detector outputs.
#' @param plot Logical.  Diagnostic map on return.  Default TRUE.
#' @param remove Logical.  If \code{TRUE} (default), return only the
#'   non-flagged rows -- the cleaned track, ready for downstream
#'   analysis.  Set to \code{FALSE} to keep all rows with the flag
#'   columns (\code{is_outlier}, \code{flagged_by_bridge},
#'   \code{flagged_by_prob}, \code{flagged_by_speed},
#'   \code{flag_iteration}, \code{block_id}) attached for inspection
#'   of what was flagged and why.
#' @param silent Logical.  If \code{FALSE} (default) the function
#'   prints a brief running narration of its inner workings: per-
#'   iteration flag counts, the block-expansion gate's decision and
#'   reason, and a final summary line.  Set to \code{TRUE} to
#'   suppress all messages (the same effect as wrapping the call in
#'   \code{suppressMessages()}).  Errors and warnings are always
#'   shown.
#' @param compact Logical.  Used together with \code{silent = FALSE}
#'   to control narration verbosity.  Default \code{FALSE} (the full
#'   per-iteration narration).  Set \code{TRUE} for a one-line-
#'   per-individual summary instead of the full per-iteration
#'   trace -- the right choice when running on a multi-individual
#'   study where the per-iteration output would otherwise scale
#'   linearly with the cohort size.  Ignored when
#'   \code{silent = TRUE}.
#'
#' @section Primitive-knob overrides:
#' The cascade's primitive calls have empirically-tuned defaults
#' that are exposed for fine-tuning at the orchestrator level
#' (\code{bridge_method}, \code{bridge_threshold_type},
#' \code{bridge_iterations}, \code{prob_threshold_type},
#' \code{detour_threshold_type}).  Each \code{NULL}-defaulted
#' argument preserves the cascade's hardcoded value byte-identically;
#' a non-\code{NULL} value forwards to the primitive's
#' \code{.fn_core} call in the iteration loop.  Use these when a
#' paper-replication or diagnostic workflow needs a non-default
#' value without rebuilding the cascade by hand.
#'
#' Two primitive knobs are deliberately NOT exposed because they
#' carry architectural rather than methodological meaning, and
#' overriding them would break the cascade's documented design:
#'
#' \itemize{
#'   \item \strong{detour \code{min_leg}} is hardcoded at \code{0}.
#'     Detour's leg gate (\code{min_leg > 0}) is a standalone-use
#'     gating mechanism that prevents the detour ratio from firing
#'     on small-displacement noise wiggles.  Inside the cascade,
#'     the conjunction rule plays the same gating role: detour
#'     contributes to flags only when it agrees with another
#'     detector under the class-aware rule.  Exposing
#'     \code{min_leg} on the cascade would create two competing
#'     gates running in parallel, with documented failure modes
#'     where a fix is gated out by \code{min_leg} despite
#'     satisfying the conjunction.  Users who specifically want
#'     standalone leg-gated detour should call
#'     \code{\link{mt_flag_outliers_detour}} directly.
#'   \item \strong{speed_cap \code{threshold_type}} is hardcoded
#'     at \code{"auto"}.  The cascade decouples two roles for
#'     speed-based detection: (i) the conjunction's speed flag,
#'     which fires when an implied step speed is anomalous
#'     \emph{relative to the local distribution} (inside a rest
#'     segment, 14 m/s is anomalous even though it sits below the
#'     gull's 36 m/s physiological cap); and (ii) the block-
#'     expansion cap, which uses the user-supplied \code{v_max} (or
#'     \code{(mass, mode)} allometric estimate) as an \emph{absolute}
#'     physiological cap for the boundary-edge component graph.
#'     Auto-threshold for (i) catches state-anomalous fixes that an
#'     absolute cap misses at sparse sampling; absolute \code{v_max}
#'     for (ii) anchors block expansion against true physiological
#'     impossibility.  Collapsing (i) and (ii) into a single
#'     user-overridable threshold reintroduces the over-flag failure
#'     mode documented in the 2026-04-29 stratified Movebank audit
#'     (homing pigeon racing flight over-flagged as anomaly relative
#'     to the resting baseline).  Users who specifically want the
#'     conjunction's speed flag to use the absolute cap should call
#'     \code{\link{mt_flag_speed_cap}(threshold_type = "hard",
#'     v_max = ...)} alongside \code{\link{mt_clean_track}} and
#'     combine the outputs as their workflow requires.
#' }
#'
#' @return The input \code{move2} object with added columns:
#'   \describe{
#'     \item{\code{is_outlier}}{Logical; TRUE where flagged.}
#'     \item{\code{flagged_by_bridge}, \code{flagged_by_prob},
#'           \code{flagged_by_speed}, \code{flagged_by_detour}}{Logical;
#'           per-detector history (useful for diagnosing which signal
#'           caught what).  \code{flagged_by_detour} is identically
#'           \code{FALSE} when \code{use_detour = FALSE}.}
#'     \item{\code{flag_iteration}}{Integer; iteration at which the
#'           fix was flagged (NA otherwise).}
#'     \item{\code{block_id}}{Integer; same value for all fixes in
#'           an expanded-block flag, NA otherwise.}
#'     \item{\code{combined_evidence}}{Numeric; the calibrated combined
#'           evidence behind each fix's decision, attached \emph{only}
#'           under the evidence consensus modes
#'           (\code{"evidence_corroborated"} -- the default --
#'           \code{"weighted_evidence"}, \code{"evidence_or_class"}).
#'           Each fix records its evidence at the moment it was flagged,
#'           or at convergence if kept; higher means more suspicious.
#'           This is the same one-lever sensitivity score
#'           \code{\link{mt_flag_consensus}} returns, exposed so a user
#'           can rank fixes by suspicion or set their own threshold.}
#'     \item{\code{loglr_bridge}, \code{loglr_prob}, \code{loglr_speed},
#'           \code{loglr_detour}}{Numeric; the per-detector signed
#'           log-likelihood-ratios behind \code{combined_evidence},
#'           attached under the same evidence consensus modes. Exposed so
#'           the cascade output can be re-scored by
#'           \code{\link{mt_flag_consensus}} (its \code{evidence_cols}
#'           default) and read by \code{\link{mt_diagnose_flags}}.
#'           \code{NA} under the Boolean rules, which do not compute them.}
#'     \item{\code{error_class}}{Character; categorical interpretation
#'           of why each flagged fix was caught.  \code{NA} where
#'           \code{is_outlier == FALSE}.  Categories:
#'           \itemize{
#'             \item \code{"block"} -- fix is a member of a
#'                   topologically-isolated component (block expansion).
#'                   The mechanistic interpretation is a coherent
#'                   multi-fix error cluster (typical of GPS spoofs,
#'                   timestamp glitches, or systematic tag-data
#'                   contamination) whose interior is locally
#'                   consistent and only the boundary transitions are
#'                   suspect to per-fix scoring.
#'             \item \code{"consensus"} -- two or more of the per-fix
#'                   detectors (bridge, probability, speed) agree on
#'                   the same fix.  Highest confidence in the flag;
#'                   the fix is geometrically, kinematically and/or
#'                   physically anomalous.
#'             \item \code{"speed_cap"} -- only the speed-cap detector
#'                   fired.  The implied step speed exceeds the
#'                   threshold but the fix is geometrically and
#'                   kinematically plausible -- typical of an
#'                   isolated transition that violates a physiological
#'                   cap (e.g. a single jump-and-stay error).
#'             \item \code{"jitter"} -- only the bridge detector
#'                   fired.  The fix is geometrically out of place
#'                   relative to its temporal neighbours but its step
#'                   kinematics are not extreme -- typical of GPS
#'                   multipath, point jitter, or a transcription
#'                   error.
#'             \item \code{"kinematic"} -- only the probability
#'                   detector fired.  The fix's joint speed / turn /
#'                   auto-difference signature conflicts with the
#'                   animal's normal behaviour -- typical of a
#'                   subtle behavioural-state confound or an error
#'                   that is too small to disturb the bridge or
#'                   speed-cap detectors.
#'             \item \code{"detour"} -- only the detour detector
#'                   fired.  The fix's path-vs-displacement ratio is
#'                   pathological but no kinematic detector agrees --
#'                   typical of out-and-back GPS spikes at sparse
#'                   sampling whose implied step speed stays below
#'                   physiological caps.  Only seen when
#'                   \code{use_detour = TRUE}.
#'           }}
#'   }
#'
#' @examples
#' \dontrun{
#' library(move2)
#' x <- movebank_download_study(study_id = 123, ...)
#' x <- mt_filter_gps_quality(x)
#' x <- move2::mt_filter_unique(x, "first")
#' x <- dplyr::arrange(x, mt_time(x))
#'
#' ## Default: return the cleaned track, ready for downstream analysis
#' x_clean <- mt_clean_track(x, v_max = 50)
#'
#' ## Inspection mode: keep all rows, see which were flagged and why
#' x_with_flags <- mt_clean_track(x, v_max = 50, remove = FALSE)
#' table(x_with_flags$is_outlier)
#' head(x_with_flags[x_with_flags$is_outlier, ])
#' }
#'
#' @seealso The four primitives this function composes:
#'   \code{\link{mt_flag_outliers_bridge}},
#'   \code{\link{mt_flag_outliers}},
#'   \code{\link{mt_flag_outliers_detour}},
#'   \code{\link{mt_flag_speed_cap}}.  Diagnostic helpers:
#'   \code{\link{mt_suggest_speed_cap}},
#'   \code{\link{v_phys_estimate}}.  State-aware bridge
#'   primitives (standalone, not currently wired into the
#'   cascade): \code{\link{mt_flag_outliers_bridge}}'s leverage-
#'   immune isotropic/directional decomposition, plus
#'   \code{\link{mt_flag_outliers_dbgb}} and
#'   \code{\link{mt_flag_outliers_dbbmm}} for the variance-
#'   estimating variants.  Alternative strategies for advanced
#'   users (different voting schemes on the probability surface
#'   or across temporal resolutions):
#'   \code{\link{mt_sequential_outliers}},
#'   \code{\link{mt_combined_outliers}},
#'   \code{\link{mt_persistence_score}} for multi-scale persistence
#'   annotation on cascade output (use the score with the
#'   \code{error_class} column for class-aware FP filtering on
#'   \code{state_anomaly} and \code{consensus} flags).
#'
#' @importFrom move2 mt_time mt_track_id mt_distance mt_time_lags
#' @importFrom sf st_coordinates st_is_longlat st_transform
#' @importFrom graphics par plot points legend
#' @importFrom grDevices adjustcolor
#' @importFrom stats setNames
#' @export
mt_clean_track <- function(x,
                            v_max              = NULL,
                            mass               = NULL,
                            mode               = NULL,
                            state              = NULL,
                            consensus          = c("evidence_corroborated",
                                                   "class_aware", "strict",
                                                   "majority", "speed_trusted",
                                                   "weighted_evidence",
                                                   "evidence_or_class",
                                                   "any", "custom"),
                            consensus_custom   = NULL,
                            transition_buffer  = 1L,
                            use_detour         = TRUE,
                            detour_k           = 1L,
                            detour_threshold   = 8,
                            iterations         = "until_clean",
                            max_iterations     = 100L,
                            max_flag_fraction  = 0.2,
                            expand_blocks      = TRUE,
                            pre_peel_aux       = c("none", "primitives"),
                            persistence_filter = c("none", "class_aware"),
                            location_error          = NULL,
                            residual_floor     = 0,
                            step_floor         = 0,
                            pool_by            = NULL,
                            ## Item G (R7) primitive-knob overrides; NULL =
                            ## cascade's empirically-tuned hardcoded defaults
                            ## (byte-identical to no-override).  See
                            ## "Primitive-knob overrides" docstring section.
                            bridge_method         = NULL,
                            bridge_threshold_type = NULL,
                            bridge_iterations     = NULL,
                            prob_threshold_type   = NULL,
                            detour_threshold_type = NULL,
                            ## Package-wide threshold values; NULL =
                            ## leaf-formal defaults (single source of
                            ## truth).  Added 2026-05-25 to close the
                            ## propagation gap; see
                            ## audits/2026-05-25-parameter-propagation/
                            ## findings.md §1.8 and §3.3.
                            entropy_threshold            = NULL,
                            gap_threshold                = NULL,
                            persistence_filter_threshold = NULL,
                            plot               = TRUE,
                            remove             = TRUE,
                            silent             = FALSE,
                            compact            = FALSE) {

  ## ---- input validation ------------------------------------------
  ## verbosity helpers: `say` is summary-level (suppressed only by
  ## silent); `say_iter` is full per-iteration narration (suppressed
  ## by silent OR compact).
  say      <- function(...) if (!silent) message(...)
  say_iter <- function(...) if (!silent && !compact) message(...)
  if (!inherits(x, "move2")) {
    rlang::abort("`x` must be a move2 object.",
                 class = "move2utils_input_not_move2")
  }
  .reject_empty_geometry(x, "mt_clean_track()")
  ## Capture the caller's column set up front.  When `remove = TRUE` the
  ## returned object is a *shrunk* track, so the annotation columns the
  ## cascade adds (is_outlier, flagged_by_*, loglr_*, combined_evidence,
  ## block_id, error_class, ...) no longer align to the original row
  ## indices -- indexing them by a pre-removal position silently reads a
  ## shifted fix.  We therefore strip them on the remove = TRUE path and
  ## return only the caller's columns minus the flagged rows.  Flags are
  ## available intact via `remove = FALSE`.  Captured before any dispatch
  ## branch or the iteration loop touches `x`.
  orig_cols         <- names(x)
  consensus         <- match.arg(consensus)
  if (consensus == "custom" && !is.function(consensus_custom)) {
    rlang::abort(paste0(
      "consensus = \"custom\" requires `consensus_custom` to be a ",
      "function with signature `f(by_bridge, by_prob, by_speed, ",
      "by_detour)` returning a logical vector."),
      class = "move2utils_mt_clean_track_bad_consensus_custom")
  }
  pre_peel_aux      <- match.arg(pre_peel_aux)
  persistence_filter <- match.arg(persistence_filter)
  if (!is.numeric(transition_buffer) || length(transition_buffer) != 1L ||
      is.na(transition_buffer) || transition_buffer < 0 ||
      transition_buffer != as.integer(transition_buffer)) {
    rlang::abort("`transition_buffer` must be a non-negative integer.",
                 class = "move2utils_mt_clean_track_bad_transition_buffer")
  }
  transition_buffer <- as.integer(transition_buffer)
  if (!is.null(v_max) && (!is.numeric(v_max) || length(v_max) != 1 ||
                           is.na(v_max) || v_max <= 0)) {
    rlang::abort("`v_max` must be NULL or a positive scalar.",
                 class = "move2utils_mt_clean_track_bad_v_max")
  }
  ## (mass, mode) is the allometric alternative to v_max.  Both can be
  ## scalar (apply uniformly) or named per-track (one entry per id).
  if (xor(is.null(mass), is.null(mode))) {
    rlang::abort("Provide both `mass` and `mode`, or neither.",
                 class = "move2utils_mt_clean_track_mass_mode_xor")
  }
  if (!is.null(mass) && !is.null(v_max)) {
    rlang::abort("Supply either `v_max` or `(mass, mode)`, not both.",
                 class = "move2utils_mt_clean_track_v_max_with_mass_mode")
  }
  if (!is.null(pool_by)) {
    ## `.resolve_pool_groups` validates shape (length 1 or 2),
    ## column existence, distinctness (for length-2), and the strict-
    ## nesting requirement.  Bad inputs error early.
    pool_maps_check <- .resolve_pool_groups(x, pool_by, silent = TRUE)
    ## state x pool_by vocab check: only meaningful when state is also
    ## supplied AND the inner pool unit (operating scope) spans more
    ## than one track per group.  Per the design (2026-05-12 +
    ## nested-pool-by extension): error if any pair of tracks in an
    ## inner pool group has entirely disjoint state vocabularies.
    ## Outer-group state inconsistency is not enforced -- outer is a
    ## fit source, not a union target.
    if (!is.null(state)) {
      state_resolved <- .resolve_state(state, x)
      if (!is.null(state_resolved)) {
        ids_event <- as.character(move2::mt_track_id(x))
        state_by_track <- split(state_resolved, ids_event)
        .validate_state_pool_vocab(state_by_track, pool_maps_check$inner)
      }
    }
  }

  ## ---- Item G / R7: resolve primitive-knob overrides ----
  ## NULL preserves the cascade's empirically-tuned hardcoded default
  ## (byte-identical to pre-Item-G behaviour).  A non-NULL value is
  ## validated against the primitive's accepted choices and forwarded
  ## to its `.fn_core` call in the iteration loop.
  ##
  ## Architecturally-hardcoded knobs (deliberately NOT exposed): detour
  ## `min_leg` (conjunction-owned gating; exposing creates two
  ## competing gates) and speed_cap `threshold_type` (cascade decouples
  ## relative-speed conjunction flag from absolute v_max block-
  ## expansion cap; exposing collapses the two roles).  See docstring
  ## section "Primitive-knob overrides" for rationale.
  bridge_method_resolved         <- if (is.null(bridge_method)) {
    "combined"
  } else {
    match.arg(bridge_method,
              choices = c("combined", "isotropic", "directional"))
  }
  bridge_threshold_type_resolved <- if (is.null(bridge_threshold_type)) {
    "entropy"
  } else {
    match.arg(bridge_threshold_type, choices = c("entropy", "gap"))
  }
  bridge_iterations_resolved     <- if (is.null(bridge_iterations)) {
    3L
  } else {
    if (!is.numeric(bridge_iterations) ||
        length(bridge_iterations) != 1L ||
        is.na(bridge_iterations) || bridge_iterations < 1 ||
        bridge_iterations != as.integer(bridge_iterations)) {
      rlang::abort("`bridge_iterations` must be NULL or a positive integer.",
                   class = "move2utils_mt_clean_track_bad_bridge_iterations")
    }
    as.integer(bridge_iterations)
  }
  prob_threshold_type_resolved   <- if (is.null(prob_threshold_type)) {
    "gap"
  } else {
    match.arg(prob_threshold_type,
              choices = c("gap", "entropy", "significance", "percentile"))
  }
  detour_threshold_type_resolved <- if (is.null(detour_threshold_type)) {
    "fixed"
  } else {
    match.arg(detour_threshold_type, choices = c("fixed", "auto"))
  }
  ## Per-call threshold value forwarded to the bridge primitive.  The
  ## user can override via `entropy_threshold` / `gap_threshold`
  ## (whichever matches the resolved threshold type); NULL forwards to
  ## the leaf formal default (sweep-friendly).  Pre-2026-05-25 the
  ## cascade resolved `0.3` / `3` here and forwarded the literal to
  ## .bridge_fn_core, shadowing both the leaf and any user-supplied
  ## override.  See audits/2026-05-25-parameter-propagation/findings.md
  ## §1.8.
  bridge_threshold_value <- switch(bridge_threshold_type_resolved,
                                    entropy = entropy_threshold,
                                    gap     = gap_threshold)

  ## ---- upfront recommendation when no cap supplied ---------------
  ## The auto-cap path is reliable on tracks with discrete GPS errors
  ## in single-mode behaviour but has empirically documented failure
  ## modes on (a) block-shaped contamination (Phase 2 K02 spoof
  ## benchmark recovered 0/175 truth at tier 0) and (b) multi-state
  ## behaviour with no real outliers (Phase 3 audit subset over-
  ## flagged 5892 multi-state fixes).  We emit an upfront warning so
  ## the user knows what trade-off they are making by relying on
  ## auto-cap, NOT a silent default.  See:
  ##   benchmarks/2026-XX-XX-vs-competitors/PHASE2_FINDINGS.md
  ##   benchmarks/2026-XX-XX-vs-competitors/PHASE3_FINDINGS.md
  ##   audits/2026-04-29-stratified-movebank/findings.md
  ## The auto-cap advice is suppressed when this function is called
  ## recursively from .dispatch_by_state() -- the user already engaged
  ## state-conditional mode, so the per-segment repetition is noise.
  ## The option is set+restored inside .dispatch_by_state().
  suppress_first_run_advice <-
    isTRUE(getOption("move2utils.suppress_first_run_advice", FALSE))
  if (is.null(v_max) && is.null(mass) && is.null(mode) &&
        !suppress_first_run_advice) {
    say(paste0(
      "No physiological speed cap supplied -- running with a ",
      "data-driven cap chosen from your track.  This works well for ",
      "most cases.  If your animal has multiple behavioural states ",
      "(e.g. perched and flying) or you expect sustained-spoof ",
      "errors, supplying `v_max =` (a published top speed in m/s) ",
      "or `(mass = ..., mode = ...)` for the allometric estimate ",
      "gives sharper results.  See `?v_phys_estimate` for the ",
      "allometric helper; `?mt_clean_track` documents the failure ",
      "modes of the auto-cap in detail."))
  }

  v_max_per_track <- NULL  # named numeric, when set the dispatch picks per id
  if (!is.null(mass)) {
    resolved <- .resolve_mass_mode(mass, mode, x, silent = silent)
    v_max_per_track <- resolved$v_max_named
    if (length(unique(v_max_per_track)) == 1L) {
      ## scalar case -- collapse for the single-track path below
      v_max <- as.numeric(v_max_per_track[[1]])
    }
  }
  ## Accept the string alias "until_clean" in addition to Inf.
  if (is.character(iterations) && length(iterations) == 1L &&
      iterations == "until_clean") {
    iterations <- Inf
  }
  if (!is.numeric(iterations) || length(iterations) != 1 ||
      (!is.infinite(iterations) && iterations <= 0)) {
    rlang::abort(
      "`iterations` must be a positive integer, Inf, or \"until_clean\".",
      class = "move2utils_mt_clean_track_bad_iterations")
  }

  effective_max_iter <- if (is.infinite(iterations)) {
    as.integer(max_iterations)
  } else {
    as.integer(min(iterations, max_iterations))
  }

  ## ---- state-conditional dispatch --------------------------------
  ## When `state` is supplied, partition each track into contiguous
  ## runs of constant state and recurse with state = NULL on every
  ## segment.  Each segment then sees a single kinematic distribution,
  ## so the per-fix threshold detectors are not forced to choose
  ## between cutting-between-modes (over-flag) or cutting-beyond-modes
  ## (under-flag) on multi-state tracks.
  state_vec <- .resolve_state(state, x)
  if (!is.null(state_vec)) {
    return(.dispatch_by_state(
      x = x, state_vec = state_vec,
      v_max_per_track = v_max_per_track, v_max = v_max,
      consensus = consensus,
      consensus_custom = consensus_custom,
      transition_buffer = transition_buffer,
      use_detour = use_detour, detour_k = detour_k,
      detour_threshold = detour_threshold,
      iterations = iterations,
      max_iterations = max_iterations,
      max_flag_fraction = max_flag_fraction,
      expand_blocks = expand_blocks,
      pre_peel_aux = pre_peel_aux,
      persistence_filter = persistence_filter,
      location_error = location_error,
      residual_floor = residual_floor, step_floor = step_floor,
      pool_by = pool_by,
      bridge_method         = bridge_method_resolved,
      bridge_threshold_type = bridge_threshold_type_resolved,
      bridge_iterations     = bridge_iterations_resolved,
      prob_threshold_type   = prob_threshold_type_resolved,
      detour_threshold_type = detour_threshold_type_resolved,
      entropy_threshold            = entropy_threshold,
      gap_threshold                = gap_threshold,
      persistence_filter_threshold = persistence_filter_threshold,
      plot = plot, remove = remove,
      silent = silent, compact = compact))
  }

  ## ---- multi-individual dispatch ---------------------------------
  ids <- move2::mt_track_id(x)
  unique_ids <- unique(ids)
  if (length(unique_ids) > 1L) {
    say("Processing ", length(unique_ids),
                          " individuals separately...")
    ## Suppress the auto-cap newcomer advice on the per-individual
    ## recursive calls -- the outer call (this one) already emitted it
    ## if applicable, so repeating it N times is noise.  Restored on
    ## exit; same mechanism as .dispatch_by_state().
    old_advice_opt <- getOption("move2utils.suppress_first_run_advice", FALSE)
    options(move2utils.suppress_first_run_advice = TRUE)
    on.exit(options(move2utils.suppress_first_run_advice = old_advice_opt),
            add = TRUE)
    results <- lapply(unique_ids, function(id) {
      xi <- x[ids == id, ]
      say(sprintf("--- %s (%d locations) ---",
                                     id, nrow(xi)))
      ## Slice location_error if it was passed as a per-fix numeric vector;
      ## scalar / column-name / "auto" / NULL pass through unchanged.
      oei <- if (is.numeric(location_error) && length(location_error) > 1L) {
        location_error[ids == id]
      } else location_error
      ## Pick this track's v_max if a per-track allometric resolver
      ## produced one; otherwise inherit the scalar v_max above.
      vmi <- if (!is.null(v_max_per_track)) {
        as.numeric(v_max_per_track[[as.character(id)]])
      } else v_max
      mt_clean_track(xi,
                     v_max              = vmi,
                     mass               = NULL,    # already resolved into v_max
                     mode               = NULL,
                     consensus          = consensus,
                     consensus_custom   = consensus_custom,
                     transition_buffer  = transition_buffer,
                     use_detour         = use_detour,
                     detour_k           = detour_k,
                     detour_threshold   = detour_threshold,
                     iterations         = iterations,
                     max_iterations     = max_iterations,
                     max_flag_fraction  = max_flag_fraction,
                     expand_blocks      = expand_blocks,
                     pre_peel_aux       = pre_peel_aux,
                     persistence_filter = persistence_filter,
                     location_error          = oei,
                     residual_floor     = residual_floor,
                     step_floor         = step_floor,
                     pool_by            = NULL,   # outer-only; never recurse with pool_by
                     ## Item G overrides propagate through per-track recursion
                     ## (each child sees the user's resolved values).
                     bridge_method         = bridge_method_resolved,
                     bridge_threshold_type = bridge_threshold_type_resolved,
                     bridge_iterations     = bridge_iterations_resolved,
                     prob_threshold_type   = prob_threshold_type_resolved,
                     detour_threshold_type = detour_threshold_type_resolved,
                     plot               = FALSE,
                     remove             = FALSE,
                     silent             = silent,
                     compact            = compact)
    })
    out <- do.call(rbind, results)
    ## ---- orchestrator pool sweep (multi-individual path) ----
    if (!is.null(pool_by)) {
      out <- .apply_orchestrator_pool_sweep(
                out, x, pool_by,
                location_error = location_error,
                residual_floor = residual_floor,
                detour_k       = detour_k,
                entropy_threshold = entropy_threshold,
                gap_threshold     = gap_threshold,
                silent         = silent)
    }
    if (plot)   .plot_clean_track(out)
    if (remove) out <- .strip_cascade_cols(out[!out$is_outlier, ], orig_cols)
    return(out)
  }

  ## ---- pre-project lon/lat once for the iteration loop ----------
  ## All three detectors run inside an iteration that has historically
  ## paid the cost of a fresh st_transform every pass (the bridge
  ## primitive auto-projects internally on lon/lat input).  Project
  ## once at the top of mt_clean_track to AEQD; the bridge's
  ## auto-detect then sees `st_is_longlat == FALSE` and skips its
  ## per-iteration projection, while .step_lengths_fast() and
  ## .turn_angles_fast() handle projected input directly.  We attach
  ## the resulting flag columns to the original (lon/lat) object at
  ## the end so the caller's CRS is preserved.
  orig_x <- x  # always preserved for back-attach of flag columns at end
  ## Canonicalise to a per-track local AEQD regardless of the input CRS, so
  ## all geometry downstream is computed in metres and the cleaning is
  ## invariant to whatever CRS the caller supplied (longlat, UTM, another
  ## AEQD).  Flag columns are attached back onto orig_x at the end, so the
  ## caller's original CRS is preserved in the output.  (Previously only
  ## longlat input was projected; pre-projected input was used as-is, which
  ## leaked the caller's CRS into the cleaning result -- a bug.)
  x <- .to_canonical_aeqd(x)

  ## ---- resolve location_error + strip user-attached metadata cols
  ## (Phase A of memory-architecture refactor; cf.
  ## CASCADE_AUDIT_2026-05-11.md Section 8.4)
  ##
  ## Movebank-downloaded tracks routinely carry 20+ user-attached
  ## metadata columns (gps_dop, satellite_count, height_above_ellipsoid,
  ## sensor_type, ...).  None of these are read inside the iteration
  ## loop, but each per-iteration `x[active_idx, ]` slice copies them
  ## anyway.  Empirical: 22-39% size reduction on Movebank-typical
  ## tracks just by stripping to the essential columns (time +
  ## track_id + geometry).
  ##
  ## Three care points:
  ##  1. Resolve location_error to a numeric vector NOW (while the
  ##     quality columns are still attached, in case the user passed
  ##     `"auto"` or a column name).  Subsequent primitive calls then
  ##     receive a per-fix numeric -- no further column lookups needed.
  ##  2. Always preserve `orig_x` (regardless of orig_was_longlat) so
  ##     the user gets back their original move2 with just the flag
  ##     columns attached.  Stripped metadata are never lost from the
  ##     caller's perspective.
  ##  3. The lite `x` carries only the columns the cascade itself
  ##     reads (time + track_id + geometry).  Per-iteration slicing
  ##     of this lite object is the structural memory win.
  location_error <- .resolve_location_error(location_error, x, nrow(x))
  keep_cols <- intersect(c(attr(x, "time_column"), attr(x, "track_id_column")),
                          names(x))
  sf_col <- attr(x, "sf_column"); if (is.null(sf_col)) sf_col <- "geometry"
  if (length(keep_cols) && !setequal(keep_cols, setdiff(names(x), sf_col))) {
    x <- x[, keep_cols, drop = FALSE]
  }

  n <- nrow(x)

  ## ---- Phase B+C: cache iteration-invariant inputs ----
  ## Extract coords + time-seconds once.  Each iteration's primitive
  ## calls then consume (cc_all, t_all_s, active_idx, ...) directly via
  ## the `.<primitive>_fn_core` entry points -- no per-iter sf-class
  ## slicing, no per-iter st_coordinates / mt_time / mt_distance /
  ## mt_time_lags / mt_track_id re-extraction.  was_longlat is FALSE
  ## here: when Phase A's AEQD projection ran (lon/lat input), `x` is
  ## projected; when the user already supplied a projected CRS, the
  ## planar math applies directly.
  cc_all  <- if (n > 0L) sf::st_coordinates(x) else matrix(numeric(0), 0L, 2L)
  t_all_s <- if (n > 0L) as.numeric(move2::mt_time(x), units = "secs") else numeric(0)
  was_longlat_for_fn_core <- isTRUE(sf::st_is_longlat(x))  # FALSE after Phase A
  if (n < 10) {
    rlang::warn("Too few locations for mt_clean_track (<10). Returning without flags.",
                class = "move2utils_mt_clean_track_too_few_locations")
    orig_x$is_outlier         <- rep(FALSE, n)
    orig_x$flagged_by_bridge  <- rep(FALSE, n)
    orig_x$flagged_by_prob    <- rep(FALSE, n)
    orig_x$flagged_by_speed   <- rep(FALSE, n)
    orig_x$flagged_by_detour  <- rep(FALSE, n)
    orig_x$flag_iteration     <- rep(NA_integer_, n)
    orig_x$block_id           <- rep(NA_integer_, n)
    orig_x$error_class        <- rep(NA_character_, n)
    return(orig_x)
  }

  ## ---- accumulators ----------------------------------------------
  is_outlier        <- rep(FALSE, n)
  flagged_by_bridge <- rep(FALSE, n)
  flagged_by_prob   <- rep(FALSE, n)
  flagged_by_speed  <- rep(FALSE, n)
  flagged_by_detour <- rep(FALSE, n)
  flag_iteration    <- rep(NA_integer_, n)
  block_id          <- rep(NA_integer_, n)
  ## combined calibrated evidence per fix (evidence consensus modes only):
  ## each fix records its evidence from the last iteration it was alive --
  ## i.e. the evidence at the moment it was flagged, or at convergence if
  ## kept.  Exposed so users have the same one-lever sensitivity score the
  ## standalone mt_flag_consensus() returns.
  combined_evidence <- rep(NA_real_, n)
  ## Per-detector signed log-LRs (same lifecycle as combined_evidence).
  ## Exposed on the output so it is self-describing and can be re-scored by
  ## mt_flag_consensus() under any evidence mode (the consensus evidence_cols
  ## default).  NA under consensus = "class_aware" (no LR computation there).
  loglr_bridge <- rep(NA_real_, n)
  loglr_prob   <- rep(NA_real_, n)
  loglr_speed  <- rep(NA_real_, n)
  loglr_detour <- rep(NA_real_, n)

  ## v_max: supplied or suggested on first iteration (saved for
  ## block expansion at the end)
  v_max_effective <- v_max

  ## ---- pre-step: iterative speed peel if v_max supplied ----------
  ## When the user provides a physiological cap, do the iterative
  ## peel first so multi-fix coherent clusters (spoofs, timestamp
  ## glitches) are removed through boundary-propagation, not via the
  ## per-fix conjunction loop that cannot see cluster interiors.  The
  ## main loop afterwards runs with v_max = NULL so its internal
  ## speed-cap uses auto threshold on the peeled remainder -- avoids
  ## cascading the hard cap through iteration + block expansion.
  used_peel <- FALSE
  if (!is.null(v_max)) {
    ## Optional asymmetric peel: build per-fix auxiliary score from a
    ## single bridge + detour pass; mt_peel_speed then flags only the
    ## higher-scoring endpoint per offending edge instead of both.
    ## Designed for 1-fix-spike-dominated tracks; on cluster-shape
    ## contamination (sustained spoofs) the symmetric default is the
    ## robust choice -- see `?mt_peel_speed` "Cluster-outlier caveat"
    ## and `audits/2026-05-11-cascade-orchestrator/results/03_async_pre_peel.csv`
    ## for the empirical evidence (CPF win, K02 regression).
    aux <- NULL
    if (pre_peel_aux == "primitives") {
      aux <- .build_pre_peel_aux(x,
                                  location_error  = location_error,
                                  residual_floor  = residual_floor,
                                  detour_k        = detour_k,
                                  detour_threshold = detour_threshold)
    }
    peel_res <- .peel_fn_core(
      cc_all, t_all_s, active_idx = seq_len(n), was_longlat_for_fn_core,
      v_max      = v_max,
      aux_scores = aux,
      max_iter   = 1000L,
      silent     = TRUE)
    peel_flags <- peel_res$is_outlier
    is_outlier <- is_outlier | peel_flags
    flagged_by_speed <- flagged_by_speed | peel_flags
    flag_iteration[peel_flags] <- 0L  # peel stage = iteration 0
    used_peel <- TRUE
    say(sprintf(
      "Speed peel (pre-step) at v_max = %g m/s%s: %d fix(es) removed in %d iteration(s).",
      v_max,
      if (pre_peel_aux == "primitives") " [asymmetric, aux=primitives]" else "",
      sum(peel_flags), peel_res$n_peel_iterations))
  }

  ## ---- iteration loop --------------------------------------------
  ## Progress bar is advisory -- the detector converges when no new
  ## flags are produced, typically well before effective_max_iter.
  ## Shown only in interactive sessions; suppressed during tests.
  show_pb <- interactive()
  pb <- if (show_pb) {
    utils::txtProgressBar(min = 0L, max = effective_max_iter, style = 3L)
  } else NULL

  converged_by <- "max_iterations"
  for (iter in seq_len(effective_max_iter)) {
    if (!is.null(pb)) utils::setTxtProgressBar(pb, iter)
    active_idx <- which(!is_outlier)
    if (length(active_idx) < 10L) {
      converged_by <- "active_set_too_small"
      break
    }

    ## Phase B+C: call .fn_core entry points directly with cached
    ## (cc_all, t_all_s) + the iteration's active_idx.  No per-iter
    ## sf-class slicing of `x[active_idx, ]`, no per-iter coord/time/
    ## track-id re-extraction inside each primitive.

    ## --- score 1: bridge ---
    ## Method / threshold_type / iterations are user-overridable via
    ## bridge_method / bridge_threshold_type / bridge_iterations (Item G);
    ## NULL preserves the empirically-tuned defaults shown here.
    ## Threshold value: user-overridable via entropy_threshold /
    ## gap_threshold (matched to bridge_threshold_type); NULL forwards
    ## to the leaf formal (single source of truth).  See
    ## audits/2026-05-25-parameter-propagation/findings.md §1.8.
    bridge_res <- .bridge_fn_core(
      cc_all, t_all_s, active_idx,
      method           = bridge_method_resolved,
      threshold_type   = bridge_threshold_type_resolved,
      threshold        = bridge_threshold_value,
      location_error   = if (is.numeric(location_error)) location_error else NULL,
      residual_floor   = residual_floor,
      iterations       = bridge_iterations_resolved,
      dedup_neighbours = TRUE,
      silent           = TRUE)
    bridge_flags_local <- bridge_res$is_outlier

    ## --- score 2: probabilistic ---
    ## threshold_type is user-overridable via prob_threshold_type (Item G);
    ## NULL preserves the "gap" default validated by the round-4 audit.
    ## The numeric threshold value is now user-overridable via
    ## `entropy_threshold` / `gap_threshold` (matched to the resolved
    ## prob_threshold_type); NULL forwards to the leaf formal.  Pre-
    ## 2026-05-25 the cascade hardcoded `threshold = 3` here regardless
    ## of prob_threshold_type, which (a) shadowed the leaf and (b)
    ## silently dropped a non-default user setting.  See
    ## audits/2026-05-25-parameter-propagation/findings.md §1.8.
    ## entropy path inside .prob_fn_core implements an inline KDE-
    ## valley detector (does NOT call .entropy_threshold_lower) so we
    ## must resolve NULL here.  Read the leaf formal so sweep overrides
    ## of the leaf still propagate.  Gap path forwards NULL through to
    ## .gap_threshold_lower's leaf formal (the conditional in
    ## .prob_fn_core at L763 handles NULL).  significance/percentile
    ## have no leaf; the 0.001 default stays inline.
    prob_threshold_value <- switch(prob_threshold_type_resolved,
                                    entropy      = if (is.null(entropy_threshold))
                                                     formals(.entropy_threshold_lower)$threshold
                                                   else entropy_threshold,
                                    gap          = gap_threshold,
                                    significance = 0.001,
                                    percentile   = 0.001)
    prob_res <- .prob_fn_core(
      cc_all, t_all_s, active_idx, was_longlat_for_fn_core,
      threshold      = prob_threshold_value,
      prob_type      = "joint",
      autodiff_alpha = "acf",
      acf_alpha      = TRUE,
      auto_alpha     = FALSE,
      method         = "histogram",
      time_normalize = TRUE,
      threshold_type = prob_threshold_type_resolved,
      step_transform = "none",
      step_floor     = step_floor,
      silent         = TRUE)
    prob_flags_local <- prob_res$is_outlier

    ## --- score 3: speed cap ---
    ## Always use auto threshold here.  The conjunction's speed flag
    ## should be RELATIVE to the local distribution, not absolute --
    ## inside a "rest" state segment, a 14 m/s implied step is
    ## anomalous even though it sits well below the gull's 36 m/s
    ## physiological cap.  At sparse sampling (1-h GPS) the
    ## physiological cap fires far too rarely to be a useful
    ## conjunction confirmer, defeating the rule on geometric spikes
    ## that detour catches.  The user-supplied physiological cap
    ## continues to gate block expansion (below) -- decoupled from
    ## the conjunction's speed flag.
    speed_res <- .speed_cap_fn_core(
      cc_all, t_all_s, active_idx, was_longlat_for_fn_core,
      v_max                 = NULL,
      threshold_type        = "auto",
      threshold             = NULL,
      jitter                = NULL,
      physiological_ceiling = NULL,
      silent                = TRUE)
    speed_flags_local <- speed_res$is_outlier
    ## v_max_effective tracks the cap used by block expansion below:
    ## prefer the user-supplied physiological cap, fall back to the
    ## auto-cap when nothing else is available.
    if (is.null(v_max) && !is.null(speed_res$v_max_used)) {
      v_max_effective <- speed_res$v_max_used
    }
    ## ---- biological-sanity ceiling on the auto-cap (Tier A1) -----
    ## Mirror at wrapper level the warning emitted inside
    ## mt_flag_speed_cap (which is suppressed by the outer
    ## suppressMessages above).  Threshold 55 m/s is the Hirt 2017
    ## upper 95% CI of the maximum biological speed across all masses
    ## and modes (peak flying ~36.5 m/s + parameter CI extending to
    ## ~52.6 m/s).  Above this the data-driven cap cannot be biology;
    ## it is a structural break inside the outlier tail.
    if (iter == 1L && is.null(v_max) && !used_peel &&
        is.finite(v_max_effective) && v_max_effective > 55) {
      say(sprintf(paste0(
        "Auto-cap landed at %.1f m/s -- above the Hirt 2017 95%% ",
        "upper CI of the maximum biological speed (~52.6 m/s).  The ",
        "gap finder is detecting a structural break within the outlier ",
        "tail. Supply `(mass, mode)` or a hard `v_max` for a ",
        "principled physiological cap.  See `?v_phys_estimate`."),
        v_max_effective))
    }

    ## --- score 4: detour ratio (time-insensitive geometric peer to bridge) ---
    ## threshold_type is user-overridable via detour_threshold_type (Item G);
    ## NULL preserves the "fixed" default that the conjunction rule expects
    ## (permissive threshold gated by conjunction).  detour `min_leg` stays
    ## hardcoded at 0 -- gating is owned by the conjunction rule.
    if (use_detour) {
      detour_res <- .detour_fn_core(
        cc_all, active_idx, was_longlat_for_fn_core,
        k              = detour_k,
        threshold      = detour_threshold,
        threshold_type = detour_threshold_type_resolved,
        min_leg        = 0,
        silent         = TRUE)
      detour_flags_local <- detour_res$is_outlier
    } else {
      detour_flags_local <- rep(FALSE, length(active_idx))
    }

    ## --- lift to full-length indexing ---
    bridge_full <- rep(FALSE, n); bridge_full[active_idx] <- bridge_flags_local
    prob_full   <- rep(FALSE, n); prob_full  [active_idx] <- prob_flags_local
    speed_full  <- rep(FALSE, n); speed_full [active_idx] <- speed_flags_local
    detour_full <- rep(FALSE, n); detour_full[active_idx] <- detour_flags_local

    flagged_by_bridge <- flagged_by_bridge | bridge_full
    flagged_by_prob   <- flagged_by_prob   | prob_full
    flagged_by_speed  <- flagged_by_speed  | speed_full
    flagged_by_detour <- flagged_by_detour | detour_full

    ## --- per-iteration flag rule -------------------------------------
    ## Delegate to the consensus core.  `consensus` defaults to
    ## "class_aware"; other modes (strict / majority / speed_trusted /
    ## any / custom) are exposed via the consensus parameter for users
    ## who want a conservative-vs-liberal knob or fully bespoke voting.
    ## See mt_flag_consensus() for the full specification of each mode.
    ## Evidence-based modes share one computation of the per-detector
    ## log-LRs, the commensurable calibrated evidence, and the class-aware
    ## flags; they differ only in how those are turned into a flag.  All
    ## are data-driven (each detector's surprisal is centred on its OWN
    ## self-thresholded boundary; the calibration MAD is per-track),
    ## principled, and unsupervised.  The innovation lives here -- the
    ## flag->outlier classifier -- not in the (fixed, growing) detectors.
    ev_modes <- c("weighted_evidence", "evidence_or_class",
                  "evidence_corroborated")
    if (consensus %in% ev_modes) {
      ll <- function(surpr, flag) .loglr_from_surprisal(surpr, flag)
      br_score <- switch(bridge_method_resolved,
        isotropic   = bridge_res$bridge_eta,
        directional = bridge_res$bridge_eta_perp,
        pmax(bridge_res$bridge_eta, bridge_res$bridge_eta_perp, na.rm = TRUE))
      sp <- log(speed_res$step_speed)
      sp[!is.finite(speed_res$step_speed) | speed_res$step_speed <= 0] <- NA_real_
      dt <- if (use_detour) {
        d <- log(detour_res$detour_ratio)
        d[!is.finite(detour_res$detour_ratio) | detour_res$detour_ratio <= 1] <- NA_real_
        ll(d, detour_flags_local)
      } else rep(NA_real_, length(active_idx))
      ev_local <- cbind(
        bridge = ll(0.5 * br_score^2,          bridge_flags_local),
        prob   = ll(-log(prob_res$joint),       prob_flags_local),
        speed  = ll(sp,                         speed_flags_local),
        detour = dt)
      ## per-detector commensurable calibration C*tanh(loglr/MAD); E = sum.
      Cc <- 4
      cal1 <- function(v) { s <- stats::mad(v[is.finite(v)], na.rm = TRUE)
        if (!is.finite(s) || s == 0) s <- 1; o <- Cc * tanh(v / s); o[is.na(o)] <- 0; o }
      Cm <- vapply(seq_len(ncol(ev_local)), function(j) cal1(ev_local[, j]),
                   numeric(nrow(ev_local)))
      if (is.null(dim(Cm))) Cm <- matrix(Cm, nrow = nrow(ev_local))
      E_local   <- rowSums(Cm)                       # combined calibrated evidence
      combined_evidence[active_idx] <- E_local       # record for output (one-lever score)
      ## record the raw per-detector log-LRs for the output columns
      loglr_bridge[active_idx] <- ev_local[, "bridge"]
      loglr_prob[active_idx]   <- ev_local[, "prob"]
      loglr_speed[active_idx]  <- ev_local[, "speed"]
      loglr_detour[active_idx] <- ev_local[, "detour"]
      pos_local <- rowSums(ev_local > 0, na.rm = TRUE)  # corroboration count (sign)
      ## "overwhelming solo" = a SATURATED *detour* signal, i.e. detour's
      ## calibrated evidence >= 2 MADs beyond its own boundary (the LR is
      ## conspicuous in its own distribution; 2 = the one principled
      ## robust-SD constant).  Restricted to detour because an out-and-back
      ## is high-specificity (hard to produce by real movement), so a lone
      ## detour can be trusted -- whereas bridge/prob over-react to sharp
      ## but legitimate turns and behavioural-state changes, so a lone
      ## saturated bridge/prob must NOT carry; those require corroboration.
      ## Column order of ev_local/Cm is bridge, prob, speed, detour.
      sat_local <- if (use_detour) Cm[, 4L] > Cc * tanh(2) else rep(FALSE, nrow(Cm))
      class_local <- .consensus_logical(bridge_flags_local, prob_flags_local,
                                        speed_flags_local, detour_flags_local,
                                        mode = "class_aware")
      dec <- switch(consensus,
        weighted_evidence     = E_local > 0,
        evidence_or_class     = class_local | (E_local > 0),
        evidence_corroborated = (E_local > 0) & (pos_local >= 2L | sat_local))
      conj <- rep(FALSE, n); conj[active_idx] <- dec
    } else {
      conj <- .consensus_logical(bridge_full, prob_full, speed_full, detour_full,
                                  mode = consensus, custom = consensus_custom)
    }
    newly_flagged <- conj & !is_outlier
    n_new <- sum(newly_flagged)

    flag_iteration[newly_flagged] <- iter
    is_outlier <- is_outlier | conj

    say_iter(sprintf(
      "Iter %d: bridge=%d prob=%d speed=%d detour=%d (v_max=%s) | conjunction=%d | new=%d cumulative=%d",
      iter,
      sum(bridge_flags_local), sum(prob_flags_local), sum(speed_flags_local),
      sum(detour_flags_local),
      if (is.finite(v_max_effective)) sprintf("%.1f", v_max_effective) else "-",
      sum(conj), n_new, sum(is_outlier)))

    ## --- stop criteria ---
    if (n_new == 0L) {
      converged_by <- "no_new_flags"
      break
    }
    if (sum(is_outlier) > max_flag_fraction * n) {
      rlang::warn(
        sprintf("Cumulative flags (%d) exceed max_flag_fraction (%.1f%%) of track. Aborting.",
                sum(is_outlier), 100 * max_flag_fraction),
        class = "move2utils_mt_clean_track_flag_fraction_exceeded")
      converged_by <- "flag_fraction_exceeded"
      break
    }
  }
  if (!is.null(pb)) {
    utils::setTxtProgressBar(pb, effective_max_iter)
    close(pb)
  }

  ## ---- max_iterations soft-cap warning ---------------------------
  ## If we reached the iteration cap without `n_new == 0` or
  ## `flag_fraction_exceeded`, the cascade was still adding flags
  ## when it ran out of iterations.  Surface this so the user can
  ## opt up `max_iterations` or move to state-conditional dispatch.
  ## Suppressed when explicitly capped via `iterations = N` (the
  ## user asked for a specific iteration count, so hitting it is
  ## not a soft-cap event).
  if (converged_by == "max_iterations" && is.infinite(iterations)) {
    rlang::warn(
      sprintf(paste0(
        "mt_clean_track reached max_iterations = %d without convergence ",
        "(n_new still > 0 in the last iteration; cumulative flags below ",
        "max_flag_fraction).  The cascade is still flagging fixes; ",
        "consider increasing `max_iterations` (e.g. `max_iterations = 200L`) ",
        "or moving to state-conditional dispatch via `state = ...` if ",
        "the track has multiple kinematic regimes."),
        effective_max_iter),
      class = "move2utils_mt_clean_track_max_iterations_reached")
  }

  ## ---- block expansion -------------------------------------------
  ## Skip when pre-peel was used: the peel's iterative nature already
  ## walks into multi-fix clusters via boundary propagation, so
  ## expansion at the same v_max would cascade.  This is a structural
  ## argument, not a tuning choice -- pre-peel and block expansion both
  ## use v_max to define boundary edges, so running them in sequence on
  ## the same partition would double-cascade flag propagation.
  ## Empirical status (Item H, 2026-05-12, per-deployment probe at
  ## `audits/2026-05-11-cascade-orchestrator/scripts/27c_per_track_lift.R`):
  ## across CPF + cohort + K02 the counterfactual lift adds 0 flags on
  ## every (track, cap-mode) pair.  The exclusion is structurally
  ## defensible; no test-set evidence currently distinguishes it from
  ## the lifted variant.  (Round-4 #C4's +1206-flag finding on Rhino
  ## was a probe artifact: the gate was run on the FULL multi-
  ## deployment object, cutting at the track-id boundary; the real
  ## cascade dispatches per-track at line 553 and never sees that cut.
  ## Dual-deployment-aware pooling is queued as a v0.3 design item via
  ## `pool_by =` -- see FUTURE_IMPROVEMENTS.md.)
  ##
  ## Otherwise: build the component partition that the cut at
  ## v_max_effective would produce, then ask the gate (dip-test +
  ## gap on log-component-sizes) whether the cut is principled or
  ## severs continuous behaviour.  Apply expansion only if the gate
  ## allows.  See .gate_block_expansion for the rationale.
  ## Physiological-plausibility floor on v_max_effective.  The auto-cap
  ## path occasionally lands at sub-physiological values on multi-state
  ## tracks (e.g. WH17: auto-cap at 0.58 m/s -- well inside the rest-
  ## state speed distribution).  When that happens, the partition's
  ## component structure reflects rest periods, not contamination, and
  ## block expansion compounds the auto-cap failure by flagging entire
  ## rest blocks.  The 1 m/s floor sits below any flying / running /
  ## swimming species' sustained motion speed (Hirt 2017 minima: 36 / 17
  ## / 13 m/s respectively) but above any plausible rest speed; an auto-
  ## cap below it is not a physiological cap, so block expansion should
  ## not run.  This is a defence-in-depth -- the cascade already warns
  ## about high auto-caps (> 55 m/s upper-CI); this is the corresponding
  ## lower guard.  See HEURISTICS.md Group 2 entry "block_expansion
  ## v_max floor" + audits/2026-05-11-cascade-orchestrator/scripts/22*.R
  ## for the WH17 auto empirical case (9623 -> 6725 with the guard).
  block_v_max_floor <- 1.0
  block_gate <- NULL
  if (expand_blocks && is.finite(v_max_effective) &&
        v_max_effective > block_v_max_floor && sum(is_outlier) > 0) {
    part       <- .compute_component_partition(x, is_outlier,
                                                v_max = v_max_effective)
    block_gate <- .gate_block_expansion(part,
                                         max_block_fraction = max_flag_fraction)
    say(sprintf("Block-expansion gate: %s -- %s.",
                if (block_gate$allow) "ALLOWED" else "DECLINED",
                block_gate$reason))
    if (block_gate$allow) {
      exp_result <- .apply_block_expansion(part, block_gate)
      new_block_flags <- !is.na(exp_result) & !is_outlier
      n_block_expanded <- sum(new_block_flags)
      if (n_block_expanded > 0L) {
        say(sprintf("Block expansion: %d additional fixes flagged in %d block(s).",
                    n_block_expanded,
                    length(unique(stats::na.omit(exp_result)))))
        is_outlier <- is_outlier | new_block_flags
        flag_iteration[new_block_flags] <- effective_max_iter + 1L
      }
      block_id <- exp_result
    } else if (is.null(v_max)) {
      say("    The auto-cap did not yield a clean block partition; supply a physiological cap via `v_max = ...` or `(mass, mode)` for a firmer connectivity ceiling.")
    }
  }

  ## ---- write columns ---------------------------------------------
  x$is_outlier        <- is_outlier
  x$flagged_by_bridge <- flagged_by_bridge
  x$flagged_by_prob   <- flagged_by_prob
  x$flagged_by_speed  <- flagged_by_speed
  x$flagged_by_detour <- flagged_by_detour
  x$flag_iteration    <- flag_iteration
  x$block_id          <- block_id
  ## expose the combined calibrated evidence under the evidence consensus
  ## modes -- the same one-lever sensitivity score mt_flag_consensus()
  ## returns (NA under the Boolean rules, which do not compute it).
  if (consensus %in% c("weighted_evidence", "evidence_or_class",
                       "evidence_corroborated")) {
    x$combined_evidence <- combined_evidence
    ## also expose the per-detector log-LRs so the cascade output can be
    ## re-scored by mt_flag_consensus() (its evidence_cols default) and read
    ## by mt_diagnose_flags(); previously these were computed but dropped.
    x$loglr_bridge <- loglr_bridge
    x$loglr_prob   <- loglr_prob
    x$loglr_speed  <- loglr_speed
    x$loglr_detour <- loglr_detour
  }
  classified <- .classify_flags(is_outlier, flagged_by_bridge,
                                 flagged_by_prob, flagged_by_speed,
                                 block_id, flagged_by_detour,
                                 flag_iteration = flag_iteration)
  x$error_class       <- classified$error_class
  x$error_classes_all <- classified$error_classes_all
  attr(x, "convergence") <- converged_by
  attr(x, "v_max_used")  <- v_max_effective

  ## ---- optional persistence post-filter --------------------------
  ## Per CASCADE_AUDIT_2026-05-11 R2 + the class-conditional finding
  ## from 2026-05-09: mt_persistence_score gives a +37-39 pp TP-vs-FP
  ## gap on `state_anomaly` + `consensus` cascade classes at p >= 3.
  ## On other classes (geometric_spike, kinematic_confluence,
  ## physiological, block) persistence is anti-informative or neutral.
  ## With `persistence_filter = "class_aware"` we apply the filter
  ## ONLY to those two classes, demoting flagged fixes whose
  ## persistence_count < 3 to is_outlier = FALSE.  Default "none"
  ## preserves current behaviour exactly.
  if (persistence_filter == "class_aware" &&
        any(x$is_outlier, na.rm = TRUE)) {
    target_class <- !is.na(x$error_class) &
                     x$error_class %in% c("state_anomaly", "consensus")
    if (any(target_class)) {
      ## persistence_filter_threshold defaults to NULL -> leaf formal
      ## of .gap_threshold_lower (sweep-friendly).  Pre-2026-05-25 the
      ## cascade hardcoded `threshold = 3` here.  See
      ## audits/2026-05-25-parameter-propagation/findings.md §1.9.
      ps_fn_args <- list(cc_all, active_idx = seq_len(n),
                          candidate_mask = x$is_outlier,
                          scales = c(2L, 4L, 8L), n_breaks = 20L,
                          silent = TRUE)
      if (!is.null(persistence_filter_threshold))
        ps_fn_args$threshold <- persistence_filter_threshold
      ps_res <- do.call(.persistence_fn_core, ps_fn_args)
      pers <- ps_res$persistence_count
      to_drop <- target_class & !is.na(pers) & pers < 3L
      n_drop  <- sum(to_drop)
      if (n_drop > 0L) {
        x$is_outlier[to_drop]     <- FALSE
        x$error_class[to_drop]    <- NA_character_
        x$flag_iteration[to_drop] <- NA_integer_
        ## Recompute is_outlier accumulator + error_classes_all
        is_outlier <- x$is_outlier
        flag_iteration <- x$flag_iteration
        say(sprintf(
          "Persistence filter (class_aware): demoted %d flag(s) in state_anomaly/consensus with persistence_count < 3.",
          n_drop))
      }
    }
  }

  ## Always attach flag columns to the preserved orig_x.  The lite
  ## working `x` may have been stripped of user-attached metadata
  ## (Phase A) and/or projected to AEQD; the user gets back their
  ## original move2 (full columns + original CRS) with just the
  ## per-cascade flag columns appended.
  flag_cols <- c("is_outlier", "flagged_by_bridge", "flagged_by_prob",
                 "flagged_by_speed", "flagged_by_detour",
                 "flag_iteration", "block_id", "error_class",
                 "error_classes_all", "combined_evidence",
                 "loglr_bridge", "loglr_prob", "loglr_speed", "loglr_detour")
  for (col in flag_cols) {
    if (col %in% names(x)) orig_x[[col]] <- x[[col]]
  }
  attr(orig_x, "convergence") <- attr(x, "convergence")
  attr(orig_x, "v_max_used")  <- attr(x, "v_max_used")
  x <- orig_x

  ## ---- primitive-disagreement diagnostic -------------------------
  ## On a clean unimodal track the bridge primitive is silent
  ## (mt_flag_outliers_bridge default = entropy-valley threshold which
  ## returns zero flags when no valley exists).  When bridge fires but
  ## the conjunction drops most of its signal, the structural picture
  ## is "bridge sees geometric anomalies that prob and speed do not
  ## confirm" -- the canonical signature of block contamination
  ## (sustained spoof, deployment confusion) where boundary fixes are
  ## geometrically anomalous but kinematically self-consistent.  We
  ## emit the diagnostic when:
  ##   (a) no physiological cap was supplied (auto-cap path),
  ##   (b) bridge caught fixes that the final flag set did not include
  ##       (n_dropped > 0).
  ##
  ## Structurally-grounded test: warn when bridge caught MORE fixes
  ## than the conjunction kept (i.e. n_dropped > n_final).  Ratio-based
  ## so it scales with track size; the absolute n_dropped >= 3L
  ## threshold the audit used to ship was a heuristic tuned to K02 at
  ## n_dropped = 5, which became obsolete once class_aware flagging
  ## + state-conditional dispatch reduced K02 / auto's n_dropped to 0.
  ## Empirical motivation: `audits/2026-05-11-cascade-orchestrator/
  ## scripts/26_R5_magic_numbers.R` -- ratio > 1 fires on the four
  ## legitimate disagreement cases on the test set (Saline, 99696651,
  ## WH17, E28946) and skips the synthetic-noise floor (CPF_E, CPF_F)
  ## that the absolute threshold caught as false alarms.
  n_final     <- sum(is_outlier, na.rm = TRUE)
  n_bridge    <- sum(flagged_by_bridge, na.rm = TRUE)
  n_dropped   <- max(0L, n_bridge - n_final)
  no_v_max_supplied <- is.null(v_max_per_track) &&
                        (is.null(v_max) || !is.finite(v_max))
  if (no_v_max_supplied && n_dropped > n_final) {
    say(sprintf(paste0(
      "Primitive-disagreement signature: bridge caught %d fixes that ",
      "the conjunction dropped (final flags = %d). Bridge is silent ",
      "on clean unimodal data; non-trivial bridge-only flags suggest ",
      "geometric anomalies prob and speed disagree on, the canonical ",
      "block-contamination signature. Supply `(mass, mode)` or a hard ",
      "`v_max` to mt_clean_track() to engage the speed primitive on ",
      "block boundaries. See `?v_phys_estimate` and ",
      "vignette(\"OUTLIER_2_diagnose_clean_track\")."),
      n_bridge, n_final))
  }

  say(sprintf("=== mt_clean_track: %d flagged (%.3f%% of %d); stopped: %s ===",
              sum(is_outlier), 100 * sum(is_outlier) / n, n, converged_by))
  if (sum(is_outlier) == 0L && !compact) {
    say("    No outliers flagged.  This is normal for clean tracks; for ",
        "tracks with multiple behavioural states (e.g. perched + flying) ",
        "the cascade can come up empty until you switch to state-conditional ",
        "cleaning.  Run `mt_diagnose_clean_track()` on the result to confirm ",
        "the run was healthy and to see the recommended next step.")
  }
  if (!compact) {
    if (remove) {
      say(sprintf("    Returning the cleaned track (%d rows). To inspect what was flagged, re-run with remove = FALSE.",
                  n - sum(is_outlier)))
    } else {
      say("    Returning all rows with flag columns attached. To drop flagged rows, either re-run with remove = TRUE (the default) or subset: x[!x$is_outlier, ].")
    }
  }

  if (plot)   .plot_clean_track(x)
  if (remove) x <- .strip_cascade_cols(x[!x$is_outlier, ], orig_cols)
  x
}


## ---- helpers ----------------------------------------------------

## On the `remove = TRUE` return, drop every column the cascade added
## and hand back only the caller's original columns (geometry is sticky
## in sf, so it survives the intersect regardless).  Keeps the cleaned
## object free of row-misaligned flag columns.  See the `orig_cols`
## capture in mt_clean_track().
.strip_cascade_cols <- function(x, orig_cols) {
  keep <- intersect(names(x), orig_cols)
  x[, keep]
}


## Build a per-fix auxiliary outlier score for asymmetric pre-peel.
## Runs a single non-iterative bridge + detour pass on the raw input,
## then combines their rank-normalised magnitudes.  The result is a
## numeric vector aligned to nrow(x), suitable for passing as
## `aux_scores` to `mt_peel_speed()`.
##
## Higher score = more likely outlier.  Bridge contributes the
## leverage-immune perpendicular residual (`bridge_eta`); detour
## contributes the time-insensitive path/displacement ratio
## (`detour_ratio`).  Each is rank-normalised on its non-zero values
## (NA / non-finite -> 0) and summed.  Tracks with nothing to score
## return all-zero, in which case the asymmetric peel falls back to
## flagging the right (later-in-time) endpoint per edge.
##
## @keywords internal
.build_pre_peel_aux <- function(x, location_error, residual_floor,
                                  detour_k, detour_threshold) {
  br <- suppressMessages(suppressWarnings(
    mt_flag_outliers_bridge(x,
                             location_error = location_error,
                             residual_floor = residual_floor,
                             plot = FALSE)))
  dt <- suppressMessages(suppressWarnings(
    mt_flag_outliers_detour(x,
                              k = detour_k, threshold = detour_threshold,
                              min_leg = 0, plot = FALSE, silent = TRUE)))
  e <- if (!is.null(br$bridge_eta)) br$bridge_eta else rep(0, nrow(x))
  d <- if (!is.null(dt$detour_ratio)) dt$detour_ratio else rep(0, nrow(x))
  e[!is.finite(e)] <- 0
  d[!is.finite(d)] <- 0
  rank_norm <- function(v) {
    out <- numeric(length(v))
    nz  <- which(v > 0)
    if (length(nz) > 1L)
      out[nz] <- rank(v[nz], ties.method = "average") / length(nz)
    out
  }
  rank_norm(e) + rank_norm(d)
}


## Resolve the (mass, mode) allometric specification into a per-track
## named numeric of v_max values.  Accepts:
##   - mass = scalar in kg  -> single v_max applied to every track
##   - mass = named numeric -> per-track; names must be a superset of
##                              the track ids in `x`
## Mode is treated symmetrically (scalar string OR named character).
##
## Two ergonomic conveniences:
##   - if `mass` carries a `units` attribute, auto-convert to kg
##     (catches the very common Movebank `animal_mass [g]` pattern)
##   - if a bare scalar `mass > 100` we warn that it looks like
##     grams; the user almost certainly wants kg
##
## Returns a list:
##   v_max_named  named numeric of v_max in m/s, one per track id
##   summary      character lines for narration (one per track)
##
## @keywords internal
.resolve_mass_mode <- function(mass, mode, x, silent = FALSE) {
  say <- .say(silent)
  unique_ids <- as.character(unique(move2::mt_track_id(x)))

  ## --- units auto-handling ---
  if (inherits(mass, "units")) {
    u <- attr(mass, "units")
    unit_str <- if (!is.null(u)) {
      paste(format(u), collapse = "")
    } else ""
    nm <- names(mass)
    mass <- as.numeric(units::set_units(mass, "kg", mode = "standard"))
    if (!is.null(nm)) names(mass) <- nm
    if (length(unit_str) == 1L && nzchar(unit_str)) {
      say(sprintf("  mass: converting from %s to kg.", unit_str))
    }
  }
  if (!is.numeric(mass)) {
    rlang::abort("`mass` must be numeric (kg) or a units-attributed numeric.",
                 class = "move2utils_mt_clean_track_bad_mass_type")
  }

  ## --- "looks like grams" guardrail (only on bare scalar without units) ---
  if (length(mass) == 1L && is.finite(mass) && mass > 100) {
    rlang::warn(
      sprintf("`mass = %g` is unusually large for kg.  If this is grams, divide by 1000 -- the allometric helper expects kg.",
              mass),
      class = "move2utils_mt_clean_track_mass_looks_like_grams")
  }

  ## --- expand to per-track named numeric ---
  ## A named vector (any length) is treated as per-track; a bare
  ## scalar is broadcast to all tracks.  This lets a user pass
  ## `mass = c("X" = 1.5)` and get a clear "missing entries" error
  ## rather than silently broadcasting 1.5 to every individual.
  if (is.null(names(mass))) {
    if (length(mass) != 1L) {
      rlang::abort(
        "Unnamed per-track `mass` is ambiguous; supply a named numeric whose names are track ids, or a scalar to apply uniformly.",
        class = "move2utils_mt_clean_track_ambiguous_mass_names")
    }
    mass_named <- setNames(rep(mass, length(unique_ids)), unique_ids)
  } else {
    if (any(!nzchar(names(mass)))) {
      rlang::abort(
        "Per-track `mass` must have non-empty names matching track ids.",
        class = "move2utils_mt_clean_track_bad_mass_names")
    }
    miss <- setdiff(unique_ids, names(mass))
    if (length(miss)) {
      rlang::abort(
        sprintf("Per-track `mass` is missing entries for: %s",
                paste(miss, collapse = ", ")),
        class = "move2utils_mt_clean_track_missing_mass_entries")
    }
    mass_named <- mass[unique_ids]
  }

  ## --- mode: scalar OR per-track ---
  valid_modes <- c("flying", "running", "swimming")
  if (is.null(names(mode))) {
    if (length(mode) != 1L) {
      rlang::abort(
        "Unnamed per-track `mode` is ambiguous; supply a named character whose names are track ids, or a scalar to apply uniformly.",
        class = "move2utils_mt_clean_track_ambiguous_mode_names")
    }
    mode <- match.arg(mode, valid_modes)
    mode_named <- setNames(rep(mode, length(unique_ids)), unique_ids)
  } else {
    miss <- setdiff(unique_ids, names(mode))
    if (length(miss)) {
      rlang::abort(
        sprintf("Per-track `mode` is missing entries for: %s",
                paste(miss, collapse = ", ")),
        class = "move2utils_mt_clean_track_missing_mode_entries")
    }
    mode_named <- mode[unique_ids]
    bad <- !mode_named %in% valid_modes
    if (any(bad)) {
      rlang::abort(
        sprintf("Invalid mode value(s): %s.  Allowed: %s",
                paste(unique(mode_named[bad]), collapse = ", "),
                paste(valid_modes, collapse = ", ")),
        class = "move2utils_mt_clean_track_invalid_mode")
    }
  }

  ## --- compute v_max per track (silently; we narrate aggregate) ---
  v_max_named <- vapply(unique_ids, function(id) {
    as.numeric(suppressWarnings(
      v_phys_estimate(mass = mass_named[[id]], mode = mode_named[[id]])))
  }, numeric(1))
  names(v_max_named) <- unique_ids

  if (length(unique(mass_named)) == 1L && length(unique(mode_named)) == 1L) {
    say(sprintf(
      "Allometric v_max from Hirt et al. (2017): %.2f m/s (mass = %g kg, mode = %s); applied to all %d individual(s).",
      v_max_named[[1]], mass_named[[1]], mode_named[[1]],
      length(unique_ids)))
  } else {
    say("Allometric v_max from Hirt et al. (2017), per-track:")
    for (id in unique_ids) {
      say(sprintf("  %-30s mass=%g kg, mode=%-8s -> v_max = %.2f m/s",
                  substr(id, 1, 30), mass_named[[id]],
                  mode_named[[id]], v_max_named[[id]]))
    }
  }

  list(v_max_named = v_max_named,
       mass_named  = mass_named,
       mode_named  = mode_named)
}


## Map per-detector flag combinations to a categorical error class.
## Returns NA where !is_outlier.  Class precedence (high to low):
##   "physiological" (pre-peel, flag_iteration == 0)
##   "block"         (block_id non-NA)
##   "consensus"     (>= 3 of 4 detectors fired)
##   "geometric_spike" (bridge AND detour, < 3 detectors total)
##   "state_anomaly" ((bridge|detour) AND speed, no prob, < 3 total)
##   "kinematic_confluence" ((bridge|detour) AND prob, no speed)
##
## The taxonomy is independent of the consensus mode used to decide
## `is_outlier`.  A fix flagged by `consensus = "majority"` (a vote-
## tallying mode) still gets a class label here, describing the
## agreement structure that caused it to be flagged.
##
## Returns a list with two elements:
##   error_class       : single string per fix (highest-priority class)
##   error_classes_all : comma-separated string of every class that
##                       fired (or "" where !is_outlier).  Useful when
##                       downstream code wants the full agreement
##                       structure rather than just the priority winner.
##
## @keywords internal
.classify_flags <- function(is_outlier, by_bridge, by_prob, by_speed,
                              block_id, by_detour = NULL,
                              flag_iteration = NULL) {
  n <- length(is_outlier)
  out      <- rep(NA_character_, n)
  out_all  <- rep("",            n)

  if (is.null(by_detour))     by_detour     <- rep(FALSE, n)
  if (is.null(flag_iteration)) flag_iteration <- rep(NA_integer_, n)

  if (!any(is_outlier))
    return(list(error_class = out, error_classes_all = out_all))

  in_block <- !is.na(block_id)
  in_peel  <- !is.na(flag_iteration) & flag_iteration == 0L
  votes    <- as.integer(by_bridge) + as.integer(by_prob) +
              as.integer(by_speed)  + as.integer(by_detour)
  geo      <- by_bridge | by_detour

  ## class_aware taxonomy: PRIMARY (highest-priority) class for `error_class`.
  ## Assignment in REVERSE precedence so the highest-priority
  ## overwrites lower:
  fire_kin   <- is_outlier & !in_block & geo & by_prob
  fire_state <- is_outlier & !in_block & geo & by_speed
  fire_geom  <- is_outlier & !in_block & by_bridge & by_detour
  fire_cons  <- is_outlier & !in_block & votes >= 3L
  fire_block <- is_outlier &  in_block
  fire_peel  <- is_outlier &  in_peel
  out[fire_kin]   <- "kinematic_confluence"
  out[fire_state] <- "state_anomaly"
  out[fire_geom]  <- "geometric_spike"
  out[fire_cons]  <- "consensus"
  out[fire_block] <- "block"
  out[fire_peel]  <- "physiological"

  ## Comma-separated full set of classes that fired per fix.  Useful
  ## downstream when a user wants to see the full agreement structure
  ## rather than just the priority winner.  Order in the string:
  ## physiological, block, consensus, geometric_spike, state_anomaly,
  ## kinematic_confluence.
  p1 <- ifelse(fire_peel,  "physiological",        "")
  p2 <- ifelse(fire_block, "block",                "")
  p3 <- ifelse(fire_cons,  "consensus",            "")
  p4 <- ifelse(fire_geom,  "geometric_spike",      "")
  p5 <- ifelse(fire_state, "state_anomaly",        "")
  p6 <- ifelse(fire_kin,   "kinematic_confluence", "")
  combined <- paste(p1, p2, p3, p4, p5, p6, sep = ",")
  combined <- gsub(",+", ",", combined)
  combined <- gsub("^,|,$", "", combined)
  out_all <- combined

  list(error_class = out, error_classes_all = out_all)
}


## Compute the component partition of the kept-fix graph at a given
## connectivity speed cap.  Two consecutive kept fixes share a
## component iff the implied step speed between them is <= v_max.
##
## Returns a list with:
##   comp_kept    integer vector of length length(kept), the component
##                label for each kept fix
##   sizes        integer vector of component sizes
##   kept         integer indices into x of the kept fixes
##   n            nrow(x)
##
## @keywords internal
.compute_component_partition <- function(x, is_outlier, v_max) {
  n <- nrow(x)
  kept <- which(!is_outlier)
  if (length(kept) < 2L) {
    return(list(comp_kept = integer(0),
                sizes     = integer(0),
                kept      = kept,
                n         = n))
  }
  x_kept <- x[kept, ]
  step_m <- .step_lengths_fast(x_kept)
  dt_s   <- as.numeric(move2::mt_time_lags(x_kept, units = "secs"))
  ## Original-timing guard (monotonicity fix).  The connectivity cut is
  ## computed on the KEPT graph, so removing the seam fixes of a coherent
  ## boundary block widens the kept-to-kept time gap and the implied
  ## across-seam speed dilutes below v_max -- the block then re-joins the
  ## main component and block expansion cannot isolate it.  This is the
  ## mechanism by which supplying a physiological cap used to *reduce*
  ## block recovery (the `!used_peel` short-circuit only hid it).  Where
  ## flagged fixes were removed BETWEEN two kept fixes (diff(kept) > 1),
  ## cap the lag at the track's median consecutive sampling interval so a
  ## large displacement across a removed span stays severable.  Genuine
  ## single-step gaps (diff == 1, including real missing-data gaps that
  ## leave no removed fix) keep their true lag and are never over-cut --
  ## an isolated spike's kept neighbours are spatially CLOSE, so even with
  ## the cap their implied speed stays sub-v_max and they correctly
  ## re-join.  See HEURISTICS.md "block-expansion removed-gap dt cap".
  full_dt <- as.numeric(move2::mt_time_lags(x, units = "secs"))
  med_int <- stats::median(full_dt[is.finite(full_dt) & full_dt > 0],
                           na.rm = TRUE)
  removed_between <- c(diff(kept) > 1L, FALSE)   # edge i = kept[i] -> kept[i+1]
  if (is.finite(med_int) && any(removed_between)) {
    dt_s[removed_between] <- pmin(dt_s[removed_between], med_int)
  }
  speed  <- ifelse(is.finite(step_m) & is.finite(dt_s) & dt_s > 0,
                   step_m / dt_s, NA_real_)
  cut_here <- !is.finite(speed) | speed > v_max
  ## comp_kept[j] = component label of kept fix j.  cut_here has one entry
  ## per kept fix; its last entry is the NA-speed trailing edge (no
  ## successor), so trim the cumulative labelling back to length(kept) --
  ## the previous `c(1L, 1L + cumsum(cut_here))` returned length(kept)+1,
  ## leaving a phantom size-1 component that inflated the size table.
  comp_kept <- c(1L, 1L + cumsum(cut_here))[seq_along(kept)]
  list(comp_kept = comp_kept,
       sizes     = as.integer(tabulate(comp_kept)),
       kept      = kept,
       n         = n)
}


## Gate the block-expansion step.
##
## Block expansion uses v_max as a connectivity threshold: it severs
## the kept-fix graph wherever the step speed exceeds v_max, then
## flags the small resulting components as outlier "blocks".  This is
## sound only when v_max is a credible physiological cap that real
## behaviour does not legitimately exceed.  When v_max comes from a
## data-driven step-cap detector on a fast-moving species it can sit
## well within the bird's normal flight range, in which case the cut
## carves up the trajectory rather than isolating discrete error
## clusters.
##
## The gate asks two questions about the component-size distribution
## the cut would produce, both expressed in the same broken-stick /
## gap-aware language used throughout the package:
##
##   1. Dominance.  Does the largest component hold at least
##      \code{1 - max_block_fraction} of the kept fixes?  This is the
##      necessary condition for "one trajectory + a few isolated
##      clusters" to be a coherent reading of the partition.  When the
##      largest component does not dominate, the cap has fragmented
##      the trajectory itself; expansion would carve up real movement.
##
##   2. Upper-tail gap (when the partition has 4+ components).  Is
##      there a clear gap on log-component-size between the largest
##      component and the rest?  Specifically, the gap from the bulk
##      of small components to the largest must exceed the dispersion
##      of the small-component bulk by at least \code{gap_factor}
##      (default 3) -- the same broken-stick / noise-floor reasoning
##      used by \code{.gap_threshold_lower} on per-step speeds, but
##      applied to the upper tail of the size distribution.
##
## When both pass, the gate returns the second-largest component's
## size as the cut: every non-largest component is flagged as a block.
## When either fails, expansion is declined and the orchestrator
## proceeds with the per-fix consensus flags only.
##
## Returns a list:
##   allow         logical, TRUE iff expansion should proceed
##   reason        character, brief explanation suitable for messaging
##   n_components  integer
##   sizes         integer vector of component sizes (sorted desc)
##   size_break    integer, components with size strictly below this
##                 are flagged as blocks (NA when allow = FALSE)
##   largest_frac  numeric, fraction of kept fixes in the largest
##                 component
##   gap_factor    numeric, ratio of the log-gap above the bulk to
##                 the bulk's dispersion (NA when not computed)
##
## @keywords internal
## gap_factor default 3 is HEURISTIC -- reuses the "3-sigma" convention
## the broken-stick gap detector uses on per-step speeds. It is the
## ratio of the log-gap above the bulk to the bulk's robust dispersion
## that the gate requires to declare a clean separation. Plausible
## range: 2--5. Lower values are more permissive (more block-expansion
## firings, more risk of false splits); higher values are more
## conservative (block expansion declines more often, more risk of
## leaving long contaminated trains uncleaned). The benchmark's K02
## tier 0 case had ratio = 0.29 (correctly declined; speed primitive
## hadn't fired on boundaries).
.gate_block_expansion <- function(part, max_block_fraction,
                                    gap_factor = 3) {
  sizes <- part$sizes
  k <- length(sizes)
  sizes_sorted <- sort(sizes, decreasing = TRUE)

  ## trivial: 0 or 1 components => nothing to expand
  if (k < 2L) {
    return(list(allow = FALSE,
                reason = sprintf(
                  "only %d component(s); nothing to expand", k),
                n_components = k, sizes = sizes_sorted,
                size_break = NA_integer_, largest_frac = NA_real_,
                gap_factor = NA_real_))
  }

  ## --- dominance test -------------------------------------------
  largest_frac <- sizes_sorted[1] / sum(sizes)
  if (largest_frac < 1 - max_block_fraction) {
    return(list(allow = FALSE,
                reason = sprintf(
                  "largest component holds only %.1f%% of kept fixes (< %.0f%% required) -- cut severs continuous trajectory",
                  100 * largest_frac, 100 * (1 - max_block_fraction)),
                n_components = k, sizes = sizes_sorted,
                size_break = NA_integer_, largest_frac = largest_frac,
                gap_factor = NA_real_))
  }

  ## With only 2 or 3 components and the largest dominating, the cap
  ## has produced a clean "trajectory + clusters" partition and we
  ## allow without further checks (too few components for a stable
  ## dispersion estimate).
  if (k <= 3L) {
    return(list(allow = TRUE,
                reason = sprintf(
                  "%d components; largest holds %.1f%% (dominance pass)",
                  k, 100 * largest_frac),
                n_components = k, sizes = sizes_sorted,
                size_break = sizes_sorted[1],
                largest_frac = largest_frac, gap_factor = NA_real_))
  }

  ## --- upper-tail gap test on log-sizes -------------------------
  ## The bulk (the small components) sets the dispersion baseline;
  ## the gap to the largest component must exceed gap_factor * bulk
  ## dispersion to count as a genuine separation rather than the
  ## tail of a continuous size distribution.  Same broken-stick logic
  ## .gap_threshold_lower uses on per-step speeds.
  log_sizes_desc <- log(sizes_sorted)
  bulk           <- log_sizes_desc[-1]                # all but the largest
  gap_above_bulk <- log_sizes_desc[1] - bulk[1]       # log(largest) - log(second)
  bulk_disp      <- stats::mad(bulk, constant = 1)    # robust scale of bulk
  if (!is.finite(bulk_disp) || bulk_disp <= 0) {
    bulk_disp <- stats::sd(bulk)
  }
  if (!is.finite(bulk_disp) || bulk_disp <= 0) {
    bulk_disp <- 1e-6                                 # all bulk identical
  }
  ratio <- gap_above_bulk / bulk_disp

  if (!is.finite(ratio) || ratio < gap_factor) {
    return(list(allow = FALSE,
                reason = sprintf(
                  "no clean gap above bulk on log-sizes (ratio %.2f < %.1f) -- size distribution is continuous, cut severs trajectory",
                  ratio, gap_factor),
                n_components = k, sizes = sizes_sorted,
                size_break = NA_integer_, largest_frac = largest_frac,
                gap_factor = ratio))
  }

  list(allow = TRUE,
       reason = sprintf(
         "%d components; largest holds %.1f%%, log-gap above bulk = %.2fx dispersion (>= %.1f)",
         k, 100 * largest_frac, ratio, gap_factor),
       n_components = k, sizes = sizes_sorted,
       size_break = sizes_sorted[1],
       largest_frac = largest_frac, gap_factor = ratio)
}


## Apply block expansion: assign block ids to fixes whose component
## has size strictly below the gate's size_break.  Always preserves
## the largest component as the trajectory.  Returns a length-n
## integer vector with NA where a fix is not in any expanded block.
##
## @keywords internal
.apply_block_expansion <- function(part, gate) {
  block_id <- rep(NA_integer_, part$n)
  if (!isTRUE(gate$allow)) return(block_id)
  sizes <- part$sizes
  if (length(sizes) < 2L) return(block_id)
  largest <- which.max(sizes)
  small_ids <- which(sizes < gate$size_break & seq_along(sizes) != largest)
  if (length(small_ids) == 0L) return(block_id)
  block_num <- 1L
  for (cid in small_ids) {
    members <- part$kept[part$comp_kept == cid]
    block_id[members] <- block_num
    block_num <- block_num + 1L
  }
  block_id
}


## Diagnostic plot for mt_clean_track results.
##
## @keywords internal
.plot_clean_track <- function(x) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(3, 3, 3, 1))

  cc <- sf::st_coordinates(x)
  graphics::plot(cc, type = "l", col = "grey70", lwd = 0.3,
                 asp = 1,
                 xlab = "", ylab = "",
                 main = sprintf("mt_clean_track -- %d flagged (%.2f%%)",
                                sum(x$is_outlier, na.rm = TRUE),
                                100 * mean(x$is_outlier, na.rm = TRUE)))
  flagged <- which(x$is_outlier)
  if (length(flagged)) {
    col <- ifelse(is.na(x$block_id[flagged]), "red", "orange")
    graphics::points(cc[flagged, , drop = FALSE],
                     col = grDevices::adjustcolor(col, 0.6),
                     pch = 20, cex = 0.6)
    graphics::legend("topright",
                     legend = c("individual outlier", "block outlier"),
                     col = c("red", "orange"), pch = 20, bty = "n")
  }
  invisible(NULL)
}


## Orchestrator-level pool sweep.
##
## After the per-track cascade converges, runs each pool-aware
## primitive wrapper once on the ORIGINAL multi-track input with
## pool_by set, captures the pool-added flag mask via the
## `pool_added` attribute exposed by `.apply_pool_union`, and
## unions those pool-added flags into the cascade's is_outlier.
## Pool flags are tagged `error_class = "pool"` for any fix whose
## class field was previously empty / NA / "none".  Per-detector
## flagged_by_* fingerprints are also updated where the pool
## helper writes to them.
##
## Primitive configurations match the cascade's internal calls
## where pool_by is meaningful:
##   - bridge:    method = "combined", threshold_type = "entropy",
##                threshold = 0.3, residual_floor = residual_floor
##   - speed_cap: threshold_type = "auto"
##   - detour:    threshold_type = "auto" (cascade uses "fixed" so
##                its detour-pool would be a no-op; pool sweep
##                uses "auto" so users actually get pool semantics
##                on detour at the orchestrator level).
## prob (mt_flag_outliers) is NOT included in v0.3 first cut --
## cascade and wrapper configurations diverge enough that pool
## semantics on prob via the orchestrator deserves a dedicated
## design pass.  Users wanting prob pooling can call
## mt_flag_outliers() standalone with pool_by; the cascade's
## prob contribution is preserved per-track regardless.
##
## Row-order alignment: both `out` (cascade) and the pool
## primitive's output use the same input geometry.  mt_track_id()
## values are the join key when row order may drift between
## primitive dispatch and cascade dispatch; this helper uses
## `move2::mt_track_id` to align.
##
## @keywords internal
.apply_orchestrator_pool_sweep <- function(out, x_orig, pool_by,
                                            location_error, residual_floor,
                                            detour_k,
                                            entropy_threshold = NULL,
                                            gap_threshold     = NULL,
                                            silent = FALSE) {
  say <- function(...) if (!silent) message(...)
  n_in <- nrow(out)
  pool_added_total <- logical(n_in)
  pool_added_bridge <- logical(n_in)
  pool_added_speed  <- logical(n_in)
  pool_added_detour <- logical(n_in)

  ## Build a row index aligning pool primitive output rows to `out`
  ## rows via (track_id, within-track row order) since both share
  ## the same input geometry and per-track dispatch is order-stable.
  align_to_out <- function(prim_out, pool_added_prim) {
    if (is.null(pool_added_prim)) return(logical(n_in))
    if (nrow(prim_out) != n_in) {
      ## Shouldn't happen for primitives that don't drop rows, but
      ## guard against future breakage by joining on track_id +
      ## within-track index.
      ids_out  <- as.character(move2::mt_track_id(out))
      ids_prim <- as.character(move2::mt_track_id(prim_out))
      idx_out  <- stats::ave(seq_along(ids_out),  ids_out,  FUN = seq_along)
      idx_prim <- stats::ave(seq_along(ids_prim), ids_prim, FUN = seq_along)
      key_out  <- paste(ids_out,  idx_out,  sep = "::")
      key_prim <- paste(ids_prim, idx_prim, sep = "::")
      m <- match(key_out, key_prim)
      aligned <- logical(n_in)
      ok <- !is.na(m)
      aligned[ok] <- pool_added_prim[m[ok]]
      return(aligned)
    }
    ## Same nrow -- check track_id agreement; if it matches by row,
    ## the lift order was preserved.
    ids_out  <- as.character(move2::mt_track_id(out))
    ids_prim <- as.character(move2::mt_track_id(prim_out))
    if (identical(ids_out, ids_prim)) {
      return(pool_added_prim)
    }
    ## Different inter-track order: align via (track_id, within-track idx).
    idx_out  <- stats::ave(seq_along(ids_out),  ids_out,  FUN = seq_along)
    idx_prim <- stats::ave(seq_along(ids_prim), ids_prim, FUN = seq_along)
    key_out  <- paste(ids_out,  idx_out,  sep = "::")
    key_prim <- paste(ids_prim, idx_prim, sep = "::")
    m <- match(key_out, key_prim)
    aligned <- logical(n_in)
    ok <- !is.na(m)
    aligned[ok] <- pool_added_prim[m[ok]]
    aligned
  }

  ## ---- Bridge pool sweep ----
  bp <- tryCatch(
    suppressMessages(suppressWarnings(
      mt_flag_outliers_bridge(x_orig,
                                method         = "combined",
                                threshold_type = "entropy",
                                ## entropy_threshold = NULL defers to
                                ## .entropy_threshold_lower's leaf
                                ## formal (single source of truth).
                                threshold      = entropy_threshold,
                                location_error = location_error,
                                residual_floor = residual_floor,
                                pool_by        = pool_by,
                                plot           = FALSE,
                                silent         = TRUE))),
    error = function(e) { say("  pool sweep: bridge skipped (", conditionMessage(e), ")"); NULL })
  if (!is.null(bp)) {
    pa <- attr(bp, "pool_added")
    pool_added_bridge <- align_to_out(bp, pa)
  }

  ## ---- Speed-cap pool sweep ----
  sp <- tryCatch(
    suppressMessages(suppressWarnings(
      mt_flag_speed_cap(x_orig,
                          threshold_type = "auto",
                          pool_by        = pool_by,
                          plot           = FALSE,
                          silent         = TRUE))),
    error = function(e) { say("  pool sweep: speed_cap skipped (", conditionMessage(e), ")"); NULL })
  if (!is.null(sp)) {
    pa <- attr(sp, "pool_added")
    pool_added_speed <- align_to_out(sp, pa)
  }

  ## ---- Detour pool sweep (uses threshold_type = "auto") ----
  dp <- tryCatch(
    suppressMessages(suppressWarnings(
      mt_flag_outliers_detour(x_orig,
                                k              = detour_k,
                                threshold_type = "auto",
                                min_leg        = 0,
                                pool_by        = pool_by,
                                plot           = FALSE,
                                silent         = TRUE))),
    error = function(e) { say("  pool sweep: detour skipped (", conditionMessage(e), ")"); NULL })
  if (!is.null(dp)) {
    pa <- attr(dp, "pool_added")
    pool_added_detour <- align_to_out(dp, pa)
  }

  pool_added_total <- pool_added_bridge | pool_added_speed | pool_added_detour

  ## Union pool flags into is_outlier and update fingerprints.
  prior_is_outlier <- out$is_outlier
  prior_is_outlier[is.na(prior_is_outlier)] <- FALSE
  newly_added <- pool_added_total & !prior_is_outlier
  out$is_outlier <- prior_is_outlier | pool_added_total

  ## Per-detector flagged_by_* fingerprints
  if ("flagged_by_bridge" %in% names(out)) {
    out$flagged_by_bridge <- out$flagged_by_bridge | pool_added_bridge
  }
  if ("flagged_by_speed" %in% names(out)) {
    out$flagged_by_speed <- out$flagged_by_speed | pool_added_speed
  }
  if ("flagged_by_detour" %in% names(out)) {
    out$flagged_by_detour <- out$flagged_by_detour | pool_added_detour
  }

  ## Tag error_class for newly-flagged fixes.  Use "pool" when the
  ## prior error_class was missing / empty / "none"; otherwise
  ## prepend "pool+" to the existing class (preserves cascade
  ## attribution).
  if ("error_class" %in% names(out)) {
    ec     <- out$error_class
    blank  <- is.na(ec) | !nzchar(ec) | ec %in% c("none")
    newly_pool <- newly_added & blank
    if (any(newly_pool)) ec[newly_pool] <- "pool"
    out$error_class <- ec
  }

  n_total_added <- sum(newly_added)
  if (n_total_added > 0L) {
    say(sprintf(
      "Pool sweep: +%d new flag(s) (bridge=%d, speed=%d, detour=%d).",
      n_total_added, sum(pool_added_bridge & newly_added),
      sum(pool_added_speed  & newly_added),
      sum(pool_added_detour & newly_added)))
  } else if (!silent) {
    say("Pool sweep: 0 new flag(s).")
  }
  out
}


## Normalise a `state` argument into a per-fix character vector aligned
## with x.  NULL passes through.  NA entries are preserved (treated as
## their own state by the segmenter).
##
## @keywords internal
.resolve_state <- function(state, x) {
  if (is.null(state)) return(NULL)

  if (is.character(state) && length(state) == 1L) {
    if (!state %in% names(x)) {
      rlang::abort(
        sprintf("`state` column '%s' not found on x.", state),
        class = "move2utils_mt_clean_track_state_missing_column")
    }
    vals <- x[[state]]
  } else {
    if (length(state) != nrow(x)) {
      rlang::abort(
        sprintf("`state` vector length (%d) must equal nrow(x) (%d).",
                length(state), nrow(x)),
        class = "move2utils_mt_clean_track_state_bad_length")
    }
    vals <- state
  }

  ## Coerce factors / numerics / logicals to character so RLE keys are
  ## comparable.  NA is preserved.
  if (is.factor(vals)) vals <- as.character(vals)
  as.character(vals)
}


## Per-segment dispatch: split x into runs of constant (track_id, state)
## and recurse into mt_clean_track on each segment with state = NULL.
## Recombines flag columns in original row order.
##
## @keywords internal
.dispatch_by_state <- function(x, state_vec, v_max_per_track, v_max,
                                consensus, consensus_custom,
                                transition_buffer,
                                use_detour, detour_k,
                                detour_threshold,
                                iterations, max_iterations,
                                max_flag_fraction, expand_blocks,
                                pre_peel_aux, persistence_filter,
                                location_error, residual_floor, step_floor,
                                pool_by = NULL,
                                bridge_method         = NULL,
                                bridge_threshold_type = NULL,
                                bridge_iterations     = NULL,
                                prob_threshold_type   = NULL,
                                detour_threshold_type = NULL,
                                entropy_threshold            = NULL,
                                gap_threshold                = NULL,
                                persistence_filter_threshold = NULL,
                                plot, remove, silent, compact) {
  say <- .say(silent)

  ## Caller's columns, captured before flag columns are pre-allocated
  ## below -- stripped on the remove = TRUE return (see .strip_cascade_cols
  ## and the orig_cols note in mt_clean_track()).
  orig_cols <- names(x)
  ids <- as.character(move2::mt_track_id(x))
  n   <- nrow(x)

  ## NA in state is its own state value via a sentinel.
  state_key <- ifelse(is.na(state_vec), "<<NA>>", state_vec)
  combined  <- paste(ids, state_key, sep = "@@")
  is_change <- c(TRUE, utils::head(combined, -1) !=
                       utils::tail(combined, -1))
  seg_id    <- cumsum(is_change)
  unique_segs <- unique(seg_id)

  say(sprintf(
    "State dispatch: %d (track x state) segments across %d tracks.",
    length(unique_segs), length(unique(ids))))

  ## Pre-allocate flag columns on the parent.  We write each segment's
  ## flag columns back to the corresponding rows of `x` rather than
  ## rbind'ing the segment move2s -- move2's rbind rejects duplicated
  ## track ids across the bound objects, which always trips when a
  ## single track is split into multiple state segments.
  x$is_outlier        <- logical(n)
  x$flagged_by_bridge <- logical(n)
  x$flagged_by_prob   <- logical(n)
  x$flagged_by_speed  <- logical(n)
  x$flagged_by_detour <- logical(n)
  x$flag_iteration    <- rep(NA_integer_, n)
  x$block_id          <- rep(NA_integer_, n)
  x$error_class       <- rep(NA_character_, n)
  x$error_classes_all <- rep("",            n)
  x$combined_evidence <- rep(NA_real_,      n)  # populated only by evidence modes; dropped below if unused

  block_offset <- 0L

  ## Suppress the per-segment "no physiological cap" advice repetition.
  ## The user supplied state =, which is a sophisticated mode of use;
  ## emitting the auto-cap newcomer note once per segment is noise.
  ## Restored on exit so a later mt_clean_track call without state =
  ## still gets the advice as expected.
  old_advice_opt <- getOption("move2utils.suppress_first_run_advice", FALSE)
  options(move2utils.suppress_first_run_advice = TRUE)
  on.exit(options(move2utils.suppress_first_run_advice = old_advice_opt),
          add = TRUE)

  ## Aggregate count of segments skipped for size reasons; reported as
  ## a single summary at the end of the dispatch instead of one warning
  ## per segment from the inner mt_clean_track().  The inner cascade's
  ## hard minimum is n = 10; below that it warns and returns unflagged.
  n_skipped_tiny  <- 0L   # < 3 fixes  (already handled below)
  n_skipped_short <- 0L   # 3 <= n < 10
  min_for_cascade <- 10L

  for (k in seq_along(unique_segs)) {
    s    <- unique_segs[k]
    rows <- which(seg_id == s)
    xs   <- x[rows, ]
    tid  <- ids[rows[1L]]
    sval <- state_vec[rows[1L]]

    if (length(rows) < 3L) {
      n_skipped_tiny <- n_skipped_tiny + 1L
      next  # already FALSE / NA from pre-allocation
    }
    if (length(rows) < min_for_cascade) {
      ## Cascade emits its own "Too few locations" warning for n<10.
      ## Skip here to keep the warning surface clean; count for the
      ## summary below.
      n_skipped_short <- n_skipped_short + 1L
      next
    }

    ## Per-track v_max: if mass/mode produced a named vector, look it
    ## up; otherwise inherit the scalar v_max passed in.
    vmi <- if (!is.null(v_max_per_track)) {
      as.numeric(v_max_per_track[[tid]])
    } else v_max

    ## Slice location_error if given as a per-fix numeric vector.
    oei <- if (is.numeric(location_error) && length(location_error) > 1L) {
      location_error[rows]
    } else location_error

    seg <- mt_clean_track(
      xs,
      v_max             = vmi,
      mass              = NULL,
      mode              = NULL,
      state             = NULL,
      consensus         = consensus,
      consensus_custom  = consensus_custom,
      transition_buffer = 0L,   # transitions are applied at this dispatch level
      use_detour        = use_detour,
      detour_k          = detour_k,
      pool_by           = NULL, # outer-only; pool sweep runs after dispatch
      bridge_method         = bridge_method,
      bridge_threshold_type = bridge_threshold_type,
      bridge_iterations     = bridge_iterations,
      prob_threshold_type   = prob_threshold_type,
      detour_threshold_type = detour_threshold_type,
      entropy_threshold            = entropy_threshold,
      gap_threshold                = gap_threshold,
      persistence_filter_threshold = persistence_filter_threshold,
      detour_threshold  = detour_threshold,
      iterations        = iterations,
      max_iterations    = max_iterations,
      max_flag_fraction = max_flag_fraction,
      expand_blocks     = expand_blocks,
      pre_peel_aux      = pre_peel_aux,
      persistence_filter = persistence_filter,
      location_error         = oei,
      residual_floor    = residual_floor,
      step_floor        = step_floor,
      plot              = FALSE,
      remove            = FALSE,
      silent            = silent,
      compact           = compact)

    ## Renumber block_id with a global offset so segment-local block
    ## ids stay disjoint after recombination.
    bi <- seg$block_id
    if (any(!is.na(bi))) {
      bi[!is.na(bi)] <- bi[!is.na(bi)] + block_offset
      block_offset  <- max(bi, na.rm = TRUE)
    }

    x$is_outlier[rows]        <- seg$is_outlier
    x$flagged_by_bridge[rows] <- seg$flagged_by_bridge
    x$flagged_by_prob[rows]   <- seg$flagged_by_prob
    x$flagged_by_speed[rows]  <- seg$flagged_by_speed
    x$flagged_by_detour[rows] <- seg$flagged_by_detour
    x$flag_iteration[rows]    <- seg$flag_iteration
    x$block_id[rows]          <- bi
    x$error_class[rows]       <- seg$error_class
    if (!is.null(seg$error_classes_all))
      x$error_classes_all[rows] <- seg$error_classes_all
    if (!is.null(seg$combined_evidence))
      x$combined_evidence[rows] <- seg$combined_evidence
  }
  ## drop the evidence column under Boolean consensus modes (no segment
  ## populated it), matching the single-track path's behaviour.
  if (all(is.na(x$combined_evidence))) x$combined_evidence <- NULL

  ## Single-line summary for short / too-short segments.  Replaces a
  ## warning-storm of "Too few locations" messages from the inner
  ## cascade calls.
  total_skipped <- n_skipped_short + n_skipped_tiny
  if (total_skipped > 0L) {
    detail <- if (n_skipped_short && n_skipped_tiny) {
      sprintf(" (%d with <3 fixes, %d with 3-9 fixes)",
              n_skipped_tiny, n_skipped_short)
    } else if (n_skipped_tiny) {
      sprintf(" (%d with <3 fixes)", n_skipped_tiny)
    } else {
      sprintf(" (%d with 3-9 fixes)", n_skipped_short)
    }
    say(sprintf(
      "State-dispatch summary: %d of %d segments were too short for the cascade and left unflagged%s.",
      total_skipped, length(unique_segs), detail))
  }

  ## --- combine-with-global pass --------------------------------
  ## State-conditional dispatch fits each primitive's data-driven
  ## threshold to its segment's distribution.  This is right for
  ## kinematic primitives (speed, prob) -- per-state distributions
  ## are the relevant context.  But it is WRONG for geometric
  ## primitives (bridge, detour) whose "impossible" doesn't depend
  ## on state: bridge's threshold inside a migrate segment is
  ## inflated by legitimate flight residuals, so colony spikes
  ## embedded in migrate sub-trains slip through.
  ##
  ## To respect each primitive's scope, also run mt_clean_track on
  ## the FULL input in single-state mode and UNION the per-detector
  ## flags with the per-segment results.  Each mode then catches
  ## what its scope of validity allows, and the combined set is the
  ## union: geometric impossibility flagged by the global run +
  ## state-anomalous kinematics flagged by the per-segment runs.
  ## error_class is re-computed from the combined per-detector flags.
  ## Combine-with-global pass: state-conditional dispatch is right
  ## for kinematic primitives (per-state distributions are the relevant
  ## context) but wrong for geometric primitives whose "impossible"
  ## doesn't depend on state.  Run a single-state baseline alongside
  ## and UNION the per-detector flags; the consensus then operates on
  ## the union via re-running .classify_flags below.
  do_combine_with_global <- TRUE
  if (do_combine_with_global) {
    say("Combine-with-global: running single-state baseline alongside state-conditional...")
    baseline <- suppressMessages(mt_clean_track(
      x,
      v_max             = v_max,
      mass              = NULL,    # already resolved into v_max via per-track
      mode              = NULL,
      state             = NULL,
      consensus         = consensus,
      consensus_custom  = consensus_custom,
      transition_buffer = 0L,
      use_detour        = use_detour,
      detour_k          = detour_k,
      detour_threshold  = detour_threshold,
      iterations        = iterations,
      max_iterations    = max_iterations,
      max_flag_fraction = max_flag_fraction,
      expand_blocks     = expand_blocks,
      pre_peel_aux      = pre_peel_aux,
      persistence_filter = persistence_filter,
      location_error         = location_error,
      residual_floor    = residual_floor,
      step_floor        = step_floor,
      pool_by           = NULL,  # pool sweep runs once at dispatch level
      bridge_method         = bridge_method,
      bridge_threshold_type = bridge_threshold_type,
      bridge_iterations     = bridge_iterations,
      prob_threshold_type   = prob_threshold_type,
      detour_threshold_type = detour_threshold_type,
      entropy_threshold            = entropy_threshold,
      gap_threshold                = gap_threshold,
      persistence_filter_threshold = persistence_filter_threshold,
      plot              = FALSE,
      remove            = FALSE,
      silent            = silent,
      compact           = compact))

    ## Union the per-detector fingerprints
    x$flagged_by_bridge <- x$flagged_by_bridge | baseline$flagged_by_bridge
    x$flagged_by_detour <- x$flagged_by_detour | baseline$flagged_by_detour
    x$flagged_by_prob   <- x$flagged_by_prob   | baseline$flagged_by_prob
    x$flagged_by_speed  <- x$flagged_by_speed  | baseline$flagged_by_speed
    ## Combine block_id: prefer existing per-segment block, fall back to baseline
    base_block <- baseline$block_id
    new_blocks <- !is.na(base_block) & is.na(x$block_id)
    if (any(new_blocks))
      x$block_id[new_blocks] <- base_block[new_blocks] + block_offset
    ## Combine flag_iteration: keep the earliest non-NA
    fi_x <- x$flag_iteration
    fi_b <- baseline$flag_iteration
    take_b <- is.na(fi_x) & !is.na(fi_b)
    x$flag_iteration[take_b] <- fi_b[take_b]
    ## Now is_outlier is the union and error_class is recomputed
    x$is_outlier <- x$is_outlier | baseline$is_outlier
    classified <- .classify_flags(
      x$is_outlier,
      x$flagged_by_bridge, x$flagged_by_prob, x$flagged_by_speed,
      x$block_id, x$flagged_by_detour,
      flag_iteration = x$flag_iteration)
    x$error_class       <- classified$error_class
    x$error_classes_all <- classified$error_classes_all
    say(sprintf("Combine-with-global: state-only=%d  +  global-only=%d  =>  union=%d.",
                sum(!baseline$is_outlier &
                    !is.na(x$error_class)),  # rough -- some are both
                sum( baseline$is_outlier & !x$error_class %in%
                     c("geometric_spike", "consensus")),
                sum(x$is_outlier)))
  }

  ## --- transition-zone demotion ---------------------------------
  ## At state transitions, the kinematic distributions on either side
  ## may both be defensible interpretations.  A fix flagged by a
  ## class that depends on local state context (state_anomaly,
  ## kinematic_confluence) inside the transition zone is demoted to
  ## "state_transition_buffered" (kept) UNLESS it also fired a
  ## class that doesn't depend on state (consensus, geometric_spike,
  ## block, physiological).
  if (transition_buffer > 0L) {
    ## state changes: a fix at index i is at a transition if
    ## state[i] != state[i-1] for any track-internal pair.
    is_state_change <- c(FALSE,
                         state_key[-1L] != state_key[-length(state_key)] &
                         ids[-1L]       == ids[-length(ids)])
    transition_zone <- rep(FALSE, n)
    chg_idx <- which(is_state_change)
    for (offset in seq(-transition_buffer, transition_buffer)) {
      hits <- chg_idx + offset
      hits <- hits[hits >= 1L & hits <= n]
      transition_zone[hits] <- TRUE
    }
    ## state-dependent classes vs robust classes
    state_dep_classes <- c("state_anomaly", "kinematic_confluence")
    is_state_dep     <- !is.na(x$error_class) &
                        x$error_class %in% state_dep_classes
    demote <- transition_zone & is_state_dep
    if (any(demote)) {
      n_demoted <- sum(demote)
      x$is_outlier[demote]  <- FALSE
      x$error_class[demote] <- "state_transition_buffered"
      x$flag_iteration[demote] <- NA_integer_
      say(sprintf(
        "Transition buffer (k=%d): demoted %d state-dependent flag(s) at state boundaries.",
        transition_buffer, n_demoted))
    }
  }

  ## ---- orchestrator pool sweep (state-dispatch path) ----
  ## Pool sweep runs after the state-conditional cascade converges
  ## and after the transition-zone demotion -- pool-added flags
  ## bypass the state-dependent demotion since they are derived
  ## from pool-fitted (state-independent) thresholds.
  if (!is.null(pool_by)) {
    x <- .apply_orchestrator_pool_sweep(
            x, x, pool_by,
            location_error = location_error,
            residual_floor = residual_floor,
            detour_k       = detour_k,
            entropy_threshold = entropy_threshold,
            gap_threshold     = gap_threshold,
            silent         = silent)
  }

  if (plot)   .plot_clean_track(x)
  if (remove) x <- .strip_cascade_cols(x[!x$is_outlier, ], orig_cols)
  x
}
