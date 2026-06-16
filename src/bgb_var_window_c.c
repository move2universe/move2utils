/*
 * bgb_var_window_c.c -- C implementation of the dBGB per-window
 * variance estimator with breakpoint detection.
 *
 * Mirrors the structure of bm_variance_window_c (dBBMM): one .Call
 * per sliding-window position, the entire breakpoint search runs in
 * C.  The earlier R implementation (.bgb_var_break in
 * R/dbgb_variance.R) drove an outer loop over candidate breakpoints
 * with optim() per candidate, paying ~30-40 R-level optim() calls per
 * window and hundreds of R<->C boundary crossings inside each
 * optim().  On a 1748-fix track that scaled to ~200x slower than
 * dBBMM.  This kernel performs the same sequential search (no break /
 * para-only / orth-only / cross-conditional refinement) but does it
 * entirely in native code.
 *
 * Model.  The dBGB negative log-likelihood with diagonal covariance
 * Sigma = diag(s_para, s_orth) is separable across axes,
 *
 *   logL = sum_i [ -log(2*pi) - 0.5 log(s_para_i) - 0.5 log(s_orth_i)
 *                  - r_para_i^2/(2 s_para_i) - r_orth_i^2/(2 s_orth_i) ]
 *
 * so each axis sigma is the maximiser of an independent 1D log-
 * likelihood -- exactly the same 1D problem dBBMM solves with
 * Brent's method.  We drop the per-pair -log(2*pi) constant; it is
 * shared by all candidate models in the BIC comparison so cancels
 * in differences.  Absolute log-likelihood values therefore differ
 * by n_pairs * log(2*pi) from R llBGBvar reference, but model
 * selection and fitted sigmas match up to optimiser tolerance.
 *
 * Breakpoint convention.  Following the R reference we use a single
 * LOO grid over the full window and flip sigma at the break: pair
 * indices [0, k) use sigma_before, indices [k, n_pairs) use
 * sigma_after, where k = floor(b_fix / 2) for break position b_fix
 * (1-indexed in fix space).  This is *not* the dBBMM convention of
 * re-precomputing LOO on prefix and suffix segments; it is the
 * convention used historically by .bgb_var_break and matched here
 * to keep the C kernel a drop-in replacement.
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * 1D negative log-likelihood for a single axis (parallel or
 * orthogonal) given the leave-one-out quantities and the per-pair
 * squared residual on that axis, evaluated over pair indices
 * [start, end).
 *
 *   v_i = T_jump_i * alpha_i * (1 - alpha_i) * sigma2
 *       + (1 - alpha_i)^2 * le1_i^2
 *       + alpha_i^2       * le2_i^2
 *
 *   nll = sum_i [ 0.5 * log(v_i) + ztz_axis_i / (2 * v_i) ]
 */
static double neg_log_lik_axis(double var,
                               const double *T_jump,
                               const double *alpha,
                               const double *le1,
                               const double *le2,
                               const double *ztz_axis,
                               int n_pairs) {
    double nll = 0.0;
    for (int i = 0; i < n_pairs; i++) {
        double v = T_jump[i] * alpha[i] * (1.0 - alpha[i]) * var
                 + (1.0 - alpha[i]) * (1.0 - alpha[i]) * le1[i] * le1[i]
                 + alpha[i]         * alpha[i]         * le2[i] * le2[i];
        if (v <= 0.0) v = DBL_EPSILON;
        nll += 0.5 * log(v) + ztz_axis[i] / (2.0 * v);
    }
    return nll;
}

/*
 * Brent's method 1D minimiser for the per-axis negative log-
 * likelihood on [a, b].  Structurally identical to brent_min in
 * bm_variance_c.c (same algorithm, tolerances and iteration cap);
 * the only difference is the objective function it calls.
 */
static double brent_min_axis(double a, double b,
                             const double *T_jump,
                             const double *alpha,
                             const double *le1,
                             const double *le2,
                             const double *ztz_axis,
                             int n_pairs,
                             double tol, double *f_min) {
    double x, w, v, fx, fw, fv, e, d, u, fu;
    double midpoint, tol1, tol2, p, q, r;
    const double golden = 0.3819660;
    const int max_iter = 500;

    x = w = v = a + golden * (b - a);
    fx = fw = fv = neg_log_lik_axis(x, T_jump, alpha, le1, le2, ztz_axis, n_pairs);
    e = 0.0;
    d = 0.0;

    for (int iter = 0; iter < max_iter; iter++) {
        midpoint = 0.5 * (a + b);
        tol1 = tol * fabs(x) + 1e-10;
        tol2 = 2.0 * tol1;

        if (fabs(x - midpoint) <= (tol2 - 0.5 * (b - a))) {
            *f_min = fx;
            return x;
        }

        if (fabs(e) > tol1) {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p; else q = -q;
            r = e;
            e = d;

            if (fabs(p) < fabs(0.5 * q * r) && p > q * (a - x) && p < q * (b - x)) {
                d = p / q;
                u = x + d;
                if ((u - a) < tol2 || (b - u) < tol2)
                    d = (x < midpoint) ? tol1 : -tol1;
            } else {
                e = (x < midpoint) ? b - x : a - x;
                d = golden * e;
            }
        } else {
            e = (x < midpoint) ? b - x : a - x;
            d = golden * e;
        }

        u = (fabs(d) >= tol1) ? x + d : x + ((d > 0) ? tol1 : -tol1);
        fu = neg_log_lik_axis(u, T_jump, alpha, le1, le2, ztz_axis, n_pairs);

        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    *f_min = fx;
    return x;
}

/*
 * Pre-compute the leave-one-out quantities and per-pair squared
 * para/orth residuals for a track segment.  Matches the R reference
 * .bgb_var_break: pairs are at every-other index (i = 1, 3, 5, ...
 * in 0-indexed C, i.e. fix indices 2, 4, 6, ... in 1-indexed R-
 * speak), and the decomposition axis at each pair is mu -> next_fix
 * (the local travel direction) -- the same convention as
 * .delta_para_orth in dbgb_variance.R.
 *
 * Degenerate axis (mu == next_fix): the squared residual is split
 * equally between para and orth, mirroring same_loc_dir handling
 * in .delta_para_orth.
 */
static int precompute_loo_dbgb(const double *x, const double *y,
                                const double *time_lag,
                                const double *loc_err,
                                int n,
                                double *T_jump,
                                double *alpha,
                                double *ztz_para,
                                double *ztz_orth,
                                double *le1,
                                double *le2) {
    int n_pairs = 0;
    int i = 1;
    while (i < n - 1) {
        double t = time_lag[i - 1] + time_lag[i];
        T_jump[n_pairs] = t;
        double a = time_lag[i - 1] / t;
        alpha[n_pairs] = a;

        double ux = x[i - 1] + a * (x[i + 1] - x[i - 1]);
        double uy = y[i - 1] + a * (y[i + 1] - y[i - 1]);
        double rx = x[i] - ux;
        double ry = y[i] - uy;
        double rn2 = rx * rx + ry * ry;

        double axx = x[i + 1] - ux;
        double axy = y[i + 1] - uy;
        double axn2 = axx * axx + axy * axy;

        double rp2, ro2;
        if (axn2 > 0.0) {
            double dot = rx * axx + ry * axy;
            rp2 = (dot * dot) / axn2;       /* (|r . axhat|)^2 */
            double diff = rn2 - rp2;
            ro2 = (diff > 0.0) ? diff : 0.0;
        } else {
            rp2 = rn2 * 0.5;
            ro2 = rn2 * 0.5;
        }
        ztz_para[n_pairs] = rp2;
        ztz_orth[n_pairs] = ro2;

        le1[n_pairs] = loc_err[i - 1];
        le2[n_pairs] = loc_err[i + 1];
        n_pairs++;
        i += 2;
    }
    return n_pairs;
}

/*
 * .Call entry point for one window position.
 *
 * Inputs (all REALSXP / INTSXP per the caller in R):
 *   r_x, r_y         length-n numeric coordinates
 *   r_time_lag       length-n numeric (time_mins[i+1] - time_mins[i]);
 *                    the last element is unused
 *   r_loc_err        length-n numeric location error (1-sigma)
 *   r_margin_breaks  integer vector of 1-indexed candidate break
 *                    positions in fix space, pre-filtered by the
 *                    caller to the odd-1-indexed margin range
 *
 * Returns a list with elements:
 *   paraSd  numeric[n]   per-fix parallel SD (NA outside the
 *                        [min(margin_breaks), max(margin_breaks))
 *                        range, matching the R reference)
 *   orthSd  numeric[n]   per-fix orthogonal SD (same NA mask)
 *   para_break  integer  selected fix-space para break, or 0
 *   orth_break  integer  selected fix-space orth break, or 0
 *
 * The four-step search exactly mirrors the R reference
 * .bgb_var_break: (1) no break, (2) try each break as a para-only
 * break, (3) as an orth-only break, (4) cross-conditional pass --
 * conditional on the chosen para break, refine orth (and vice
 * versa).  At each step we compare BIC = 2*nll + k*log(n) where
 * k is the parameter count (2, 3, 3, or 4).
 */
SEXP bgb_var_window_c(SEXP r_x, SEXP r_y, SEXP r_time_lag, SEXP r_loc_err,
                      SEXP r_margin_breaks) {
    int n = length(r_x);
    double *x        = REAL(r_x);
    double *y        = REAL(r_y);
    double *time_lag = REAL(r_time_lag);
    double *loc_err  = REAL(r_loc_err);
    int *margin_breaks = INTEGER(r_margin_breaks);
    int n_breaks = length(r_margin_breaks);

    if (n < 3) {
        error("bgb_var_window_c: window must contain at least 3 fixes (got %d)", n);
    }
    if (n_breaks < 1) {
        error("bgb_var_window_c: margin_breaks must contain at least one candidate");
    }

    /* Single LOO precompute for the full window. */
    int max_pairs = n / 2 + 1;
    double *Tj  = (double *)R_alloc(max_pairs, sizeof(double));
    double *al  = (double *)R_alloc(max_pairs, sizeof(double));
    double *le1 = (double *)R_alloc(max_pairs, sizeof(double));
    double *le2 = (double *)R_alloc(max_pairs, sizeof(double));
    double *zp  = (double *)R_alloc(max_pairs, sizeof(double));
    double *zo  = (double *)R_alloc(max_pairs, sizeof(double));

    int n_pairs = precompute_loo_dbgb(x, y, time_lag, loc_err, n,
                                       Tj, al, zp, zo, le1, le2);
    if (n_pairs < 1) {
        error("bgb_var_window_c: window too short for variance estimation");
    }

    double logn = log((double)n);

    /* ---- Step 1: no break ---- */
    double f_p_full, f_o_full;
    double s2_p_full = brent_min_axis(0.0, 1e15, Tj, al, le1, le2,
                                       zp, n_pairs, 1e-8, &f_p_full);
    double s2_o_full = brent_min_axis(0.0, 1e15, Tj, al, le1, le2,
                                       zo, n_pairs, 1e-8, &f_o_full);
    double bic_none = 2.0 * (f_p_full + f_o_full) + 2.0 * logn;

    double best_bic    = bic_none;
    int    best_pb     = 0;          /* 0 == no break on para axis */
    int    best_ob     = 0;          /* 0 == no break on orth axis */
    double best_s2_pb  = s2_p_full;
    double best_s2_pa  = s2_p_full;
    double best_s2_ob  = s2_o_full;
    double best_s2_oa  = s2_o_full;

    /* Cache per-axis breakpoint nll's so step 4 reuses step 2/3 work. */
    double *cache_fp     = (double *)R_alloc(n_breaks, sizeof(double));
    double *cache_s2_pb  = (double *)R_alloc(n_breaks, sizeof(double));
    double *cache_s2_pa  = (double *)R_alloc(n_breaks, sizeof(double));
    double *cache_fo     = (double *)R_alloc(n_breaks, sizeof(double));
    double *cache_s2_ob  = (double *)R_alloc(n_breaks, sizeof(double));
    double *cache_s2_oa  = (double *)R_alloc(n_breaks, sizeof(double));

    /* ---- Steps 2 & 3: para-only and orth-only breaks ---- */
    for (int bi = 0; bi < n_breaks; bi++) {
        int b_fix = margin_breaks[bi];           /* 1-indexed fix position */
        int k = b_fix / 2;                       /* pair-space split */

        if (k <= 0 || k >= n_pairs) {
            cache_fp[bi] = R_PosInf;
            cache_fo[bi] = R_PosInf;
            cache_s2_pb[bi] = NA_REAL;
            cache_s2_pa[bi] = NA_REAL;
            cache_s2_ob[bi] = NA_REAL;
            cache_s2_oa[bi] = NA_REAL;
            continue;
        }

        /* Para-only break: split sigma_para at k, orth axis stays full. */
        double f_pb, f_pa;
        double s2_pb = brent_min_axis(0.0, 1e15, Tj, al, le1, le2,
                                       zp, k, 1e-8, &f_pb);
        double s2_pa = brent_min_axis(0.0, 1e15,
                                       Tj + k, al + k, le1 + k, le2 + k,
                                       zp + k, n_pairs - k, 1e-8, &f_pa);
        double f_para_break = f_pb + f_pa;
        cache_fp[bi]    = f_para_break;
        cache_s2_pb[bi] = s2_pb;
        cache_s2_pa[bi] = s2_pa;

        double bic_para_only = 2.0 * (f_para_break + f_o_full) + 3.0 * logn;
        if (bic_para_only < best_bic) {
            best_bic    = bic_para_only;
            best_pb     = b_fix;
            best_ob     = 0;
            best_s2_pb  = s2_pb;
            best_s2_pa  = s2_pa;
            best_s2_ob  = s2_o_full;
            best_s2_oa  = s2_o_full;
        }

        /* Orth-only break: split sigma_orth at k, para axis stays full. */
        double f_ob, f_oa;
        double s2_ob = brent_min_axis(0.0, 1e15, Tj, al, le1, le2,
                                       zo, k, 1e-8, &f_ob);
        double s2_oa = brent_min_axis(0.0, 1e15,
                                       Tj + k, al + k, le1 + k, le2 + k,
                                       zo + k, n_pairs - k, 1e-8, &f_oa);
        double f_orth_break = f_ob + f_oa;
        cache_fo[bi]    = f_orth_break;
        cache_s2_ob[bi] = s2_ob;
        cache_s2_oa[bi] = s2_oa;

        double bic_orth_only = 2.0 * (f_p_full + f_orth_break) + 3.0 * logn;
        if (bic_orth_only < best_bic) {
            best_bic    = bic_orth_only;
            best_pb     = 0;
            best_ob     = b_fix;
            best_s2_pb  = s2_p_full;
            best_s2_pa  = s2_p_full;
            best_s2_ob  = s2_ob;
            best_s2_oa  = s2_oa;
        }
    }

    /* ---- Step 4: cross-conditional refinement ---- */
    /* If the current best has a para break, vary orth.  The two axes
     * are independent given fixed breaks, so the joint nll is just
     * the sum of the per-axis breakpoint nlls, both already cached. */
    if (best_pb != 0) {
        int para_idx = -1;
        for (int bi = 0; bi < n_breaks; bi++) {
            if (margin_breaks[bi] == best_pb) { para_idx = bi; break; }
        }
        if (para_idx >= 0) {
            double f_para_break = cache_fp[para_idx];
            for (int bi = 0; bi < n_breaks; bi++) {
                if (!R_FINITE(cache_fo[bi])) continue;
                double bic_both = 2.0 * (f_para_break + cache_fo[bi]) + 4.0 * logn;
                if (bic_both < best_bic) {
                    best_bic    = bic_both;
                    best_ob     = margin_breaks[bi];
                    best_s2_ob  = cache_s2_ob[bi];
                    best_s2_oa  = cache_s2_oa[bi];
                }
            }
        }
    }
    if (best_ob != 0) {
        int orth_idx = -1;
        for (int bi = 0; bi < n_breaks; bi++) {
            if (margin_breaks[bi] == best_ob) { orth_idx = bi; break; }
        }
        if (orth_idx >= 0) {
            double f_orth_break = cache_fo[orth_idx];
            for (int bi = 0; bi < n_breaks; bi++) {
                if (!R_FINITE(cache_fp[bi])) continue;
                double bic_both = 2.0 * (cache_fp[bi] + f_orth_break) + 4.0 * logn;
                if (bic_both < best_bic) {
                    best_bic    = bic_both;
                    best_pb     = margin_breaks[bi];
                    best_s2_pb  = cache_s2_pb[bi];
                    best_s2_pa  = cache_s2_pa[bi];
                }
            }
        }
    }

    /* ---- Build per-fix paraSd / orthSd ---- */
    SEXP r_para_sd = PROTECT(allocVector(REALSXP, n));
    SEXP r_orth_sd = PROTECT(allocVector(REALSXP, n));
    double *para_sd = REAL(r_para_sd);
    double *orth_sd = REAL(r_orth_sd);

    int min_brk = margin_breaks[0];
    int max_brk = margin_breaks[0];
    for (int bi = 1; bi < n_breaks; bi++) {
        if (margin_breaks[bi] < min_brk) min_brk = margin_breaks[bi];
        if (margin_breaks[bi] > max_brk) max_brk = margin_breaks[bi];
    }

    for (int i = 0; i < n; i++) {
        int fix_pos = i + 1;                     /* 1-indexed */
        if (fix_pos < min_brk || fix_pos >= max_brk) {
            para_sd[i] = NA_REAL;
            orth_sd[i] = NA_REAL;
            continue;
        }
        double s2p = (best_pb == 0 || fix_pos < best_pb) ? best_s2_pb : best_s2_pa;
        double s2o = (best_ob == 0 || fix_pos < best_ob) ? best_s2_ob : best_s2_oa;
        para_sd[i] = (s2p > 0.0) ? sqrt(s2p) : 0.0;
        orth_sd[i] = (s2o > 0.0) ? sqrt(s2o) : 0.0;
    }

    SEXP result = PROTECT(allocVector(VECSXP, 4));
    SEXP names  = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, mkChar("paraSd"));
    SET_STRING_ELT(names, 1, mkChar("orthSd"));
    SET_STRING_ELT(names, 2, mkChar("para_break"));
    SET_STRING_ELT(names, 3, mkChar("orth_break"));
    setAttrib(result, R_NamesSymbol, names);
    SET_VECTOR_ELT(result, 0, r_para_sd);
    SET_VECTOR_ELT(result, 1, r_orth_sd);
    SET_VECTOR_ELT(result, 2, ScalarInteger(best_pb));
    SET_VECTOR_ELT(result, 3, ScalarInteger(best_ob));

    UNPROTECT(4);
    return result;
}

/*
 * Track-level entry point for dBGB.  Sweeps every sliding window in
 * C, runs the same 4-step breakpoint search bgb_var_window_c does
 * per window, and aggregates the per-fix paraSd/orthSd estimates by
 * RMS across windows -- the same root-mean-square aggregation the
 * R-level driver did via aggregate(... ~ seg, FUN = function(x)
 * sqrt(mean(x^2))).
 *
 * Replaces R-level lapply over windows + do.call(rbind, ...) +
 * aggregate(...) with a single .Call.  Optionally parallelised over
 * windows via OpenMP; per-thread working buffers are heap-allocated
 * (R_alloc isn't thread-safe) and the per-fix accumulators are
 * updated under #pragma omp atomic.
 *
 * Inputs:
 *   r_x, r_y         length-n track coordinates
 *   r_time_mins      length-n timestamps in minutes (the dBGB R
 *                    reference uses time_mins, not time_lag, since
 *                    its breakpoint convention reuses the same
 *                    LOO grid; we mirror that)
 *   r_loc_err        length-n location error (1-sigma, in metres)
 *   r_window_size    odd integer
 *   r_margin         odd integer
 *
 * Returns list:
 *   para_sd  numeric[n]   per-fix RMS-aggregated parallel SD,
 *                          NA on the first/last (window_size-1)/2
 *                          fixes where no window contributes
 *   orth_sd  numeric[n]   ditto orthogonal axis
 *   n_estim  integer[n]   number of windows contributing to each fix
 */
SEXP bgb_var_break_track_c(SEXP r_x, SEXP r_y, SEXP r_time_mins,
                            SEXP r_loc_err,
                            SEXP r_window_size, SEXP r_margin) {
    int n          = length(r_x);
    int window_size = asInteger(r_window_size);
    int margin     = asInteger(r_margin);

    if (window_size > n)
        error("window_size (%d) > number of locations (%d)", window_size, n);
    if (window_size < 2 * margin + 1)
        error("window_size (%d) < 2*margin+1 (%d)", window_size, 2 * margin + 1);
    if ((window_size % 2) == 0 || (margin % 2) == 0)
        error("window_size and margin must both be odd");

    const double *xx  = REAL(r_x);
    const double *xy  = REAL(r_y);
    const double *xtm = REAL(r_time_mins);
    const double *xle = REAL(r_loc_err);

    int n_windows = n - window_size + 1;

    /* Build margin_breaks once: candidate within-window break
     * positions (1-indexed) that satisfy the same constraints the R
     * driver applied:
     *   potential_breaks = 2 .. (window_size - 1)
     *   keep those with margin_breaks >= margin AND
     *                    margin_breaks <= (1 + window_size - margin)
     *                    AND value % 2 == 1
     */
    int n_potential = window_size - 2;
    int *margin_breaks = (int *)R_alloc(n_potential, sizeof(int));
    int n_breaks = 0;
    for (int p = 2; p < window_size; p++) {
        if (p >= margin && p <= (1 + window_size - margin) && (p % 2) == 1) {
            margin_breaks[n_breaks++] = p;
        }
    }
    if (n_breaks < 1)
        error("margin too large for window_size (no valid candidate breaks)");

    int min_brk = margin_breaks[0];
    int max_brk = margin_breaks[0];
    for (int bi = 1; bi < n_breaks; bi++) {
        if (margin_breaks[bi] < min_brk) min_brk = margin_breaks[bi];
        if (margin_breaks[bi] > max_brk) max_brk = margin_breaks[bi];
    }

    /* Output accumulators: sum of variances (for RMS aggregation) and
     * window contribution count, per fix. */
    SEXP r_para_sd = PROTECT(allocVector(REALSXP, n));
    SEXP r_orth_sd = PROTECT(allocVector(REALSXP, n));
    SEXP r_n_estim = PROTECT(allocVector(INTSXP,  n));
    double *para_sd = REAL(r_para_sd);
    double *orth_sd = REAL(r_orth_sd);
    int    *n_estim = INTEGER(r_n_estim);
    double *sum_var_p = (double *)R_alloc(n, sizeof(double));
    double *sum_var_o = (double *)R_alloc(n, sizeof(double));
    for (int i = 0; i < n; i++) {
        sum_var_p[i] = 0.0;
        sum_var_o[i] = 0.0;
        n_estim[i]   = 0;
    }

    int max_pairs = window_size / 2 + 1;

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
        /* Per-thread working buffers (heap, since R_alloc isn't
         * thread-safe).  The cache_* arrays mirror the per-axis
         * breakpoint nll caches the single-window kernel uses to
         * avoid recomputing in step 4. */
        double *Tj  = (double *)malloc(max_pairs * sizeof(double));
        double *al  = (double *)malloc(max_pairs * sizeof(double));
        double *le1 = (double *)malloc(max_pairs * sizeof(double));
        double *le2 = (double *)malloc(max_pairs * sizeof(double));
        double *zp  = (double *)malloc(max_pairs * sizeof(double));
        double *zo  = (double *)malloc(max_pairs * sizeof(double));
        double *time_lag_w = (double *)malloc(window_size * sizeof(double));
        double *cache_fp    = (double *)malloc(n_breaks * sizeof(double));
        double *cache_s2_pb = (double *)malloc(n_breaks * sizeof(double));
        double *cache_s2_pa = (double *)malloc(n_breaks * sizeof(double));
        double *cache_fo    = (double *)malloc(n_breaks * sizeof(double));
        double *cache_s2_ob = (double *)malloc(n_breaks * sizeof(double));
        double *cache_s2_oa = (double *)malloc(n_breaks * sizeof(double));

        #ifdef _OPENMP
        #pragma omp for schedule(static)
        #endif
        for (int w = 0; w < n_windows; w++) {
            const double *wx  = xx  + w;
            const double *wy  = xy  + w;
            const double *wtm = xtm + w;
            const double *wle = xle + w;

            /* time_lag for this window: time_mins[i+1] - time_mins[i].
             * Last entry is unused by precompute_loo_dbgb but written
             * for correctness. */
            for (int i = 0; i < window_size - 1; i++)
                time_lag_w[i] = wtm[i + 1] - wtm[i];
            time_lag_w[window_size - 1] = 0.0;

            int n_pairs = precompute_loo_dbgb(wx, wy, time_lag_w, wle,
                                              window_size,
                                              Tj, al, zp, zo, le1, le2);
            if (n_pairs < 1) continue;

            double logn = log((double)window_size);

            /* Step 1: no break */
            double f_p_full, f_o_full;
            double s2_p_full = brent_min_axis(0.0, 1e15, Tj, al, le1, le2,
                                              zp, n_pairs, 1e-8, &f_p_full);
            double s2_o_full = brent_min_axis(0.0, 1e15, Tj, al, le1, le2,
                                              zo, n_pairs, 1e-8, &f_o_full);
            double bic_none = 2.0 * (f_p_full + f_o_full) + 2.0 * logn;

            double best_bic    = bic_none;
            int    best_pb     = 0, best_ob = 0;
            double best_s2_pb  = s2_p_full, best_s2_pa = s2_p_full;
            double best_s2_ob  = s2_o_full, best_s2_oa = s2_o_full;

            /* Steps 2 & 3: para-only and orth-only breaks */
            for (int bi = 0; bi < n_breaks; bi++) {
                int b_fix = margin_breaks[bi];
                int k = b_fix / 2;
                if (k <= 0 || k >= n_pairs) {
                    cache_fp[bi]    = R_PosInf;
                    cache_fo[bi]    = R_PosInf;
                    cache_s2_pb[bi] = NA_REAL; cache_s2_pa[bi] = NA_REAL;
                    cache_s2_ob[bi] = NA_REAL; cache_s2_oa[bi] = NA_REAL;
                    continue;
                }

                double f_pb, f_pa;
                double s2_pb = brent_min_axis(0.0, 1e15, Tj, al, le1, le2,
                                              zp, k, 1e-8, &f_pb);
                double s2_pa = brent_min_axis(0.0, 1e15,
                                              Tj + k, al + k, le1 + k, le2 + k,
                                              zp + k, n_pairs - k, 1e-8, &f_pa);
                cache_fp[bi]    = f_pb + f_pa;
                cache_s2_pb[bi] = s2_pb; cache_s2_pa[bi] = s2_pa;
                double bic_para_only = 2.0 * (cache_fp[bi] + f_o_full) + 3.0 * logn;
                if (bic_para_only < best_bic) {
                    best_bic = bic_para_only;
                    best_pb  = b_fix; best_ob = 0;
                    best_s2_pb = s2_pb; best_s2_pa = s2_pa;
                    best_s2_ob = s2_o_full; best_s2_oa = s2_o_full;
                }

                double f_ob, f_oa;
                double s2_ob = brent_min_axis(0.0, 1e15, Tj, al, le1, le2,
                                              zo, k, 1e-8, &f_ob);
                double s2_oa = brent_min_axis(0.0, 1e15,
                                              Tj + k, al + k, le1 + k, le2 + k,
                                              zo + k, n_pairs - k, 1e-8, &f_oa);
                cache_fo[bi]    = f_ob + f_oa;
                cache_s2_ob[bi] = s2_ob; cache_s2_oa[bi] = s2_oa;
                double bic_orth_only = 2.0 * (f_p_full + cache_fo[bi]) + 3.0 * logn;
                if (bic_orth_only < best_bic) {
                    best_bic = bic_orth_only;
                    best_pb  = 0; best_ob = b_fix;
                    best_s2_pb = s2_p_full; best_s2_pa = s2_p_full;
                    best_s2_ob = s2_ob; best_s2_oa = s2_oa;
                }
            }

            /* Step 4: cross-conditional refinement */
            if (best_pb != 0) {
                int para_idx = -1;
                for (int bi = 0; bi < n_breaks; bi++)
                    if (margin_breaks[bi] == best_pb) { para_idx = bi; break; }
                if (para_idx >= 0) {
                    double f_para_break = cache_fp[para_idx];
                    for (int bi = 0; bi < n_breaks; bi++) {
                        if (!R_FINITE(cache_fo[bi])) continue;
                        double bic_both = 2.0 * (f_para_break + cache_fo[bi]) + 4.0 * logn;
                        if (bic_both < best_bic) {
                            best_bic = bic_both;
                            best_ob  = margin_breaks[bi];
                            best_s2_ob = cache_s2_ob[bi];
                            best_s2_oa = cache_s2_oa[bi];
                        }
                    }
                }
            }
            if (best_ob != 0) {
                int orth_idx = -1;
                for (int bi = 0; bi < n_breaks; bi++)
                    if (margin_breaks[bi] == best_ob) { orth_idx = bi; break; }
                if (orth_idx >= 0) {
                    double f_orth_break = cache_fo[orth_idx];
                    for (int bi = 0; bi < n_breaks; bi++) {
                        if (!R_FINITE(cache_fp[bi])) continue;
                        double bic_both = 2.0 * (cache_fp[bi] + f_orth_break) + 4.0 * logn;
                        if (bic_both < best_bic) {
                            best_bic = bic_both;
                            best_pb  = margin_breaks[bi];
                            best_s2_pb = cache_s2_pb[bi];
                            best_s2_pa = cache_s2_pa[bi];
                        }
                    }
                }
            }

            /* Accumulate per-fix variance.  Window contributes for
             * fix_pos in [min_brk, max_brk) (1-indexed within window),
             * matching the NA mask in bgb_var_window_c. */
            for (int fix_pos = min_brk; fix_pos < max_brk; fix_pos++) {
                int loc = w + fix_pos - 1;          /* 0-indexed track position */
                double s2p = (best_pb == 0 || fix_pos < best_pb) ? best_s2_pb : best_s2_pa;
                double s2o = (best_ob == 0 || fix_pos < best_ob) ? best_s2_ob : best_s2_oa;
                #ifdef _OPENMP
                #pragma omp atomic
                #endif
                sum_var_p[loc] += s2p;
                #ifdef _OPENMP
                #pragma omp atomic
                #endif
                sum_var_o[loc] += s2o;
                #ifdef _OPENMP
                #pragma omp atomic
                #endif
                n_estim[loc] += 1;
            }
        }

        free(Tj); free(al); free(le1); free(le2); free(zp); free(zo);
        free(time_lag_w);
        free(cache_fp); free(cache_s2_pb); free(cache_s2_pa);
        free(cache_fo); free(cache_s2_ob); free(cache_s2_oa);
    }

    /* RMS aggregation:  sd[i] = sqrt(sum_var[i] / n_estim[i]) */
    for (int i = 0; i < n; i++) {
        if (n_estim[i] > 0) {
            double mean_var_p = sum_var_p[i] / n_estim[i];
            double mean_var_o = sum_var_o[i] / n_estim[i];
            para_sd[i] = (mean_var_p > 0.0) ? sqrt(mean_var_p) : 0.0;
            orth_sd[i] = (mean_var_o > 0.0) ? sqrt(mean_var_o) : 0.0;
        } else {
            para_sd[i] = NA_REAL;
            orth_sd[i] = NA_REAL;
        }
    }

    SEXP result = PROTECT(allocVector(VECSXP, 3));
    SEXP names  = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, mkChar("para_sd"));
    SET_STRING_ELT(names, 1, mkChar("orth_sd"));
    SET_STRING_ELT(names, 2, mkChar("n_estim"));
    setAttrib(result, R_NamesSymbol, names);
    SET_VECTOR_ELT(result, 0, r_para_sd);
    SET_VECTOR_ELT(result, 1, r_orth_sd);
    SET_VECTOR_ELT(result, 2, r_n_estim);

    UNPROTECT(5);
    return result;
}
