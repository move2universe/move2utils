/*
 * bm_variance_c.c — C implementation of the BM variance estimation
 *
 * Replaces the pure R .bm_variance() function with a C version that:
 * 1. Pre-computes the leave-one-out quantities (T_jump, alpha, ztz, errors)
 * 2. Uses Brent's method for 1D optimization (avoiding R's optimize() overhead)
 * 3. Processes an entire sliding window's worth of breakpoint tests in one .Call()
 *
 * This is the main performance bottleneck — each window position calls
 * .bm_variance() once for the whole window + once per potential breakpoint.
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

/* Negative log-likelihood for BM variance given pre-computed quantities */
static double neg_log_lik(double var, const double *T_jump, const double *alpha,
                          const double *loc_err_1, const double *loc_err_2,
                          const double *ztz, int n_pairs) {
    double nll = 0.0;
    for (int i = 0; i < n_pairs; i++) {
        double v = T_jump[i] * alpha[i] * (1.0 - alpha[i]) * var +
                   (1.0 - alpha[i]) * (1.0 - alpha[i]) * loc_err_1[i] * loc_err_1[i] +
                   alpha[i] * alpha[i] * loc_err_2[i] * loc_err_2[i];
        if (v <= 0.0) v = DBL_EPSILON;
        nll -= log(1.0 / (2.0 * M_PI * v) * exp(-ztz[i] / (2.0 * v)));
    }
    return nll;
}

/*
 * Brent's method for 1D minimization on [a, b].
 * Returns the minimizer; stores the minimum value in *f_min.
 */
static double brent_min(double a, double b,
                        const double *T_jump, const double *alpha,
                        const double *loc_err_1, const double *loc_err_2,
                        const double *ztz, int n_pairs,
                        double tol, double *f_min) {
    double x, w, v, fx, fw, fv, e, d, u, fu;
    double midpoint, tol1, tol2, p, q, r;
    const double golden = 0.3819660;
    const int max_iter = 500;

    x = w = v = a + golden * (b - a);
    fx = fw = fv = neg_log_lik(x, T_jump, alpha, loc_err_1, loc_err_2, ztz, n_pairs);
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

        /* Try parabolic interpolation */
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
        fu = neg_log_lik(u, T_jump, alpha, loc_err_1, loc_err_2, ztz, n_pairs);

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
 * Pre-compute the leave-one-out quantities for a track segment.
 * Takes x, y, time_lag, location_error arrays of length n.
 * Fills T_jump, alpha, ztz, loc_err_1, loc_err_2 arrays.
 * Returns number of pairs computed.
 */
static int precompute_loo(const double *x, const double *y,
                          const double *time_lag, const double *loc_err,
                          int n,
                          double *T_jump, double *alpha, double *ztz,
                          double *loc_err_1, double *loc_err_2) {
    int n_pairs = 0;
    int i = 1;  /* 0-indexed: start at second element */
    while (i < n - 1) {
        double t = time_lag[i - 1] + time_lag[i];
        T_jump[n_pairs] = t;
        double a = time_lag[i - 1] / t;
        alpha[n_pairs] = a;
        double ux = x[i - 1] + a * (x[i + 1] - x[i - 1]);
        double uy = y[i - 1] + a * (y[i + 1] - y[i - 1]);
        double dx = x[i] - ux;
        double dy = y[i] - uy;
        ztz[n_pairs] = dx * dx + dy * dy;
        loc_err_1[n_pairs] = loc_err[i - 1];
        loc_err_2[n_pairs] = loc_err[i + 1];
        n_pairs++;
        i += 2;
    }
    return n_pairs;
}

/*
 * C entry point: compute BM variance for a single track segment.
 * Called from R as .Call("bm_variance_c", x, y, time_lag, loc_err)
 * Returns list(BMvar, cll)
 */
SEXP bm_variance_c(SEXP x, SEXP y, SEXP time_lag, SEXP loc_err) {
    int n = length(x);
    double *xx = REAL(x);
    double *xy = REAL(y);
    double *xtl = REAL(time_lag);
    double *xle = REAL(loc_err);

    /* Allocate working arrays (max n/2 pairs) */
    int max_pairs = n / 2;
    double *T_jump = (double *)R_alloc(max_pairs, sizeof(double));
    double *alpha_arr = (double *)R_alloc(max_pairs, sizeof(double));
    double *ztz = (double *)R_alloc(max_pairs, sizeof(double));
    double *le1 = (double *)R_alloc(max_pairs, sizeof(double));
    double *le2 = (double *)R_alloc(max_pairs, sizeof(double));

    int n_pairs = precompute_loo(xx, xy, xtl, xle, n,
                                  T_jump, alpha_arr, ztz, le1, le2);

    if (n_pairs == 0) {
        error("Not enough locations for variance estimation");
    }

    double f_min;
    double var_opt = brent_min(0.0, 1e15, T_jump, alpha_arr, le1, le2, ztz,
                                n_pairs, 1e-8, &f_min);

    /* Return list(BMvar=var_opt, cll=-f_min) */
    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("BMvar"));
    SET_STRING_ELT(names, 1, mkChar("cll"));
    setAttrib(result, R_NamesSymbol, names);

    SET_VECTOR_ELT(result, 0, ScalarReal(var_opt));
    SET_VECTOR_ELT(result, 1, ScalarReal(-f_min));

    UNPROTECT(2);
    return result;
}

/*
 * C entry point: process an entire sliding window position.
 * Tests the whole window + all breakpoint positions in one call.
 * Returns list(window_vars, break_pos) where break_pos is 0 if no break.
 */
SEXP bm_variance_window_c(SEXP x, SEXP y, SEXP time_lag, SEXP loc_err,
                           SEXP r_uneven_breaks, SEXP r_breaks_range) {
    int n = length(x);
    double *xx = REAL(x);
    double *xy = REAL(y);
    double *xtl = REAL(time_lag);
    double *xle = REAL(loc_err);
    int *uneven_breaks = INTEGER(r_uneven_breaks);
    int n_breaks = length(r_uneven_breaks);
    int *breaks_range = INTEGER(r_breaks_range);
    int n_range = length(r_breaks_range);

    int max_pairs = n / 2;
    double *T_jump = (double *)R_alloc(max_pairs, sizeof(double));
    double *alpha_arr = (double *)R_alloc(max_pairs, sizeof(double));
    double *ztz = (double *)R_alloc(max_pairs, sizeof(double));
    double *le1 = (double *)R_alloc(max_pairs, sizeof(double));
    double *le2 = (double *)R_alloc(max_pairs, sizeof(double));

    /* Whole window variance */
    int n_pairs = precompute_loo(xx, xy, xtl, xle, n,
                                  T_jump, alpha_arr, ztz, le1, le2);
    double f_whole;
    double var_whole = brent_min(0.0, 1e15, T_jump, alpha_arr, le1, le2, ztz,
                                  n_pairs, 1e-8, &f_whole);
    double bic_whole = 2.0 * f_whole + log((double)n);

    /* Try each breakpoint */
    double best_bic_break = R_PosInf;
    int best_b = 0;
    double best_var_before = var_whole;
    double best_var_after = var_whole;

    for (int bi = 0; bi < n_breaks; bi++) {
        int b = uneven_breaks[bi] - 1;  /* Convert from 1-indexed R to 0-indexed C */

        /* Before section: indices 0..b */
        int np_before = precompute_loo(xx, xy, xtl, xle, b + 1,
                                        T_jump, alpha_arr, ztz, le1, le2);
        if (np_before == 0) continue;
        double f_before;
        double var_before = brent_min(0.0, 1e15, T_jump, alpha_arr, le1, le2, ztz,
                                       np_before, 1e-8, &f_before);

        /* After section: indices b..n-1 */
        int np_after = precompute_loo(xx + b, xy + b, xtl + b, xle + b, n - b,
                                       T_jump, alpha_arr, ztz, le1, le2);
        if (np_after == 0) continue;
        double f_after;
        double var_after = brent_min(0.0, 1e15, T_jump, alpha_arr, le1, le2, ztz,
                                      np_after, 1e-8, &f_after);

        double bic_break = 2.0 * (f_before + f_after) + 2.0 * log((double)n);
        if (bic_break < best_bic_break) {
            best_bic_break = bic_break;
            best_b = uneven_breaks[bi];  /* Keep 1-indexed for R */
            best_var_before = var_before;
            best_var_after = var_after;
        }
    }

    /* Build output: variance per segment position */
    int n_output = n_range - 1;
    SEXP result = PROTECT(allocVector(VECSXP, 3));
    SEXP names = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, mkChar("vars"));
    SET_STRING_ELT(names, 1, mkChar("break_pos"));
    SET_STRING_ELT(names, 2, mkChar("has_break"));
    setAttrib(result, R_NamesSymbol, names);

    SEXP vars = PROTECT(allocVector(REALSXP, n_output));
    double *xvars = REAL(vars);

    if (best_bic_break < bic_whole) {
        /* Break found */
        for (int i = 0; i < n_output; i++) {
            xvars[i] = (breaks_range[i] < best_b) ? best_var_before : best_var_after;
        }
        SET_VECTOR_ELT(result, 1, ScalarInteger(best_b));
        SET_VECTOR_ELT(result, 2, ScalarLogical(1));
    } else {
        /* No break */
        for (int i = 0; i < n_output; i++) {
            xvars[i] = var_whole;
        }
        SET_VECTOR_ELT(result, 1, ScalarInteger(0));
        SET_VECTOR_ELT(result, 2, ScalarLogical(0));
    }
    SET_VECTOR_ELT(result, 0, vars);

    UNPROTECT(3);
    return result;
}

/*
 * Track-level entry point: sweep all sliding windows in C and aggregate the
 * per-fix variance estimates internally.  Replaces the R-level
 *
 *     for (w in 1:n_windows) .Call("bm_variance_window_c", ...)
 *     do.call(rbind, ...) ; aggregate(BMvar ~ loc, ...)
 *
 * pattern with a single .Call.  Optionally parallelized over windows via
 * OpenMP; each window's optimisation is fully independent so the only
 * shared writes are the per-fix accumulator (handled with #pragma omp
 * atomic).  Numerics are identical to bm_variance_window_c per window,
 * verified by tests/testthat/test-dbbmm-variance-dyn.R.
 *
 * Inputs (full track):
 *   x, y, time_lag, loc_err  : length n vectors
 *   window_size              : odd integer
 *   margin                   : odd integer (>=1, <=(window_size-1)/2)
 *
 * Returns list:
 *   variance   : numeric[n]   per-fix mean BM variance, NA outside the
 *                              valid window-coverage region
 *   n_estim    : integer[n]   number of windows contributing to each fix
 *   break_pos  : integer[*]   1-indexed track positions where a window's
 *                              BIC found a breakpoint (deduplicated;
 *                              order not guaranteed under OpenMP)
 */
SEXP bm_variance_track_c(SEXP r_x, SEXP r_y, SEXP r_time_lag, SEXP r_loc_err,
                          SEXP r_window_size, SEXP r_margin) {
    int n = length(r_x);
    int window_size = asInteger(r_window_size);
    int margin = asInteger(r_margin);

    if (window_size > n)
        error("window_size (%d) > number of locations (%d)", window_size, n);
    if (window_size < 2 * margin + 1)
        error("window_size (%d) < 2*margin+1 (%d)", window_size, 2 * margin + 1);
    if ((window_size % 2) == 0 || (margin % 2) == 0)
        error("window_size and margin must both be odd");

    const double *xx  = REAL(r_x);
    const double *xy  = REAL(r_y);
    const double *xtl = REAL(r_time_lag);
    const double *xle = REAL(r_loc_err);

    int n_windows = n - window_size + 1;

    /* breaks_range: margin..(window_size - margin + 1)  (1-indexed positions
     * inside a window).  uneven_breaks: subset whose value is odd.
     * Both stable across windows; precompute once. */
    int n_range  = window_size - 2 * margin + 2;
    int n_output = n_range - 1;

    int *uneven_breaks = (int *)R_alloc(n_range, sizeof(int));
    int n_uneven = 0;
    for (int i = 0; i < n_range; i++) {
        int v = margin + i;
        if ((v % 2) == 1) uneven_breaks[n_uneven++] = v;
    }

    /* Allocate output accumulators inside R-level SEXPs.  Initialise to
     * 0 / NA; the OpenMP region writes via atomic updates. */
    SEXP r_var      = PROTECT(allocVector(REALSXP, n));
    SEXP r_n_estim  = PROTECT(allocVector(INTSXP,  n));
    double *variance = REAL(r_var);
    int    *n_estim  = INTEGER(r_n_estim);
    for (int i = 0; i < n; i++) { variance[i] = 0.0; n_estim[i] = 0; }

    /* Per-window break flags allocated outside the parallel region; each
     * window writes its own slot (no race) and we collect afterwards. */
    int *break_per_win = (int *)R_alloc(n_windows, sizeof(int));
    for (int i = 0; i < n_windows; i++) break_per_win[i] = 0;

    int max_pairs = window_size / 2 + 1;

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
        /* Per-thread working buffers.  malloc/free, NOT R_alloc, since
         * R_alloc isn't thread-safe.  These are small (5 * max_pairs *
         * sizeof(double) ~= 1 KB for window_size = 31). */
        double *T_jump    = (double *)malloc(max_pairs * sizeof(double));
        double *alpha_arr = (double *)malloc(max_pairs * sizeof(double));
        double *ztz       = (double *)malloc(max_pairs * sizeof(double));
        double *le1       = (double *)malloc(max_pairs * sizeof(double));
        double *le2       = (double *)malloc(max_pairs * sizeof(double));

        #ifdef _OPENMP
        #pragma omp for schedule(static)
        #endif
        for (int w = 0; w < n_windows; w++) {
            const double *wx  = xx  + w;
            const double *wy  = xy  + w;
            const double *wtl = xtl + w;
            const double *wle = xle + w;
            int wn = window_size;

            /* Whole-window variance */
            int n_pairs = precompute_loo(wx, wy, wtl, wle, wn,
                                         T_jump, alpha_arr, ztz, le1, le2);
            if (n_pairs == 0) continue;
            double f_whole;
            double var_whole = brent_min(0.0, 1e15, T_jump, alpha_arr, le1, le2, ztz,
                                         n_pairs, 1e-8, &f_whole);
            double bic_whole = 2.0 * f_whole + log((double)wn);

            /* Try each candidate breakpoint */
            double best_bic = R_PosInf;
            int    best_b   = 0;
            double var_before = var_whole, var_after = var_whole;
            for (int bi = 0; bi < n_uneven; bi++) {
                int b = uneven_breaks[bi] - 1;
                int np_b = precompute_loo(wx, wy, wtl, wle, b + 1,
                                          T_jump, alpha_arr, ztz, le1, le2);
                if (np_b == 0) continue;
                double f_b;
                double v_b = brent_min(0.0, 1e15, T_jump, alpha_arr, le1, le2, ztz,
                                       np_b, 1e-8, &f_b);
                int np_a = precompute_loo(wx + b, wy + b, wtl + b, wle + b, wn - b,
                                          T_jump, alpha_arr, ztz, le1, le2);
                if (np_a == 0) continue;
                double f_a;
                double v_a = brent_min(0.0, 1e15, T_jump, alpha_arr, le1, le2, ztz,
                                       np_a, 1e-8, &f_a);
                double bic = 2.0 * (f_b + f_a) + 2.0 * log((double)wn);
                if (bic < best_bic) {
                    best_bic   = bic;
                    best_b     = uneven_breaks[bi];
                    var_before = v_b;
                    var_after  = v_a;
                }
            }

            int has_break = (best_bic < bic_whole);
            if (has_break) break_per_win[w] = w + best_b;  /* 1-indexed track position */

            /* Accumulate per-fix variance.  Each window contributes to
             * positions w + (margin..(window_size - margin)) (0-indexed). */
            for (int i = 0; i < n_output; i++) {
                int loc = w + margin + i - 1;  /* 0-indexed */
                double v = has_break
                    ? ((margin + i) < best_b ? var_before : var_after)
                    : var_whole;
                #ifdef _OPENMP
                #pragma omp atomic
                #endif
                variance[loc] += v;
                #ifdef _OPENMP
                #pragma omp atomic
                #endif
                n_estim[loc] += 1;
            }
        }

        free(T_jump); free(alpha_arr); free(ztz); free(le1); free(le2);
    }

    /* Convert accumulators to means; mark uncovered fixes NA.
     * Aggregate the break-per-window vector into a deduplicated 1-indexed
     * list (matches the legacy R-level `breaks_found` semantics). */
    for (int i = 0; i < n; i++) {
        variance[i] = (n_estim[i] > 0) ? variance[i] / n_estim[i] : NA_REAL;
    }

    int n_breaks = 0;
    for (int w = 0; w < n_windows; w++) if (break_per_win[w] > 0) n_breaks++;
    SEXP r_breaks = PROTECT(allocVector(INTSXP, n_breaks));
    int *bp = INTEGER(r_breaks);
    int j = 0;
    for (int w = 0; w < n_windows; w++) {
        if (break_per_win[w] > 0) bp[j++] = break_per_win[w];
    }

    SEXP result = PROTECT(allocVector(VECSXP, 3));
    SEXP names  = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, mkChar("variance"));
    SET_STRING_ELT(names, 1, mkChar("n_estim"));
    SET_STRING_ELT(names, 2, mkChar("break_pos"));
    setAttrib(result, R_NamesSymbol, names);
    SET_VECTOR_ELT(result, 0, r_var);
    SET_VECTOR_ELT(result, 1, r_n_estim);
    SET_VECTOR_ELT(result, 2, r_breaks);

    UNPROTECT(5);
    return result;
}
