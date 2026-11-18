/*
 * ud_grid_probe.c -- the grid bounds of the dBBMM and dBGB UD kernels, without
 * the densities.
 *
 * dbbmm2_omp() and bgb_omp() (bgb_bbmm_omp.c, ported from move) spread every
 * bridge time step over a window of `ext` standard deviations around the bridge
 * mean, and stop with an error when that window leaves the grid.  This routine
 * walks exactly the same time steps with exactly the same arithmetic, so the R
 * side can decide, before calling a kernel, whether a grid is large enough and,
 * if it is not, how large it has to be.  The kernels themselves are unchanged.
 *
 * The two kernels test `xEnd > nXGrid` and `yEnd > nYGrid` but loop to `<=`, so
 * at equality they read past the grid-coordinate vectors and write outside the
 * grid.  The y window is also one row lower than the x logic would suggest:
 * rows floor(v)-1 .. ceil(v')-1.  Both are inherited from move and kept; the R
 * side treats a grid as safe only when the kernel never reaches those indices.
 *
 * Arguments mirror the kernels' own: mode 0 = dBBMM (s1 = per-fix variance,
 * `interest` gates segments), mode 1 = dBGB (s1, s2 = parallel and orthogonal
 * standard deviations; the caller has subset to the core fixes, and
 * `interest` gates segments as for dBBMM).
 *
 * Returns a double vector of length 8:
 *   [0..3] the window envelope in map units: min(mux - w), max(mux + w),
 *          min(muy - w), max(muy + w), with w the kernel's window half-width;
 *          +/-Inf when no time step is evaluated;
 *   [4..7] the index extremes on the given grid: min xStart, max xEnd,
 *          max yEnd, min yStart (as the kernels compute them).
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <limits.h>

SEXP ud_grid_probe(SEXP mode_, SEXP x, SEXP y, SEXP s1, SEXP s2, SEXP t,
                   SEXP locEr, SEXP xGrid, SEXP yGrid, SEXP tStep, SEXP ext2,
                   SEXP interest)
{
    int mode = asInteger(mode_);
    PROTECT(x = coerceVector(x, REALSXP));
    PROTECT(y = coerceVector(y, REALSXP));
    PROTECT(s1 = coerceVector(s1, REALSXP));
    PROTECT(s2 = coerceVector(s2, REALSXP));
    PROTECT(t = coerceVector(t, REALSXP));
    PROTECT(locEr = coerceVector(locEr, REALSXP));
    PROTECT(xGrid = coerceVector(xGrid, REALSXP));
    PROTECT(yGrid = coerceVector(yGrid, REALSXP));
    PROTECT(tStep = coerceVector(tStep, REALSXP));
    PROTECT(ext2 = coerceVector(ext2, REALSXP));
    PROTECT(interest = coerceVector(interest, LGLSXP));

    double ext = REAL(ext2)[0];
    double dT = REAL(tStep)[0];
    R_len_t nLoc = length(x);
    R_len_t nXGrid = length(xGrid);
    R_len_t nYGrid = length(yGrid);

    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, 8));
    double *a = REAL(ans);
    a[0] = R_PosInf; a[1] = R_NegInf; a[2] = R_PosInf; a[3] = R_NegInf;
    a[4] = INT_MAX;  a[5] = INT_MIN;  a[6] = INT_MIN;  a[7] = INT_MAX;

    if (nLoc < 2 || nXGrid < 2 || nYGrid < 2 || !(dT > 0)) {
        UNPROTECT(12);
        return ans;
    }

    double *xt = REAL(t);
    double *xlocEr = REAL(locEr);
    double *xx = REAL(x);
    double *xy = REAL(y);
    double *xxGrid = REAL(xGrid);
    double *xyGrid = REAL(yGrid);
    double *xs1 = REAL(s1);
    double *xs2 = REAL(s2);
    int *intr = LOGICAL(interest);
    double x0 = xxGrid[0];
    double y0 = xyGrid[0];
    double xRes = xxGrid[1] - xxGrid[0];
    double yRes = xyGrid[1] - xyGrid[0];

    double exmin = R_PosInf, exmax = R_NegInf, eymin = R_PosInf, eymax = R_NegInf;
    int minXS = INT_MAX, maxXE = INT_MIN, maxYE = INT_MIN, minYS = INT_MAX;

    R_len_t k = 0;
    int interrupt_counter = 0;

    /* Same loop header and segment walk as both kernels. */
    for (double ti = xt[0] + (fmod((xt[nLoc-1] - xt[0]), dT) / 2.0);
         ti <= xt[nLoc-1]; ti += dT)
    {
        while (xt[k+1] < ti) {
            k++;
        }

        if (++interrupt_counter % 10000 == 0) {
            R_CheckUserInterrupt();
        }

        double alpha = (ti - xt[k]) / (xt[k+1] - xt[k]);
        double mux = xx[k] + (xx[k+1] - xx[k]) * alpha;
        double muy = xy[k] + (xy[k+1] - xy[k]) * alpha;
        double w;
        int xStart, xEnd, yEnd, yStart;

        if (mode == 0) {
            /* dbbmm2_omp */
            if (!intr[k]) continue;
            double sigma = (xt[k+1] - xt[k]) * alpha * (1 - alpha) * xs1[k] +
                            pow(1 - alpha, 2) * pow(xlocEr[k], 2) +
                            pow(alpha, 2) * pow(xlocEr[k+1], 2);
            double sigma_ext = sqrt(sigma) * ext;
            w = sigma_ext;
            xStart = (int)floor((mux - x0) / xRes - sigma_ext / xRes);
            xEnd   = (int)ceil((mux - x0) / xRes + sigma_ext / xRes);
            yEnd   = (int)nYGrid - floor((muy - y0) / yRes - sigma_ext / yRes);
            yStart = (int)nYGrid - ceil((muy - y0) / yRes + sigma_ext / yRes);
        } else {
            /* bgb_omp */
            if (!intr[k]) continue;
            double sigmaPara = sqrt((xt[k+1] - xt[k]) * alpha * (1 - alpha) * pow(xs1[k], 2) +
                                    pow(1 - alpha, 2) * pow(xlocEr[k], 2) +
                                    pow(alpha, 2) * pow(xlocEr[k+1], 2));
            double sigmaOrth = sqrt((xt[k+1] - xt[k]) * alpha * (1 - alpha) * pow(xs2[k], 2) +
                                    pow(1 - alpha, 2) * pow(xlocEr[k], 2) +
                                    pow(alpha, 2) * pow(xlocEr[k+1], 2));
            double maxSigma = fmax2(sigmaPara, sigmaOrth);
            w = maxSigma * ext;
            xStart = (int)floor((mux - x0) / xRes - maxSigma * ext / xRes);
            xEnd   = (int)ceil((mux - x0) / xRes + maxSigma * ext / xRes);
            yEnd   = (int)nYGrid - floor((muy - y0) / yRes - maxSigma * ext / yRes);
            yStart = (int)nYGrid - ceil((muy - y0) / yRes + maxSigma * ext / yRes);
        }

        if (mux - w < exmin) exmin = mux - w;
        if (mux + w > exmax) exmax = mux + w;
        if (muy - w < eymin) eymin = muy - w;
        if (muy + w > eymax) eymax = muy + w;
        if (xStart < minXS) minXS = xStart;
        if (xEnd > maxXE) maxXE = xEnd;
        if (yEnd > maxYE) maxYE = yEnd;
        if (yStart < minYS) minYS = yStart;
    }

    a[0] = exmin; a[1] = exmax; a[2] = eymin; a[3] = eymax;
    a[4] = minXS; a[5] = maxXE; a[6] = maxYE; a[7] = minYS;
    UNPROTECT(12);
    return ans;
}
