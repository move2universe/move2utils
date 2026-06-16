/*
 * bgb_bbmm_omp.c — OpenMP-parallelised versions of dbbmm2 and bgb.
 *
 * These are drop-in replacements for the original functions in bgb_bbmm.c.
 * The outer time loop remains sequential (segment index k depends on
 * previous iteration). The inner grid loops are parallelised with OpenMP.
 *
 * Compiles and runs correctly without OpenMP — the pragmas are guarded
 * with #ifdef _OPENMP. Without OpenMP, behaviour is identical to the
 * original functions.
 */

#include <R.h>
#include <math.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#ifdef _OPENMP
#include <omp.h>
#endif

SEXP dbbmm2_omp(SEXP x, SEXP y, SEXP s, SEXP t, SEXP locEr,
                 SEXP xGrid, SEXP yGrid, SEXP tStep, SEXP ext2, SEXP interest)
{
    PROTECT(ext2 = coerceVector(ext2, REALSXP));
    PROTECT(interest = coerceVector(interest, LGLSXP));
    PROTECT(xGrid = coerceVector(xGrid, REALSXP));
    PROTECT(yGrid = coerceVector(yGrid, REALSXP));
    PROTECT(x = coerceVector(x, REALSXP));
    PROTECT(y = coerceVector(y, REALSXP));
    PROTECT(s = coerceVector(s, REALSXP));
    PROTECT(t = coerceVector(t, REALSXP));
    PROTECT(tStep = coerceVector(tStep, REALSXP));
    PROTECT(locEr = coerceVector(locEr, REALSXP));

    double ext = REAL(ext2)[0];
    double dT = REAL(tStep)[0];
    R_len_t nLoc = length(x);
    R_len_t nXGrid = length(xGrid);
    R_len_t nYGrid = length(yGrid);

    double x0 = REAL(xGrid)[0];
    double y0 = REAL(yGrid)[0];
    double *xt = REAL(t);
    double *xlocEr = REAL(locEr);
    double *xx = REAL(x);
    double *xy = REAL(y);
    double *xyGrid = REAL(yGrid);
    double *xxGrid = REAL(xGrid);
    double xRes = xxGrid[1] - xxGrid[0];
    double yRes = xyGrid[1] - xyGrid[0];
    double *xs = REAL(s);

    SEXP ans;
    PROTECT(ans = allocMatrix(REALSXP, nYGrid, nXGrid));
    double *rans = REAL(ans);

    /* Zero the output matrix */
    R_len_t total_cells = nXGrid * nYGrid;
    for (R_len_t c = 0; c < total_cells; c++) {
        rans[c] = 0.0;
    }

    R_len_t k = 0;
    int interrupt_counter = 0;

    for (double ti = xt[0] + (fmod((xt[nLoc-1] - xt[0]), dT) / 2.0);
         ti <= xt[nLoc-1]; ti += dT)
    {
        while (xt[k+1] < ti) {
            k++;
        }

        if (!LOGICAL(interest)[k]) continue;

        double alpha = (ti - xt[k]) / (xt[k+1] - xt[k]);
        double mux = xx[k] + (xx[k+1] - xx[k]) * alpha;
        double muy = xy[k] + (xy[k+1] - xy[k]) * alpha;
        double sigma = (xt[k+1] - xt[k]) * alpha * (1 - alpha) * xs[k] +
                        pow(1 - alpha, 2) * pow(xlocEr[k], 2) +
                        pow(alpha, 2) * pow(xlocEr[k+1], 2);

        double sigma_ext = sqrt(sigma) * ext;
        int xStart = (int)floor((mux - x0) / xRes - sigma_ext / xRes);
        int xEnd   = (int)ceil((mux - x0) / xRes + sigma_ext / xRes);
        int yEnd   = (int)nYGrid - floor((muy - y0) / yRes - sigma_ext / yRes);
        int yStart = (int)nYGrid - ceil((muy - y0) / yRes + sigma_ext / yRes);

        if (xStart < 0)
            error("Lower x grid not large enough, consider extending the raster in that direction or enlarging the ext argument");
        if (xEnd > nXGrid)
            error("Higher x grid not large enough, consider extending the raster in that direction or enlarging the ext argument");
        if (yEnd > nYGrid)
            error("Lower y grid not large enough, consider extending the raster in that direction or enlarging the ext argument");
        if (yStart < 0)
            error("Higher y grid not large enough, consider extending the raster in that direction or enlarging the ext argument");

        /* Check for user interrupt every 100 time steps */
        if (++interrupt_counter % 100 == 0) {
            R_CheckUserInterrupt();
        }

        double inv_2pi_sigma = 1.0 / (2.0 * M_PI * sigma);
        double inv_2sigma = -1.0 / (2.0 * sigma);
        double cell_area = xRes * yRes;

        /* OpenMP parallelise the inner grid loop */
        #ifdef _OPENMP
        #pragma omp parallel for collapse(2) schedule(static)
        #endif
        for (int ii = xStart; ii <= xEnd; ii++) {
            for (int jj = yStart; jj <= yEnd; jj++) {
                double dx = xxGrid[ii] - mux;
                double dy = xyGrid[nYGrid - jj - 1] - muy;
                double ZTZ = dx * dx + dy * dy;
                rans[ii * nYGrid + jj] += inv_2pi_sigma * exp(inv_2sigma * ZTZ) * cell_area;
            }
        }
    }

    UNPROTECT(11);
    return ans;
}


SEXP bgb_omp(SEXP x, SEXP y, SEXP sPara, SEXP sOrth, SEXP t, SEXP locEr,
              SEXP xGrid, SEXP yGrid, SEXP dTT, SEXP ext2)
{
    PROTECT(ext2 = coerceVector(ext2, REALSXP));
    PROTECT(dTT = coerceVector(dTT, REALSXP));
    R_len_t nLoc = length(x);
    R_len_t nXGrid = length(xGrid);
    R_len_t nYGrid = length(yGrid);
    PROTECT(xGrid = coerceVector(xGrid, REALSXP));
    PROTECT(yGrid = coerceVector(yGrid, REALSXP));
    PROTECT(x = coerceVector(x, REALSXP));
    PROTECT(y = coerceVector(y, REALSXP));
    PROTECT(sPara = coerceVector(sPara, REALSXP));
    PROTECT(sOrth = coerceVector(sOrth, REALSXP));
    PROTECT(t = coerceVector(t, REALSXP));
    PROTECT(locEr = coerceVector(locEr, REALSXP));

    double ext = REAL(ext2)[0];
    double dT = REAL(dTT)[0];
    double x0 = REAL(xGrid)[0];
    double y0 = REAL(yGrid)[0];
    double *xt = REAL(t);
    double *xlocEr = REAL(locEr);
    double *xx = REAL(x);
    double *xy = REAL(y);
    double *xyGrid = REAL(yGrid);
    double *xxGrid = REAL(xGrid);
    double xRes = xxGrid[1] - xxGrid[0];
    double yRes = xyGrid[1] - xyGrid[0];
    double *xsPara = REAL(sPara);
    double *xsOrth = REAL(sOrth);

    SEXP ans;
    PROTECT(ans = allocMatrix(REALSXP, nYGrid, nXGrid));
    double *rans = REAL(ans);

    R_len_t total_cells = nXGrid * nYGrid;
    for (R_len_t c = 0; c < total_cells; c++) {
        rans[c] = 0.0;
    }

    R_len_t k = 0;
    int interrupt_counter = 0;

    for (double ti = xt[0] + (fmod((xt[nLoc-1] - xt[0]), dT) / 2.0);
         ti <= xt[nLoc-1]; ti += dT)
    {
        while (xt[k+1] < ti) {
            k++;
        }

        /* Check for user interrupt every 100 time steps */
        if (++interrupt_counter % 100 == 0) {
            R_CheckUserInterrupt();
        }

        double alpha = (ti - xt[k]) / (xt[k+1] - xt[k]);
        double mux = xx[k] + (xx[k+1] - xx[k]) * alpha;
        double muy = xy[k] + (xy[k+1] - xy[k]) * alpha;
        double sigmaPara = sqrt((xt[k+1] - xt[k]) * alpha * (1 - alpha) * pow(xsPara[k], 2) +
                                pow(1 - alpha, 2) * pow(xlocEr[k], 2) +
                                pow(alpha, 2) * pow(xlocEr[k+1], 2));
        double sigmaOrth = sqrt((xt[k+1] - xt[k]) * alpha * (1 - alpha) * pow(xsOrth[k], 2) +
                                pow(1 - alpha, 2) * pow(xlocEr[k], 2) +
                                pow(alpha, 2) * pow(xlocEr[k+1], 2));

        double maxSigma = fmax2(sigmaPara, sigmaOrth);
        int xStart = (int)floor((mux - x0) / xRes - maxSigma * ext / xRes);
        int xEnd   = (int)ceil((mux - x0) / xRes + maxSigma * ext / xRes);
        int yEnd   = (int)nYGrid - floor((muy - y0) / yRes - maxSigma * ext / yRes);
        int yStart = (int)nYGrid - ceil((muy - y0) / yRes + maxSigma * ext / yRes);

        if (xStart < 0)
            error("The raster does not extent far enough in the X dimension towards the left side");
        if (xEnd > nXGrid)
            error("The raster does not extent far enough in the X dimension towards the right side");
        if (yEnd > nYGrid)
            error("The raster does not extent far enough in the Y dimension towards the lower side");
        if (yStart < 0)
            error("The raster does not extent far enough in the Y dimension towards the upper side");

        /* Pre-compute direction vector C (constant for this time step) */
        double C = sqrt(pow(mux - xx[k+1], 2) + pow(muy - xy[k+1], 2));
        double inv_2pi_sp_so = 1.0 / (2.0 * M_PI * sigmaPara * sigmaOrth);
        double inv_sp2 = -0.5 / (sigmaPara * sigmaPara);
        double inv_so2 = -0.5 / (sigmaOrth * sigmaOrth);
        double cell_area = xRes * yRes;

        /* OpenMP parallelise the inner grid loop */
        #ifdef _OPENMP
        #pragma omp parallel for collapse(2) schedule(static)
        #endif
        for (int ii = xStart; ii <= xEnd; ii++) {
            for (int jj = yStart; jj <= yEnd; jj++) {
                double gx = xxGrid[ii];
                double gy = xyGrid[nYGrid - jj - 1];
                double B = sqrt(pow(gx - xx[k+1], 2) + pow(gy - xy[k+1], 2));
                double A = sqrt(pow(mux - gx, 2) + pow(muy - gy, 2));

                double deltaPara, deltaOrth;
                if (A == 0.0) {
                    deltaPara = 0.0;
                    deltaOrth = 0.0;
                } else if (C == 0.0) {
                    double s = sqrt(0.5);
                    deltaOrth = A * sin(acos(s));
                    deltaPara = A * cos(acos(s));
                } else {
                    double costheta = (B*B - A*A - C*C) / (2.0 * A * C);
                    if (costheta > 1.0) costheta = 1.0;
                    if (costheta < -1.0) costheta = -1.0;
                    double theta = acos(costheta);
                    deltaOrth = A * sin(theta);
                    deltaPara = A * cos(theta);
                }

                double val = inv_2pi_sp_so *
                    exp(inv_sp2 * deltaPara * deltaPara + inv_so2 * deltaOrth * deltaOrth) *
                    cell_area;
                rans[ii * nYGrid + jj] += val;
            }
        }
    }

    UNPROTECT(11);
    return ans;
}
