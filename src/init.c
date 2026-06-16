#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* Original kernels from move package */
extern SEXP bgb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dbbmm2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP llBGBvar(SEXP, SEXP);

/* C implementation of BM variance estimation */
extern SEXP bm_variance_c(SEXP, SEXP, SEXP, SEXP);
extern SEXP bm_variance_window_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bm_variance_track_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* C implementation of dBGB per-window variance estimation */
extern SEXP bgb_var_window_c(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bgb_var_break_track_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* OpenMP-parallelised grid kernels */
extern SEXP dbbmm2_omp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bgb_omp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"bgb",                  (DL_FUNC) &bgb,                  10},
    {"dbbmm2",               (DL_FUNC) &dbbmm2,               10},
    {"llBGBvar",             (DL_FUNC) &llBGBvar,              2},
    {"bm_variance_c",        (DL_FUNC) &bm_variance_c,         4},
    {"bm_variance_window_c", (DL_FUNC) &bm_variance_window_c,  6},
    {"bm_variance_track_c",  (DL_FUNC) &bm_variance_track_c,   6},
    {"bgb_var_window_c",     (DL_FUNC) &bgb_var_window_c,      5},
    {"bgb_var_break_track_c",(DL_FUNC) &bgb_var_break_track_c, 6},
    {"dbbmm2_omp",           (DL_FUNC) &dbbmm2_omp,           10},
    {"bgb_omp",              (DL_FUNC) &bgb_omp,              10},
    {NULL, NULL, 0}
};

void R_init_move2utils(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
