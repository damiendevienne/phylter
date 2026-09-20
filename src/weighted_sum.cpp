#include <R.h>
#include <Rinternals.h>

// Internal kernel: retain the original multiply-then-add order, and never
// mutate an input matrix (R callers may share these with the initial state).
extern "C" SEXP phylter_weighted_sum(SEXP matrices, SEXP weights) {
    const R_xlen_t genes = XLENGTH(matrices);
    if (TYPEOF(matrices) != VECSXP || genes == 0 ||
        TYPEOF(weights) != REALSXP || XLENGTH(weights) != genes)
        Rf_error("Invalid matrices or weights for compromise");
    SEXP first = VECTOR_ELT(matrices, 0);
    const R_xlen_t size = XLENGTH(first);
    for (R_xlen_t g = 0; g < genes; ++g) {
        SEXP mat = VECTOR_ELT(matrices, g);
        if (TYPEOF(mat) != REALSXP || XLENGTH(mat) != size)
            Rf_error("Compromise matrices must be double matrices of equal size");
    }
    SEXP result = PROTECT(Rf_allocVector(REALSXP, size));
    Rf_setAttrib(result, R_DimSymbol, Rf_getAttrib(first, R_DimSymbol));
    double* out = REAL(result);
    const double* w = REAL(weights);
    const double* initial = REAL(first);
    for (R_xlen_t j = 0; j < size; ++j) out[j] = initial[j] * w[0];
    for (R_xlen_t g = 1; g < genes; ++g) {
        R_CheckUserInterrupt();
        const double* x = REAL(VECTOR_ELT(matrices, g));
        for (R_xlen_t j = 0; j < size; ++j) out[j] += x[j] * w[g];
    }
    UNPROTECT(1);
    return result;
}
