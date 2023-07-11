// extract from the package mrfDepth - 07/2023

#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void medcoupleC(void *, void *, void *, void *);

/* .Call calls */
extern SEXP dirOutl_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"medcoupleC",       (DL_FUNC) &medcoupleC,        4},
    {NULL, NULL, 0}
};

void R_init_phylter(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
