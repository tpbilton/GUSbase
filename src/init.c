#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "GUSbase.h"

/* .Call calls */
extern SEXP pest_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pest_em_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gest_em_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pest_em_ep_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"pest_c", (DL_FUNC) &pest_c, 7},
    {"pest_em_c", (DL_FUNC) &pest_em_c, 9},
    {"gest_em_c", (DL_FUNC) &gest_em_c, 9},
    {"pest_em_ep_c", (DL_FUNC) &pest_em_ep_c, 8},
    {NULL, NULL, 0}
};

void R_init_GUSbase(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
