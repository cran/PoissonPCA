#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

//Generated automatically by R tools::package_native_routine_registration_skeleton

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void get_scores_log(void *, void *, void *, void *, void *, void *, void *, void *);
extern void log_ECVar(void *, void *, void *, void *, void *);
extern void log_ECVar_unbiassed(void *, void *, void *, void *, void *, void *, void *);
extern void log_g(void *, void *, void *, void *);
extern void precalculate(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"get_scores_log_cpp",      (DL_FUNC) &get_scores_log,      8},
    {"log_ECVar",           (DL_FUNC) &log_ECVar,           5},
    {"log_ECVar_unbiassed", (DL_FUNC) &log_ECVar_unbiassed, 7},
    {"log_g",               (DL_FUNC) &log_g,               4},
    {"precalculate",        (DL_FUNC) &precalculate,        4},
    {NULL, NULL, 0}
};

void R_init_PoissonPCA(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
