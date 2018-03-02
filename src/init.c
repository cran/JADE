#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void FG(void *, void *, void *, void *, void *, void *);
extern void rjdc(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"FG",   (DL_FUNC) &FG,   6},
    {"rjdc", (DL_FUNC) &rjdc, 5},
    {NULL, NULL, 0}
};

void R_init_JADE(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
