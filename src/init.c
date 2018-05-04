#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP EstFuncC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MatchLoopC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
		       SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MatchLoopCfast(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
			   SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP VarCalcMatchC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
			  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FastMatchC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FasterMatchC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
  {"EstFuncC", (DL_FUNC) &EstFuncC, 7},
  {"MatchLoopC", (DL_FUNC) &MatchLoopC, 18},
  {"MatchLoopCfast", (DL_FUNC) &MatchLoopC, 17},
  {"VarCalcMatch", (DL_FUNC) &MatchLoopC, 17},
  {"FastMatchC", (DL_FUNC) &MatchLoopC, 9},
  {"FasterMatchC", (DL_FUNC) &MatchLoopC, 8},
  {NULL, NULL, 0}
};

void R_init_EstFuncC(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

void R_init_MatchLoopC(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

void R_init_MatchLoopCfast(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

void R_init_VarCalcMatch(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

void R_init_FastMatchC(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

void R_init_FasterMatchC(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
