/*=============================================================================+
| Shared AMGX constants used by both C/C++ and Fortran (via cpp/fpp)           |
+-----------------------------------------------------------------------------*/

/* Verbosity */
#define AMGX_QUIET  0
#define AMGX_LOUD   1

/* Call semantics */
#define AMGX_FIRST_CALL   0
#define AMGX_FINAL_CALL  -1

/* Matrix update flags */
#define AMGX_MATRIX_SAME     0
#define AMGX_MATRIX_CHANGED  1

