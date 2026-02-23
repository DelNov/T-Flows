#include <stdio.h>
#if T_FLOWS_PETSC == 1
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#endif

// C++ sources need this if they are to be linked with Fortran
extern "C" {

/*-----------------------------------------------------------------------------+
|                                                                              |
|  PETSc is defined                                                            |
|                                                                              |
+-----------------------------------------------------------------------------*/
#if T_FLOWS_PETSC == 1
  PetscInt       one  = 1;
  const PetscInt zero = 0;
  PetscErrorCode err;

  static char help[] = "This is to initialize PETSc from T-Flows!\n";

  /*---------------------------------------------------------------------------+
  |  PetscInitialize (and PetscIntialized)                                     |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Sys/PetscInitialize/                |
  |  https://petsc.org/release/manualpages/Sys/PetscInitialized/               |
  +---------------------------------------------------------------------------*/
  void c_petsc_initialize_() {

    PetscBool initialized;
    err = PetscInitialized(&initialized);

    if(!initialized) {
      err = PetscInitialize(0, NULL, (char*)0, help);
    }
  }

  /*---------------------------------------------------------------------------+
  |  PetscLogDefaultBegin                                                      |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Profiling/PetscLogDefaultBegin/     |
  +---------------------------------------------------------------------------*/
  void c_petsc_log_default_begin_() {

    err = PetscLogDefaultBegin();
  }

  /*---------------------------------------------------------------------------+
  |  PetscOptionsSetValue                                                      |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Sys/PetscOptionsSetValue/           |
  +---------------------------------------------------------------------------*/
  void c_petsc_options_set_value_(const char name[], const char value[]) {

    err = PetscOptionsSetValue(NULL, name, value);
  }

  /*---------------------------------------------------------------------------+
  |  PetscLogView                                                              |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Profiling/PetscLogView/             |
  +---------------------------------------------------------------------------*/
  void c_petsc_log_view_() {

    err = PetscLogView(PETSC_VIEWER_STDOUT_WORLD);
  }

  /*---------------------------------------------------------------------------+
  |  PetscFinalize (and PetscInitialized)                                      |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Sys/PetscFinalize/                  |
  |  https://petsc.org/release/manualpages/Sys/PetscFinalized/                 |
  +---------------------------------------------------------------------------*/
  void c_petsc_finalize_() {

    PetscBool initialized;
    err = PetscInitialized(&initialized);

    if(initialized) {
      err = PetscFinalize();
    }
  }

  /*---------------------------------------------------------------------------+
  |                                                                            |
  |  Matrix routines                                                           |
  |                                                                            |
  +---------------------------------------------------------------------------*/

  /*---------------------------------------------------------------------------+
  |  MatCreate                                                                 |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatCreate/                      |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_create_(Mat * A) {

    err = MatCreate(MPI_COMM_WORLD, A);
  }

  /*---------------------------------------------------------------------------+
  |  MatSetSizes                                                               |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatSetSizes/                    |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_set_sizes_(Mat * A, PetscInt * m, PetscInt * M) {

    err = MatSetSizes(*A, *m, *m, *M, *M);
  }

  /*---------------------------------------------------------------------------+
  |  MatSetType (to MATAIJ)                                                    |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatSetType/                     |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_set_type_to_mat_aij_(Mat * A) {

    err = MatSetType(*A, MATAIJ);
  }

  /*---------------------------------------------------------------------------+
  |  MatMPIAIJSetPreallocation                                                 |
  |  MatSeqAIJSetPreallocation                                                 |
  |                                                                            |
  |  Combines two calls which seems to be important or necessary               |
  |  https://petsc.org/release/manualpages/Mat/MATAIJ/                         |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatMPIAIJSetPreallocation/      |
  |  https://petsc.org/release/manualpages/Mat/MatSeqAIJSetPreallocation/      |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_aij_set_preallocation_(Mat      * A,
                                          PetscInt * d_nnz,
                                          PetscInt * o_nnz) {

    err = MatMPIAIJSetPreallocation(*A, 0, d_nnz, 0, o_nnz);
    err = MatSeqAIJSetPreallocation(*A, 0, d_nnz);
  }

  /*---------------------------------------------------------------------------+
  |  MatSetValue                                                               |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatSetValue/                    |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_set_value_(Mat         * A,
                              PetscInt    * row,
                              PetscInt    * col,
                              PetscScalar * value) {

    err = MatSetValue(*A, *row, *col, *value, INSERT_VALUES);
  }

  /*---------------------------------------------------------------------------+
  |  MatSetValues                                                              |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatSetValues/                   |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_set_values_(Mat         * A,
                               PetscInt    * m,
                               PetscInt    * row,
                               PetscInt    * n,
                               PetscInt    * col,
                               PetscScalar * values) {

    err = MatSetValues(*A, *m, row, *n, col, values, INSERT_VALUES);
  }

  /*---------------------------------------------------------------------------+
  |  MatSetValue's sister; this one adds, rather than inserts values           |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatSetValue/                    |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_add_value_(Mat         * A,
                              PetscInt    * row,
                              PetscInt    * col,
                              PetscScalar * value) {

    err = MatSetValue(*A, *row, *col, *value, ADD_VALUES);
  }

  /*---------------------------------------------------------------------------+
  |  MatSetValues' sister; this one adds, rather than inserts values           |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatSetValues/                   |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_add_values_(Mat         * A,
                               PetscInt    * m,
                               PetscInt    * row,
                               PetscInt    * n,
                               PetscInt    * col,
                               PetscScalar * values) {

    err = MatSetValues(*A, *m, row, *n, col, values, ADD_VALUES);
  }

  /*---------------------------------------------------------------------------+
  |  MatAssemblyBegin and MatAssemblyEnd                                       |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatAssemblyBegin/               |
  |  https://petsc.org/release/manualpages/Mat/MatAssemblyEnd/                 |
  |  https://petsc.org/release/manualpages/Mat/MatSetOption/                   |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_assemble_(Mat * A) {

    /* These two always go in pairs */
    err = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
    err = MatAssemblyEnd  (*A, MAT_FINAL_ASSEMBLY);

    /*----------------------------------------------------------------+
    |  Uncomment the block below to check if matrix is symmetric      |
    |  It must be called after MatAssemblyBegin and MatAssemblyEnd    |
    |  Ideally, this should become a separate function in the future  |
    |                                                                 |
    |  https://petsc.org/release/manualpages/Mat/MatIsSymmetric/      |
    +----------------------------------------------------------------*/
    /* PetscBool symmetric;
       MatIsSymmetric(*A, 0.0, &symmetric);
       if(!symmetric) {printf("Matrix is not symmetric\n");}
       else           {printf("Matrix is symmetric\n");} */

    /* If one wishes to repeatedly assemble matrices that retain the
       same nonzero pattern (such as within a time-dependent problem),
       the option should be specified after the first matrix has been
       fully assembled. It didn't lead to any improvement in speed. */
    err = MatSetOption(*A, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
  }

  /*---------------------------------------------------------------------------+
  |  MatSetNullSpace                                                           |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatNullSpaceCreate/             |
  |  https://petsc.org/release/manualpages/Mat/MatSetNullSpace/                |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_set_null_space_(Mat * A) {

    MatNullSpace nullspace;

    MatNullSpaceCreate(MPI_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
    MatSetNullSpace(*A, nullspace);
    MatNullSpaceDestroy(&nullspace);
  }

  /*---------------------------------------------------------------------------+
  |  MatSetNullSpace                                                           |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatSetNullSpace/                |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_remove_null_space_(Mat * A) {

    MatSetNullSpace(*A, NULL);
  }

  /*---------------------------------------------------------------------------+
  |  MatZeroEntries                                                            |
  |                                                                            |
  |  https://petsc.org/main/manualpages/Mat/MatZeroEntries/                    |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_zero_entries_(Mat * A) {

    MatZeroEntries(* A);
  }

  /*---------------------------------------------------------------------------+
  |  MatDestroy                                                                |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatDestroy/                     |
  +---------------------------------------------------------------------------*/
  void c_petsc_mat_destroy_(Mat * A) {

    err = MatDestroy(A);
  }

  /*---------------------------------------------------------------------------+
  |                                                                            |
  |  Vector routines                                                           |
  |                                                                            |
  +---------------------------------------------------------------------------*/

  /*---------------------------------------------------------------------------+
  |  VecCreate                                                                 |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Vec/VecCreate/                      |
  +---------------------------------------------------------------------------*/
  void c_petsc_vec_create_(Vec * v) {

    err = VecCreate(MPI_COMM_WORLD, v);
  }

  /*---------------------------------------------------------------------------+
  |  VecCreateMPI                                                              |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Vec/VecCreateMPI/                   |
  +---------------------------------------------------------------------------*/
  void c_petsc_vec_create_mpi_(Vec * v, PetscInt * m, PetscInt * M) {

    err = VecCreateMPI(MPI_COMM_WORLD, *m, *M, v);
  }

  /*---------------------------------------------------------------------------+
  |  VecSetSizes                                                               |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Vec/VecSetSizes/                    |
  +---------------------------------------------------------------------------*/
  void c_petsc_vec_set_sizes_(Vec * v, PetscInt * m, PetscInt * M) {

    err = VecSetSizes(*v, *m, *M);
  }

  /*---------------------------------------------------------------------------+
  |  VecSetType (to VECSTANDARD)                                               |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Vec/VecSetType/                     |
  +---------------------------------------------------------------------------*/
  void c_petsc_vec_set_type_to_standard_(Vec * v) {

    err = VecSetType(*v, VECSTANDARD);
  }

  /*---------------------------------------------------------------------------+
  |  VecSetValue                                                               |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Vec/VecSetValue/                    |
  +---------------------------------------------------------------------------*/
  void c_petsc_vec_set_value_(Vec * v, PetscInt * row, PetscScalar * value) {

    err = VecSetValue(*v, *row, *value, INSERT_VALUES);
  }

  /*---------------------------------------------------------------------------+
  |  VecSetValue's sister; this one adds, rather than inserts values           |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Vec/VecSetValue/                    |
  +---------------------------------------------------------------------------*/
  void c_petsc_vec_add_value_(Vec * v, PetscInt * row, PetscScalar * value) {

    err = VecSetValue(*v, *row, *value, ADD_VALUES);
  }

  /*---------------------------------------------------------------------------+
  |  VecAssemblyBegin and VecAssemblyEnd                                       |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Vec/VecAssemblyBegin/               |
  |  https://petsc.org/release/manualpages/Vec/VecAssemblyEnd/                 |
  +---------------------------------------------------------------------------*/
  void c_petsc_vec_assemble_(Vec * v) {

    err = VecAssemblyBegin(*v);
    err = VecAssemblyEnd(*v);
  }

  /*---------------------------------------------------------------------------+
  |  VecGetValues                                                              |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Vec/VecGetValues/                   |
  +---------------------------------------------------------------------------*/
  void c_petsc_vec_get_values_(Vec         * v,
                               PetscInt    * ni,
                               PetscInt    * row,
                               PetscScalar * value) {

    err = VecGetValues(*v, *ni, row, value);
  }

  /*---------------------------------------------------------------------------+
  |  VecDestroy                                                                |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Vec/VecDestroy/                     |
  +---------------------------------------------------------------------------*/
  void c_petsc_vec_destroy_(Vec * v) {

    err = VecDestroy(v);
  }

  /*---------------------------------------------------------------------------+
  |                                                                            |
  |  Solver routines                                                           |
  |                                                                            |
  +---------------------------------------------------------------------------*/

  /*---------------------------------------------------------------------------+
  |  KSPCreate                                                                 |
  |                                                                            |
  |  "To solve a linear system with KSP, one must first create a solver con-   |
  |  text with the command KSPCreate" (https://petsc.org/release/manual/ksp/)  |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPCreate/                      |
  +---------------------------------------------------------------------------*/
  void c_petsc_ksp_create_(KSP * ksp) {

    err = KSPCreate(MPI_COMM_WORLD, ksp);
  }

  /*---------------------------------------------------------------------------+
  |  KSPSetOperators                                                           |
  |                                                                            |
  |  "Before actually solving a linear system with KSP, the user must call     |
  |  the following routine to set the matrices associated with the linear      |
  |  system: KSPSetOperators" (https://petsc.org/release/manual/ksp/)          |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPSetOperators/                |
  +---------------------------------------------------------------------------*/
  void c_petsc_ksp_set_operators_(KSP * ksp,
                                  Mat * A) {

    /*--------------------------------------------------------------------+
    |  Here the matrix that defines the linear system also serves as the  |
    |  preconditioning matrix. Since all the matrices will have the same  |
    |  nonzero pattern here, we indicate this so the linear solvers can   |
    |  take advantage of this (How on Earth?*)                            |
    +--------------------------------------------------------------------*/

    /* Set precondioning matrix to be A */
    err = KSPSetOperators(*ksp, *A, *A);
  }

  /*---------------------------------------------------------------------------+
  |  KSPSetType                                                                |
  |                                                                            |
  |  "The Krylov subspace methods accept a number of options, many of which    |
  |  are discussed below. First, to set the Krylov subspace method that is to  |
  |  be used, one calls the command KSPSetType.  The type can be one of:       |
  |  KSPRICHARDSON, KSPCHEBYSHEV, KSPCG, KSPGMRES, KSPTCQMR, KSPBCGS, KSPCGS,  |
  |  KSPTFQMR, KSPCR, KSPLSQR, KSPBICG, KSPPREONLY (or equivalent KSPNONE), or |
  |  others; see KSP Objects or the KSPType man page for more. The KSP method  |
  |  can also be set with the options database command -ksp_type, followed by  |
  |  one of the options: richardson, chebyshev, cg, gmres, tcqmr, bcgs, cgs,   |
  |  tfqmr, cr, lsqr, bicg ..." (https://petsc.org/release/manual/ksp/)        |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPSetType/                     |
  +---------------------------------------------------------------------------*/
  void c_petsc_ksp_set_type_(KSP  * ksp,
                             char * sol) {

    /* Set solver */
    err = KSPSetType(*ksp, sol);
  }

  /*---------------------------------------------------------------------------+
  |  KSP routines to set preconditioner                                        |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPGetPC/                       |
  |  https://petsc.org/release/manualpages/PC/PCSetType/                       |
  |  https://petsc.org/release/manualpages/KSP/KSPSetFromOptions/              |
  |  https://petsc.org/release/manualpages/KSP/KSPSetUp/                       |
  +---------------------------------------------------------------------------*/
  void c_petsc_ksp_set_preconditioner_(KSP  * ksp,
                                       PC   * pc,
                                       char * prec) {
    /* Set preconditioner */
    err = KSPGetPC(*ksp, pc);
    err = PCSetType(*pc, prec);

    /* These two lines are needed to finish the setup */
    err = KSPSetFromOptions(*ksp);
    err = KSPSetUp         (*ksp);
  }

  /*---------------------------------------------------------------------------+
  |  And please don't start from zero - for the sake of God                    |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPSetInitialGuessNonzero/      |v
  +---------------------------------------------------------------------------*/
  void c_petsc_ksp_set_initial_guess_nonzero_(KSP * ksp) {

    err = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE);
  }

  /*---------------------------------------------------------------------------+
  |  KSPSetTolerances                                                          |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPSetTolerances/               |
  +---------------------------------------------------------------------------*/
  void c_petsc_ksp_set_tolerances_(KSP         * ksp,
                                   PetscScalar * rtol,
                                   PetscScalar * atol,
                                   PetscInt    * maxits) {

    err = KSPSetTolerances(*ksp, *rtol, *atol, 1.0e+3, *maxits);
  }

  /*---------------------------------------------------------------------------+
  |  KSPSetTolerances                                                          |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPSolve/                       |
  +---------------------------------------------------------------------------*/
  void c_petsc_ksp_solve_(KSP * ksp, Vec * b, Vec * x) {

    err = KSPSolve(*ksp, *b, *x);
  }

  /*---------------------------------------------------------------------------+
  |  KSPGetConvergedReason                                                     |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPGetConvergedReason/          |
  +---------------------------------------------------------------------------*/
  void c_petsc_ksp_converged_reason_(KSP * ksp, KSPConvergedReason * reason) {

    err = KSPGetConvergedReason(*ksp, reason);
  }

  /*---------------------------------------------------------------------------+
  |  KSPGetIterationNumber                                                     |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPGetIterationNumber/          |
  +---------------------------------------------------------------------------*/
  void c_petsc_ksp_get_iteration_number_(KSP * ksp, PetscInt * its) {

    err = KSPGetIterationNumber(*ksp, its);
  }

  /*---------------------------------------------------------------------------+
  |  KSPGetResidualNorm                                                        |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPGetResidualNorm/             |
  +---------------------------------------------------------------------------*/
  void c_petsc_ksp_get_residual_norm_(KSP * ksp, PetscScalar * rnorm) {

    err = KSPGetResidualNorm(*ksp, rnorm);
  }

  /*---------------------------------------------------------------------------+
  |  KSPDestroy                                                                |
  |                                                                            |
  |  "Once the KSP context is no longer needed, it should be destroyed with    |
  |  the command KSPDestroy" (https://petsc.org/release/manual/ksp/)           |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPDestroy/                     |
  +---------------------------------------------------------------------------*/
  void c_petsc_ksp_destroy_(KSP * ksp) {

    err = KSPDestroy(ksp);
  }

/*-----------------------------------------------------------------------------+
|                                                                              |
|  PETSc is not defined                                                        |
|                                                                              |
+-----------------------------------------------------------------------------*/
#else

  void c_petsc_initialize_() {}

  void c_petsc_finalize_() {}

  void c_petsc_mat_set_null_space_(float * A) {}

  void c_petsc_mat_remove_null_space_(float * A) {}

#endif

} // extern "C"
