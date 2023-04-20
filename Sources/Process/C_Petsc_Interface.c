#include <stdio.h>
#if T_FLOWS_PETSC == 1
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#endif

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
  |  PetscInitialize (and PetcIntialized)                                      |
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
  |  PetscOptionsSetValue                                                      |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Sys/PetscOptionsSetValue/           |
  +---------------------------------------------------------------------------*/
  void c_petsc_options_value_(const char name[], const char value[]) {

    err = PetscOptionsSetValue(NULL, name, value);
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
  |  MatAssemblyBegin and MatAssemblyEnd                                       |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Mat/MatAssemblyBegin/               |
  |  https://petsc.org/release/manualpages/Mat/MatAssemblyEnd/                 |
  |  https://petsc.org/release/manualpages/Mat/MatSetOption/                   |
  +---------------------------------------------------------------------------*/
  void c_petsc_assemble_mat_(Mat * A) {

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
  |  VecAssemblyBegin and VecAssemblyEnd                                       |
  |                                                                            |
  |  https://petsc.org/release/manualpages/Vec/VecAssemblyBegin/               |
  |  https://petsc.org/release/manualpages/Vec/VecAssemblyEnd/                 |
  +---------------------------------------------------------------------------*/
  void c_petsc_assemble_vec_(Vec * v) {

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
  |                                                                            |
  |  Solver routines                                                           |
  |                                                                            |
  +---------------------------------------------------------------------------*/

  /*---------------------------------------------------------------------------+
  |  KSPCreate                                                                 |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPCreate/                      |
  +---------------------------------------------------------------------------*/
  void c_petsc_ksp_create_(KSP * ksp) {

    err = KSPCreate(MPI_COMM_WORLD, ksp);
  }

  /*---------------------------------------------------------------------------+
  |  KSP routines to set solver and preconditioner                             |
  |                                                                            |
  |  https://petsc.org/release/manualpages/KSP/KSPSetOperators/                |
  |  https://petsc.org/release/manualpages/KSP/KSPSetType/                     |
  |  https://petsc.org/release/manualpages/KSP/KSPGetPC/                       |
  |  https://petsc.org/release/manualpages/PC/PCSetType/                       |
  |  https://petsc.org/release/manualpages/KSP/KSPSetFromOptions/              |
  |  https://petsc.org/release/manualpages/KSP/KSPSetUp/                       |
  |  https://petsc.org/release/manualpages/KSP/KSPSetInitialGuessNonzero/      |v
  +---------------------------------------------------------------------------*/
  void c_petsc_set_solver_and_preconditioner_(KSP  * ksp,
                                              PC   * pc,
                                              Mat  * A,
                                              char * sol,
                                              char * prec) {

    /*---------------------------------------------------------------------+
    |  Set operators. Here the matrix that defines the linear system       |
    |  also serves as the preconditioning matrix. Since all the matrices   |
    |  will have the same nonzero pattern here, we indicate this so the    |
    |  linear solvers can take advantage of this.                          |
    +---------------------------------------------------------------------*/

    /* Set precondioning matrix to be A */
    err = KSPSetOperators(*ksp, *A, *A);

    /* Set solver */
    err = KSPSetType(*ksp, sol);

    /* Set preconditioner */
    err = KSPGetPC(*ksp, pc);
    err = PCSetType(*pc, prec);

    /* These two lines are needed to finish the setup */
    err = KSPSetFromOptions(*ksp);
    err = KSPSetUp         (*ksp);

    /*----------------------------------------------------------+
    |  And please don't start from zero - for the sake of God   |
    +----------------------------------------------------------*/
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
