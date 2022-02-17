#include <stdio.h>
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

PetscInt       one  = 1;
const PetscInt zero = 0;
PetscErrorCode err;

static char help[] = "This is to initialize PETSc from T-Flows.\n";

/*-----------------------------------------------------------------------------+
|  PetscInitialize                                                             |
+-----------------------------------------------------------------------------*/
void c_petsc_initialize_() {

  /* Issue PETSc call */
  err = PetscInitialize(0, NULL, (char*)0, help);
}

/*-----------------------------------------------------------------------------+
|                                                                              |
|  Matrix routines                                                             |
|                                                                              |
+-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------+
|  MatCreate                                                                   |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/Mat/MatCreate.html               |
+-----------------------------------------------------------------------------*/
void c_petsc_mat_create_(Mat * A) {

  err = MatCreate(MPI_COMM_WORLD, A);
}

/*-----------------------------------------------------------------------------+
|  MatSetSizes                                                                 |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/Mat/MatSetSizes.html             |
+-----------------------------------------------------------------------------*/
void c_petsc_mat_set_sizes_(Mat * A, PetscInt * m_lower, PetscInt * m_upper) {

  err = MatSetSizes(*A, *m_lower, *m_lower, *m_upper, *m_upper);
}

/*-----------------------------------------------------------------------------+
|  MatSetType (to MATAIJ)                                                      |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/Mat/MatSetType.html              |
+-----------------------------------------------------------------------------*/
void c_petsc_mat_set_type_to_mat_aij_(Mat * A) {

  err = MatSetType(*A, MATAIJ);
}

/*-----------------------------------------------------------------------------+
|  MatMPIAIJSetPreallocation                                                   |
|  MatSeqAIJSetPreallocation                                                   |
|                                                                              |
|  Combines two calls which seems to be important or necessary                 |
|  https://petsc.org/release/docs/manualpages/Mat/MATAIJ.html#MATAIJ           |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/Mat/MatMPIAIJSetPreallocation.html
|  https://petsc.org/release/docs/manualpages/Mat/MatSeqAIJSetPreallocation.html
+-----------------------------------------------------------------------------*/
void c_petsc_mat_aij_set_preallocation_(Mat      * A,
                                        PetscInt * m_lower,
                                        PetscInt * d_nnz,
                                        PetscInt * m_upper,
                                        PetscInt * o_nnz) {

  err = MatMPIAIJSetPreallocation(*A, *m_lower, d_nnz, *m_upper, o_nnz);
  err = MatSeqAIJSetPreallocation(*A, *m_lower, d_nnz);
}

/*-----------------------------------------------------------------------------+
|  MatSetValue                                                                 |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/Mat/MatSetValue.html             |
+-----------------------------------------------------------------------------*/
void c_petsc_mat_set_value_(Mat      * A,
                            PetscInt * row,
                            PetscInt * col,
                            double   * value) {

  err = MatSetValue(*A, *row, *col, *value, INSERT_VALUES);
}

/*-----------------------------------------------------------------------------+
|  MatAssemblyBegin and MatAssemblyEnd                                         |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/Mat/MatAssemblyBegin.html        |
|  https://petsc.org/release/docs/manualpages/Mat/MatAssemblyEnd.html          |
+-----------------------------------------------------------------------------*/
void c_petsc_assemble_mat_(Mat * A) {

  err = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
  err = MatAssemblyEnd  (*A, MAT_FINAL_ASSEMBLY);
}

/*-----------------------------------------------------------------------------+
|                                                                              |
|  Vector routines                                                             |
|                                                                              |
+-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------+
|  VecCreateMPI                                                                |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/Vec/VecCreateMPI.html            |
+-----------------------------------------------------------------------------*/
void c_petsc_vec_create_mpi_(Vec * v, PetscInt * m_lower, PetscInt * m_upper) {

  err = VecCreateMPI(MPI_COMM_WORLD, *m_lower, *m_upper, v);
}


/*-----------------------------------------------------------------------------+
|  VecSetValue                                                                 |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/Vec/VecSetValue.html             |
+-----------------------------------------------------------------------------*/
void c_petsc_vec_set_value_(Vec      * v,
                            PetscInt * row,
                            double   * value) {

  err = VecSetValue(*v, *row, *value, INSERT_VALUES);
}

/*-----------------------------------------------------------------------------+
|  VecAssemblyBegin and VecAssemblyEnd                                         |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/Vec/VecAssemblyBegin.html        |
|  https://petsc.org/release/docs/manualpages/Vec/VecAssemblyEnd.html          |
+-----------------------------------------------------------------------------*/
void c_petsc_assemble_vec_(Vec * v) {

  err = VecAssemblyBegin(*v);
  err = VecAssemblyEnd(*v);
}

/*-----------------------------------------------------------------------------+
|  VecGetValues                                                                |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/Vec/VecGetValues.html            |
+-----------------------------------------------------------------------------*/
void c_petsc_vec_get_values_(Vec      * v,
                             PetscInt * ni,
                             PetscInt * row,
                             double   * value) {

  err = VecGetValues(*v, *ni, row, value);
}

/*-----------------------------------------------------------------------------+
|                                                                              |
|  Solver routines                                                             |
|                                                                              |
+-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------+
|  KSPCreate                                                                   |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/KSP/KSPCreate.html               |
+-----------------------------------------------------------------------------*/
void c_petsc_ksp_create_(KSP * ksp) {

  err = KSPCreate(MPI_COMM_WORLD, ksp);
}

/*-----------------------------------------------------------------------------+
|  KSP routines to set solver and preconditioner                               |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/KSP/KSPSetOperators.html         |
|  https://petsc.org/release/docs/manualpages/KSP/KSPSetType.html              |
|  https://petsc.org/release/docs/manualpages/KSP/KSPGetPC.html                |
|  https://petsc.org/release/docs/manualpages/PC/PCSetType.html                |
|  https://petsc.org/release/docs/manualpages/KSP/KSPSetFromOptions.html       |
|  https://petsc.org/release/docs/manualpages/KSP/KSPSetUp.html                |
|  https://petsc.org/release/docs/manualpages/KSP/KSPSetInitialGuessNonzero.html
+-----------------------------------------------------------------------------*/
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

/*-----------------------------------------------------------------------------+
|  KSPSetTolerances                                                            |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/KSP/KSPSetTolerances.html        |
+-----------------------------------------------------------------------------*/
void c_petsc_ksp_set_tolerances_(KSP      * ksp,
                                 double   * rtol,
                                 double   * atol,
                                 PetscInt * maxits) {

  err = KSPSetTolerances(*ksp, *rtol, *atol, 1.0e+3, *maxits);
}

/*-----------------------------------------------------------------------------+
|  KSPSetTolerances                                                            |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/KSP/KSPSolve.html                |
+-----------------------------------------------------------------------------*/
void c_petsc_ksp_solve_(KSP * ksp,
                        Vec * b,
                        Vec * x) {

  /* Issue PETSc call */
  err = KSPSolve(*ksp, *b, *x);
}

/*-----------------------------------------------------------------------------+
|  KSPGetIterationNumber                                                       |
|                                                                              |
|  https://petsc.org/release/docs/manualpages/KSP/KSPGetIterationNumber.html   |
+-----------------------------------------------------------------------------*/
void c_petsc_ksp_get_iteration_number_(KSP      * ksp,
                                       PetscInt * its) {

  /* Issue PETSc call */
  err = KSPGetIterationNumber(*ksp, its);
}
