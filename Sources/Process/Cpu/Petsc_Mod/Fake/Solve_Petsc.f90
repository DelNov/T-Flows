!==============================================================================!
  subroutine Solve_Petsc(Pet,                      &
                         solver, prec, prec_opts,  &
                         A, x, b,                  &
                         miter, niter,             &
                         tol, fin_res,             &
                         blend_matrix)
!------------------------------------------------------------------------------!
!>  This subroutine serves as placeholders or dummy variant when PETSc is
!>  not available in the system.  The functionality here is intentionally
!>  minimal, primarily serving to maintain code structure and compatibility
!>  in environments where PETSc is not used.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)          :: Pet      !! parent object of the Petsc_Type
  character(*),  intent(in)  :: solver   !! name of the solver to use
  character(*),  intent(in)  :: prec     !! name of the preconditioner to use
  character(SL), intent(in)  :: prec_opts(MAX_STRING_ITEMS)
    !! list of options passed to preconditioner from the control file
  type(Matrix_Type)          :: A        !! matrix in T-Flows format
  real                       :: x(-Pet % pnt_grid % n_bnd_cells :  &
                                   Pet % pnt_grid % n_cells)  !! unknown vector
  real                       :: b( Pet % pnt_grid % n_cells)  !! source vector
  integer,       intent(in)  :: miter    !! maximum number of iterations
  integer,       intent(out) :: niter    !! performed number of iterations
  real,          intent(in)  :: tol      !! target solver tolerance
  real,          intent(out) :: fin_res  !! final residual after linear solver
  logical,       intent(in)  :: blend_matrix
    !! flag indicating if the system matrix is blended with upwind terms
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Pet)
  Unused(solver)
  Unused(prec)
  Unused(prec_opts)
  Unused(A)
  Unused(x)
  Unused(b)
  Unused(miter)
  Unused(niter)
  Unused(tol)
  Unused(fin_res)
  Unused(blend_matrix)
!==============================================================================!

  end subroutine

