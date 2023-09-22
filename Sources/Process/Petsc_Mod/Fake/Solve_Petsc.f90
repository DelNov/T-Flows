!==============================================================================!
  subroutine Solve_Petsc(Pet,                      &
                         solver, prec, prec_opts,  &
                         A, x, b,                  &
                         miter, niter,             &
                         tol, fin_res,             &
                         blend_matrix)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)          :: Pet
  character(*),  intent(in)  :: solver                       ! solver
  character(*),  intent(in)  :: prec                         ! preconditioner
  character(SL), intent(in)  :: prec_opts(MAX_STRING_ITEMS)  ! prec. options
  type(Matrix_Type)          :: A
  real                       :: x(-Pet % pnt_grid % n_bnd_cells :  &
                                   Pet % pnt_grid % n_cells)
  real                       :: b( Pet % pnt_grid % n_cells)
  integer,       intent(in)  :: miter
  integer,       intent(out) :: niter
  real,          intent(in)  :: tol      ! tolerance
  real,          intent(out) :: fin_res  ! final residual
  logical,       intent(in)  :: blend_matrix
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

  ! Just to avoid compiler's warning
  niter   = 0
  fin_res = 0.0

  call Message % Error(60,                                                   &
                  'This version was compiled without PETSc, and yet '    //  &
                  'they are specified in the control file.  This error ' //  &
                  'is critical, exiting. Eithrt fix the control file '   //  &
                  'by setting LINEAR_SOLVERS to native, or install   '   //  &
                  'PETSc on your system',                                    &
                  file=__FILE__, line=__LINE__, one_proc=.true.)

  end subroutine

