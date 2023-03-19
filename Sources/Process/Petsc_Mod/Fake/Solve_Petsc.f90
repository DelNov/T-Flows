!==============================================================================!
  subroutine Solve_Petsc(Pet,                      &
                         solver, prec, prec_opts,  &
                         A, x, b,                  &
                         miter, niter,             &
                         tol, fin_res,             &
                         norm)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)    :: Pet
  character(*)         :: solver          ! solver
  character(*)         :: prec            ! preconditioner
  character(SL)        :: prec_opts(MSI)  ! preconditioner options
  type(Matrix_Type)    :: A
  real                 :: x(-Pet % pnt_grid % n_bnd_cells :  &
                             Pet % pnt_grid % n_cells)
  real                 :: b( Pet % pnt_grid % n_cells)
  integer, intent(in)  :: miter
  integer, intent(out) :: niter
  real,    intent(in)  :: tol      ! tolerance
  real,    intent(out) :: fin_res  ! final residual
  real,    optional    :: norm     ! normalization
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

