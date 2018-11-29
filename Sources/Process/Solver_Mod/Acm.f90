!==============================================================================!
  subroutine Acm(sol, x, r1, prec, niter, tol, ini_res, fin_res, norm)
!------------------------------------------------------------------------------!
!   This is the nucleus of the Additive Correction Multigrid method.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Matrix_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Solver_Type), target :: sol
  real              :: x(-sol % pnt_grid % n_bnd_cells :  &
                          sol % pnt_grid % n_cells)
  real              :: r1(sol % pnt_grid % n_cells)      ! [A]{x}={r1}
  character(len=80) :: prec                              ! preconditioner
  integer           :: niter                             ! number of iterations
  real              :: tol                               ! tolerance
  real              :: ini_res, fin_res                  ! residual
  real, optional    :: norm                              ! normalization
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Matrix_Type), pointer :: a
  type(Matrix_Type), pointer :: a_lev(:)
  type(Matrix_Type), pointer :: d_lev(:)
  type(Vector_Type), pointer :: x_lev(:)
  type(Vector_Type), pointer :: b_lev(:)
  integer                    :: lev
!==============================================================================!

  grid  => sol % pnt_grid
  a     => sol % a
  a_lev => sol % a_lev
  d_lev => sol % d_lev
  x_lev => sol % x_lev
  b_lev => sol % b_lev

  !-------------------------------------------!
  !   Coarsen system matrix over all levels   !
  !-------------------------------------------!
  call Acm_Coarsen_Matrix(sol)

  !---------------------------------!
  !   Simply copy the first level   !
  !---------------------------------!
  x_lev(1) % val(1:) = x(1:)
  b_lev(1) % val(1:) = r1(1:)

  lev = 1
  call Cg_Level(lev,               &  ! level
                a_lev(lev),        &  ! system matrix
                d_lev(lev),        &  ! preconditioning matrix
                x_lev(lev) % val,  &  ! solution
                b_lev(lev) % val,  &  ! right hand side
                niter,             &  ! niter (for now)
                tol,               &  ! tolerance
                ini_res,           &  ! initial residual (unused)
                fin_res)              ! final residual
STOP

  end subroutine
