!==============================================================================!
  subroutine Acm(sol, x, b, prec, niter, tol, ini_res, fin_res, norm)
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
  real              :: b( sol % pnt_grid % n_cells)      ! [A]{x}={b}
  character(len=80) :: prec                              ! preconditioner
  integer           :: niter                             ! number of iterations
  real              :: tol                               ! tolerance
  real              :: ini_res                           ! (unused)
  real              :: fin_res                           ! final residual
  real, optional    :: norm                              ! normalization
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Matrix_Type), pointer :: a
  type(Matrix_Type), pointer :: a_lev(:)
  type(Matrix_Type), pointer :: d_lev(:)
  type(Vector_Type), pointer :: x_lev(:)
  type(Vector_Type), pointer :: b_lev(:)
  type(Vector_Type), pointer :: r_lev(:)
  integer                    :: lev, c, c_lev
  real                       :: summ
!==============================================================================!

  grid  => sol % pnt_grid
  a     => sol % a
  a_lev => sol % a_lev
  d_lev => sol % d_lev
  x_lev => sol % x_lev
  b_lev => sol % b_lev
  r_lev => sol % r_lev

  !-------------------------------------------!
  !   Coarsen system matrix over all levels   !
  !-------------------------------------------!
  call Acm_Coarsen_Matrix(sol)

  !---------------------------------!
  !   Simply copy the first level   !
  !---------------------------------!
  x_lev(1) % val(1:) = x(1:)
  b_lev(1) % val(1:) = b(1:)

  lev = 1
  call Cg_Level(lev,               &  ! level
                a_lev(lev),        &  ! system matrix
                d_lev(lev),        &  ! preconditioning matrix
                x_lev(lev) % val,  &  ! solution
                b_lev(lev) % val,  &  ! right hand side
                prec,              &  ! preconditioner
                niter,             &  ! niter (for now)
                tol,               &  ! tolerance
                0.1,               &  ! residual ratio
                fin_res)              ! final residual

  !---------------------------------------------------------------!
  !   Compute residual {b} = {b} - [A]{x} and store it in b_lev   !
  !    (b_lev which was holding r.h.s. will be over-written))     !
  !---------------------------------------------------------------!
  call Residual_Vector(grid % level(lev) % n_cells,  &
                       r_lev(lev) % val,             &
                       b_lev(lev) % val,             &
                       a_lev(lev),                   &
                       x_lev(lev) % val)

  ! Check RHS
  summ = 0.0
  do c = 1, grid % level(lev) % n_cells
    summ = summ + r_lev(lev) % val(c)
  end do
  PRINT '(a,i2,es15.3)', 'SUMM @ ', lev, summ

  lev = 2
  b_lev(lev) % val(:) = 0.0
  do c = 1, grid % level(1) % n_cells
    c_lev = grid % level(lev) % cell(c)
    b_lev(lev) % val(c_lev) = b_lev(lev) % val(c_lev) + r_lev(1) % val(c)
  end do

  call Cg_Level(lev,               &  ! level
                a_lev(lev),        &  ! system matrix
                d_lev(lev),        &  ! preconditioning matrix
                x_lev(lev) % val,  &  ! solution
                b_lev(lev) % val,  &  ! right hand side
                prec,              &  ! preconditioner
                niter,             &  ! niter (for now)
                tol,               &  ! tolerance
                0.1,               &  ! residual ratio
                fin_res)              ! final residual

  lev = 1
  call Cg_Level(lev,               &  ! level
                a_lev(lev),        &  ! system matrix
                d_lev(lev),        &  ! preconditioning matrix
                x_lev(lev) % val,  &  ! solution
                b_lev(lev) % val,  &  ! right hand side
                prec,              &  ! preconditioner
                niter,             &  ! niter (for now)
                tol,               &  ! tolerance
                0.1,               &  ! residual ratio
                fin_res)              ! final residual

STOP

  end subroutine
