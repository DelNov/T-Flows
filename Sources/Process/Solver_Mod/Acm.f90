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
  !                                           !
  !   Coarsen system matrix over all levels   !
  !                                           !
  !-------------------------------------------!
  call Acm_Coarsen_Matrix(sol)

  !---------------------------------!
  !   Simply copy the first level   !
  !---------------------------------!
  x_lev(1) % val(1:) = x(1:)
  b_lev(1) % val(1:) = b(1:)

  lev = 1
  !---------------------------------!
  !   Solve    [A]{x} = b and       !
  !   compute  {r} = {b} - [A]{x}   !
  !---------------------------------!
  call Cg_Level(lev,               &  ! level
                a_lev(lev),        &  ! system matrix
                d_lev(lev),        &  ! preconditioning matrix
                x_lev(lev) % val,  &  ! solution
                b_lev(lev) % val,  &  ! right hand side
                r_lev(lev) % val,  &  ! residual vector
                prec,              &  ! preconditioner
                20,                &  ! max iterations
                tol,               &  ! tolerance
                0.1,               &  ! residual ratio
                fin_res)              ! final residual

  ! Check RHS
  summ = 0.0
  do c = 1, grid % level(lev) % n_cells
    summ = summ + r_lev(lev) % val(c)
  end do
  PRINT '(a,i2,es15.3)', 'SUMM @ ', lev, summ

  !-------------------------!
  !                         !
  !   Go down the V cycle   !
  !                         !
  !-------------------------!
  do lev = 2, 5

    !-----------------+
    !   Restriction   !                    !
    !-----------------+--------------------!
    !     +---+---+          +-------+     !
    !     |   |   |          |       |     !
    !     +---+---+   =-->   |       |     !
    !     |   |   |          |       |     !
    !     +---+---+          +-------+     !
    !       lev-1               lev        !
    !--------------------------------------!
    b_lev(lev) % val(:) = 0.0
    do c = 1, grid % level(lev-1) % n_cells       ! through finer cells
      c_lev = grid % level(lev-1) % coarser_c(c)  ! get coarse cell
      b_lev(lev) % val(c_lev) = b_lev(lev) % val(c_lev) + r_lev(lev-1) % val(c)
    end do

    !---------------------------------!
    !   Solve    [A]{x} = b and       !
    !   compute  {r} = {b} - [A]{x}   !
    !---------------------------------!
    call Cg_Level(lev,               &  ! level
                  a_lev(lev),        &  ! system matrix
                  d_lev(lev),        &  ! preconditioning matrix
                  x_lev(lev) % val,  &  ! solution
                  b_lev(lev) % val,  &  ! right hand side
                  r_lev(lev) % val,  &  ! residual vector
                  prec,              &  ! preconditioner
                  20,                &  ! max iterations
                  tol,               &  ! tolerance
                  0.1,               &  ! residual ratio
                  fin_res)              ! final residual

  end do

  lev = 1
  !---------------------------------!
  !   Solve    [A]{x} = b and       !
  !   compute  {r} = {b} - [A]{x}   !
  !---------------------------------!
  call Cg_Level(lev,               &  ! level
                a_lev(lev),        &  ! system matrix
                d_lev(lev),        &  ! preconditioning matrix
                x_lev(lev) % val,  &  ! solution
                b_lev(lev) % val,  &  ! right hand side
                r_lev(lev) % val,  &  ! residual vector
                prec,              &  ! preconditioner
                20,                &  ! niter (for now)
                tol,               &  ! tolerance
                0.1,               &  ! residual ratio
                fin_res)              ! final residual

STOP

  end subroutine
