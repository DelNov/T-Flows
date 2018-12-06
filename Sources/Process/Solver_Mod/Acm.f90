!==============================================================================!
  subroutine Acm(sol, x, b, prec, n_cycles, tol, ini_res, fin_res, norm)
!------------------------------------------------------------------------------!
!   This is the nucleus of the Additive Correction Multigrid method.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Matrix_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Solver_Type), target :: sol
  real              :: x(-sol % pnt_grid % n_bnd_cells :  &
                          sol % pnt_grid % n_cells)
  real              :: b( sol % pnt_grid % n_cells)      ! [A]{x}={b}
  character(len=80) :: prec                              ! preconditioner
  integer           :: n_cycles                          ! number of cycles
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
  type(Vector_Type), pointer :: p_lev(:)
  type(Vector_Type), pointer :: r_lev(:)
  integer                    :: lev, c, c_lev, lev_max, cyc
  integer                    :: n_iter, p_iter  ! requested and performed iters
  real                       :: res_rat         ! requested residual ratio
!==============================================================================!

  ! Take aliases first
  grid  => sol % pnt_grid
  a     => sol % a
  a_lev => sol % a_lev
  d_lev => sol % d_lev
  x_lev => sol % x_lev
  b_lev => sol % b_lev
  p_lev => sol % p_lev
  r_lev => sol % r_lev

  ! Set maximum number of levels to perform
  call Control_Mod_V_Cycle_Max_Grid_Levels(lev_max)
  lev_max = min(lev_max, grid % n_levels)

  ! Set number of smoothing iterations
  n_iter = 10  ! 10 seems to work the best for IC preconditioner
  call Control_Mod_V_Cycle_Number_Of_Smoothing_Iterations(n_iter)

  ! Set requested residual ration in the level
  res_rat = 0.0001
  call Control_Mod_V_Cycle_Residual_Ratio(res_rat)

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
  r_lev(1) % val(1:) = 0.0

  !------------------------------------!
  !                                    !
  !                                    !
  !   Solve on the finest level once   !
  !                                    !
  !    Solve    [A]{x} = b and         !
  !    Compute  {r} = {b} - [A]{x}     !
  !                                    !
  !                                    !
  !------------------------------------!
  call Cg_Level(1,               &  ! level
                a_lev(1),        &  ! system matrix
                d_lev(1),        &  ! preconditioning matrix
                x_lev(1) % val,  &  ! solution
                b_lev(1) % val,  &  ! right hand side
                p_lev(1) % val,  &  ! p1 vector from cg algorith
                r_lev(1) % val,  &  ! residual vector
                prec,            &  ! preconditioner
                1,               &  ! cycle
                n_iter,          &  ! max iterations
                p_iter,          &  ! performed iterations
                tol,             &  ! tolerance
                res_rat,         &  ! residual ratio
                fin_res)            ! final residual

  !--------------------!
  !                    !
  !                    !
  !   Start V cycles   !
  !                    !
  !                    !
  !--------------------!
  do cyc = 1, n_cycles

    !-------------------------!
    !                         !
    !   Go down the V cycle   !
    !                         !
    !-------------------------!
    do lev = 2, lev_max

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
      if(p_iter > 0) then  ! restriction only if it did something in the solver
        b_lev(lev) % val(:) = 0.0
        do c = 1, grid % level(lev-1) % n_cells       ! through finer cells
          c_lev = grid % level(lev-1) % coarser_c(c)  ! get coarse cell
          b_lev(lev) % val(c_lev) =  &
          b_lev(lev) % val(c_lev) + r_lev(lev-1) % val(c)
        end do
      end if

      !---------------------------------!
      !   Solve    [A]{x} = b and       !
      !   compute  {r} = {b} - [A]{x}   !
      !---------------------------------!
      call Cg_Level(lev,               &  ! level
                    a_lev(lev),        &  ! system matrix
                    d_lev(lev),        &  ! preconditioning matrix
                    x_lev(lev) % val,  &  ! solution
                    b_lev(lev) % val,  &  ! right hand side
                    p_lev(lev) % val,  &  ! p1 vector from cg algorith
                    r_lev(lev) % val,  &  ! residual vector
                    prec,              &  ! preconditioner
                    1,                 &
                    n_iter,            &  ! max iterations
                    p_iter,            &  ! performed iterations
                    tol,               &  ! tolerance
                    res_rat,           &  ! residual ratio
                    fin_res)              ! final residual

    end do  ! end of going down

                    !----------------------------------!
                    !                                  !
                    !    You reched the rock bottom    !
                    !   solved at the coarsest level   !
                    !                                  !
                    !----------------------------------!

    !-----------------------!
    !                       !
    !   Go up the V cycle   !
    !                       !
    !-----------------------!
    do lev = lev_max-1, 1, -1

      !------------------+
      !   Prolongation   !
      !------------------+-------------------!
      !     +-------+          +---+---+     !
      !     |       |          |   |   |     !
      !     |       |   =-->   +---+---+     !
      !     |       |          |   |   |     !
      !     +-------+          +---+---+     !
      !       lev+1               lev        !
      !--------------------------------------!
      if(p_iter > 0) then  ! prolongation only if it did something in the solver
        do c = 1, grid % level(lev) % n_cells       ! through finer cells
          c_lev = grid % level(lev) % coarser_c(c)  ! get coarse cell
          x_lev(lev) % val(c) =  &
          x_lev(lev) % val(c) + x_lev(lev+1) % val(c_lev)
        end do
      end if

      !---------------------------------!
      !   Solve    [A]{x} = b and       !
      !   compute  {r} = {b} - [A]{x}   !
      !---------------------------------!
      call Cg_Level(lev,               &  ! level
                    a_lev(lev),        &  ! system matrix
                    d_lev(lev),        &  ! preconditioning matrix
                    x_lev(lev) % val,  &  ! solution
                    b_lev(lev) % val,  &  ! right hand side
                    p_lev(lev) % val,  &  ! p1 vector from cg algorith
                    r_lev(lev) % val,  &  ! residual vector
                    prec,              &  ! preconditioner
                    1,                 &  ! cycle
                    n_iter,            &  ! max iterations
                    p_iter,            &  ! performed iterations
                    tol,               &  ! tolerance
                    res_rat,           &  ! residual ratio
                    fin_res)              ! final residual
    end do
  end do

  end subroutine
