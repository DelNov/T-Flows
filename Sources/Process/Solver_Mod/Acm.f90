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
  integer                    :: s, c1, c2, lev, s_lev, c1_lev, c2_lev
  integer                    :: dia_11, dia_22, pos_12, pos_21
  real                       :: a_12, a_21
  integer                    :: c_lev, dia_pos, i, j
!==============================================================================!

  grid  => sol % pnt_grid
  a     => sol % a
  a_lev => sol % a_lev

  !---------------------------------!
  !                                 !
  !   Simply copy the first level   !
  !                                 !
  !---------------------------------!
  a_lev(1) % val(:) = a % val(:)

  !--------------------------!
  !                          !
  !   Coarsen other levels   !
  !                          !
  !--------------------------!
  do lev = 2, grid % n_levels

    do s = 1, grid % level(1) % n_faces    ! go through finest faces

      s_lev = grid % level(lev) % face(s)  ! check if face expands to lev

      if(s_lev > 0) then

        ! Matrix entries in the finest level
        a_12 = a_lev(1) % val(a_lev(1) % pos(1, s))
        a_21 = a_lev(1) % val(a_lev(1) % pos(2, s))

        c1_lev = grid % level(lev) % faces_c(1, s_lev)   ! cell 1 at level lev
        c2_lev = grid % level(lev) % faces_c(2, s_lev)   ! cell 2 at level lev

        dia_11 = a_lev(lev) % dia(c1_lev)  ! cell position in the matrix
        dia_22 = a_lev(lev) % dia(c2_lev)  ! cell position in the matrix

        pos_12 = a_lev(lev) % pos(1, s_lev)  ! face position in the matrix
        pos_21 = a_lev(lev) % pos(2, s_lev)  ! face position in the matrix

        a_lev(lev) % val(dia_11) = a_lev(lev) % val(dia_11) - a_12
        a_lev(lev) % val(dia_22) = a_lev(lev) % val(dia_22) - a_21

        a_lev(lev) % val(pos_12) = a_lev(lev) % val(pos_12) + a_12
        a_lev(lev) % val(pos_21) = a_lev(lev) % val(pos_21) + a_21
      end if

    end do   ! face
  end do  ! level

  !----------------------------------!
  !                                  !
  !   Call smoother just for kicks   !
  !                                  !
  !----------------------------------!
  do lev = 1, grid % n_levels
    call Cg_Level(a_lev(lev),        &  ! system matrix
                  sol % d_lev(lev),        &  ! preconditioning matrix
                  sol % x_lev(lev) % val,  &  ! solution
                  sol % b_lev(lev) % val,  &  ! right hand side
                  lev,                     &  ! mg level
                  niter,                   &  ! niter (for now)
                  tol,                     &  ! tolerance
                  ini_res,                 &  ! initial residual (unused)
                  fin_res,                 &  ! final residual
                  norm)                       ! not sure if really needed
    print *, '- residual at level ', lev, ':', fin_res
  end do

STOP

  end subroutine
