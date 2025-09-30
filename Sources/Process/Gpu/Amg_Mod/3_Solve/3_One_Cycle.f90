!==============================================================================!
  subroutine One_Cycle(Amg, finest_level, coarsest_level,  &
                       iter, msel, fac, levels)
!------------------------------------------------------------------------------!
!   Performs one Amg cycle with grid finest_level as finest grid
!
!   During first cycle: initialize parameters controlling yale-smp
!   factorization and adaptive determination of coarsest grid:
!   ifac=1: on next call of yale-smp factorize matrix
!   ifi =1: on first return to next to coarsest grid after coarse
!           grid solution compare residual with the residual on the
!           same grid before coarse grid solution. if reduction of
!           residual not satisfying, reduce number of grids used in
!           cycling by one and repete the process with the now
!           coarsest grid.
!   nptsf:  number of grid points on finest grid used in cycle,
!           divided by 10. only grids with less then nptsf points
!           are allowed to become coarsest grid during adaptive
!           coarse grid determination.
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Arguments]----------------------------------!
  class(Amg_Type) :: Amg
  integer         :: finest_level, coarsest_level, iter, msel
  real            :: fac
  integer         :: levels
!-----------------------------------[Locals]-----------------------------------!
  real    :: res
  integer :: ng(AMG_MAX_LEVELS)
  integer :: ifac, ifi, ivstar, level, mink, n, nl, nptsf
  logical :: moredown
!------------------------------------[Save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(iter .eq. 1) ifac = 1
  if(msel .ne. 2) then
    ifi   = 0
    nptsf = 0
  else
    mink  = AMG_BIG_INTEGER
    ifi   = 1
    nptsf = (Amg % imax(finest_level) - Amg % imin(finest_level)+1) / 10
  end if

  !-------------------!
  !                   !
  !   One grid only   !
  !                   !
  !-------------------!
  if(finest_level .ge. coarsest_level) then
    call Amg % Solve_On_Coarsest_Level(coarsest_level, ifac)
    call Amg % Normalize_U(finest_level)
    return
  end if

  !-------------------------------------------------------------!
  !                                                             !
  !   More then one grid (level seems to be the grid counter)   !
  !                                                             !
  !-------------------------------------------------------------!

  ! Initialize "ng", what is it, counter for levels of a sort?
  do level = finest_level, coarsest_level
    ng(level) = 0
  end do

  ! OK, let's see the following line will transform:
  ! AMG_V_CYCLE      -> AMG_V_STAR_CYCLE
  ! AMG_V_STAR_CYCLE -> AMG_V_CYCLE
  ! AMG_F_CYCLE      ->  0 (not defined)
  ! AMG_W_CYCLE      -> -1 =:-o
  ivstar = 3 - Amg % cycle % type

  ! Start from grid "finest_level"
  level = finest_level

  !------------------------------------------------------!
  !                                                      !
  !   Relax (downwards) (going towards coarser grids?)   !
  !                                                      !
  !------------------------------------------------------!
  downward: do
    moredown = .false.
      ! This might not cut it for fine(r) grids
      if(Amg % fine_solver .eq. AMG_SOLVER_GS) then
        do n = 1, Amg % n_relax_down
          do nl = 1, Amg % nrdlen
            call Amg % Gauss_Seidel_Sweep(level, Amg % type_relax_down(nl))
          end do
        end do
      else if(Amg % fine_solver .eq. AMG_SOLVER_CG) then
        call Amg % Cg_On_Level(level, Amg % n_relax_down)

      else if(Amg % fine_solver .eq. AMG_SOLVER_BICG) then
        call Amg % Bicg_On_Level(level, Amg % n_relax_down)

      end if
      if(ifi .eq. 1 .and. Amg % imax(level) - Amg % imin(level) .lt. nptsf) then

       !-------------------------------------------------!
       !   mink: lowest grid number for which residual   !
       !   is stored during first downwards relaxation   !
       !-------------------------------------------------!
        if(mink .eq. AMG_BIG_INTEGER) mink = level
        call Amg % Calculate_Residual(level, Amg % resi(level))
      end if

      ng(level) = ng(level)+1

      ! Increase the grid counter - going for coarser still
      level = level + 1
      call Amg % Set_U_To_Zero(level)
      call Amg % Restrict_Residuals(level, level-1)
      if(level .lt. coarsest_level) cycle downward

      !----------------------------!
      !                            !
      !   Solve on coarsest grid   !
      !                            !
      !----------------------------!
      call Amg % Solve_On_Coarsest_Level(coarsest_level, ifac)

      !---------------------!
      !                     !
      !   Relax (upwards)   !
      !                     !
      !---------------------!
      upward: do
        call Amg % Scale_Solution(level, ivstar)
        level = level-1
        call Amg % Interpolate_Correction(level)

        ! This might not cut it for fine(r) grids
        if(Amg % fine_solver .eq. AMG_SOLVER_GS) then
          do n = 1, Amg % n_relax_up
            do nl = 1, Amg % nrulen
              call Amg % Gauss_Seidel_Sweep(level, Amg % type_relax_up(nl))
            end do
          end do
        else if(Amg % fine_solver .eq. AMG_SOLVER_CG) then
          call Amg % Cg_On_Level(level, Amg % n_relax_up)

        else if(Amg % fine_solver .eq. AMG_SOLVER_BICG) then
          call Amg % Bicg_On_Level(level, Amg % n_relax_up)

        end if
        if(ifi .eq. 1 .and. level .ge. mink) then

          !----------------------------------------------!
          !   On first return to next to coarsest grid   !
          !    compare residual with the previous one    !
          !----------------------------------------------!
          call Amg % Calculate_Residual(level, res)

          ! If residual reduction satisfying:
          !   coarse grid adaption finished
          if(res .lt. Amg % resi(level)*fac) then
            ifi  = 0

          ! If residual reduction not satisfying:
          ! continue with coarse grid adaptation
          else
            ! Yale: if(Amg % coarse_solver .eq. 2) levels = level
            coarsest_level = level
            ifac = 1
            call Amg % Set_U_To_Zero(level)
            call Amg % Restrict_Residuals(level, level-1)
            call Amg % Solve_On_Coarsest_Level(coarsest_level, ifac)
            cycle upward
          end if
        end if

        if(level .eq. finest_level) then
          call Amg % Normalize_U(finest_level)
          return
        end if

        !------------------------------------------------!
        !   Grid switching corresponding to cycle type   !
        !------------------------------------------------!
        if(Amg % cycle % type .lt. AMG_F_CYCLE) then

          if(Amg % cycle % type .ge. 0) cycle upward

          if(level .eq. finest_level+1 .and.  &
             ng(level) .lt. abs(Amg % cycle % type)) then
            moredown = .true.
            exit upward
          end if
          cycle upward
        end if

        if(ng(level) .lt. 2) then
          moredown = .true.
          exit upward
        end if

        if(Amg % cycle % type .eq. AMG_W_CYCLE)  ng(level) = 0
      end do upward

    if(.not. moredown) exit downward
  end do downward

  end subroutine
