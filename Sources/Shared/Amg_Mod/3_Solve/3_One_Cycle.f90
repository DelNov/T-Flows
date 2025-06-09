!==============================================================================!
  subroutine One_Cycle(Amg, l, igam,     &
                       a, u, f, ia, ja,  &
                       iw, ifg, icg,     &
                       m, iter, msel, fac, levels)
!------------------------------------------------------------------------------!
!   Performs one Amg cycle with grid l as finest grid
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
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: l
  integer         :: igam
  real            :: a(:), u(:), f(:)
  integer         :: ia(:), ja(:), iw(:), ifg(:), icg(:)
  integer         :: m, iter, msel
  real            :: fac
  integer         :: levels
!-----------------------------------[locals]-----------------------------------!
  real    :: res
  integer :: ng(AMG_MAX_LEVELS)
  integer :: ifac, ifi, ivstar, level, mink, n, nl, nptsf
  integer :: nda
  logical :: moredown
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  nda = size(a, 1)

  if(iter .eq. 1) ifac = 1
  if(msel .ne. 2) then
    ifi   = 0
    nptsf = 0
  else
    mink  = AMG_BIG_INTEGER
    ifi   = 1
    nptsf = (Amg % imax(l) - Amg % imin(l)+1) / 10
  endif

  !-------------------!
  !                   !
  !   One grid only   !
  !                   !
  !-------------------!
  if(l .ge. m) then
    call Amg % Solve_On_Coarsest_Level(m, ifac,          &
                                       a, u, f, ia, ja,  &
                                       iw, icg)
    call Amg % Normalize_U(l, u)
    return
  end if

  !---------------------------------------------------------!
  !                                                         !
  !   More then one grid (level seems to be the grid counter)   !
  !                                                         !
  !---------------------------------------------------------!

  ! Initialize "ng", what is it, counter for levels of a sort?
  do level = l, m
    ng(level) = 0
  end do
  ivstar = 3 - igam

  ! Start from grid "l"
  level = l

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
            call Amg % Gauss_Seidel_Sweep(level, Amg % type_relax_down(nl),  &
                                          a, u, f, ia, ja,                   &
                                          iw, icg)
          end do
        end do
      else if(Amg % fine_solver .eq. AMG_SOLVER_CG) then
        do n = 1, Amg % n_relax_down
          call Amg % Cg_On_Level(level, 2,         &
                                 a, u, f, ia, ja,  &
                                 iw,icg)
        end do
      else if(Amg % fine_solver .eq. AMG_SOLVER_BICG) then
        do n = 1, Amg % n_relax_down
          call Amg % Bicg_On_Level(level, 2,         &
                                   a, u, f, ia, ja,  &
                                   iw,icg)
        end do
      end if
      if(ifi .eq. 1 .and. Amg % imax(level) - Amg % imin(level) .lt. nptsf) then

       !-------------------------------------------------!
       !   mink: lowest grid number for which residual   !
       !   is stored during first downwards relaxation   !
       !-------------------------------------------------!
        if(mink .eq. AMG_BIG_INTEGER) mink = level
        call Amg % timer_start()
        call Amg % Calculate_Residual(level, Amg % resi(level),  &
                                      a, u, f, ia, ja, iw)
        call Amg % timer_stop(15)
      end if

      ng(level) = ng(level)+1

      ! Increase the grid counter - going for coarser still
      level = level + 1
      call Amg % Set_U_To_Zero(level,u)
      call Amg % Restrict_Residuals(level, a, u, f, ia, ja, iw, ifg)
      if(level .lt. m) cycle downward

      !-------------------------------------------------------!
      !                                                       !
      !   Solve on coarsest grid - "m" is the coarsest grid   !
      !                                                       !
      !-------------------------------------------------------!
      call Amg % Solve_On_Coarsest_Level(m, ifac,          &
                                         a, u, f, ia, ja,  &
                                         iw, icg)

      !---------------------!
      !                     !
      !   Relax (upwards)   !
      !                     !
      !---------------------!
      upward: do
        call Amg % Scale_Solution(level, ivstar,  &
                                  a, u, f, ia, ja, iw)
        level = level-1
        call Amg % Interpolate_Correction(level, a, u, ia, ja, iw, ifg)

        ! This might not cut it for fine(r) grids
        if(Amg % fine_solver .eq. AMG_SOLVER_GS) then
          do n = 1, Amg % n_relax_up
            do nl = 1, Amg % nrulen
              call Amg % Gauss_Seidel_Sweep(level, Amg % type_relax_up(nl),  &
                                            a, u, f, ia, ja,                 &
                                            iw,icg)
            end do
          end do
        else if(Amg % fine_solver .eq. AMG_SOLVER_CG) then
          do n = 1, Amg % n_relax_up
            call Amg % Cg_On_Level(level, 2,         &
                                   a, u, f, ia, ja,  &
                                   iw,icg)
          end do
        else if(Amg % fine_solver .eq. AMG_SOLVER_BICG) then
          do n = 1, Amg % n_relax_up
            call Amg % Bicg_On_Level(level, 2,         &
                                     a, u, f, ia, ja,  &
                                     iw,icg)
          end do
        end if
        if(ifi .eq. 1 .and. level .ge. mink) then

          !----------------------------------------------!
          !   On first return to next to coarsest grid   !
          !    compare residual with the previous one    !
          !----------------------------------------------!
          call Amg % timer_start()
          call Amg % Calculate_Residual(level,res,  &
                                   a,u,f,ia,ja,iw)
          call Amg % timer_stop(15)

          ! If residual reduction satisfying:
          !   coarse grid adaption finished
          if(res .lt. Amg % resi(level)*fac) then
            ifi  = 0

          ! If residual reduction not satisfying:
          ! continue with coarse grid adaptation
          else
            ! Yale: if(Amg % coarse_solver .eq. 2) levels = level
            m    = level
            ifac = 1
            call Amg % Set_U_To_Zero(level,u)
            call Amg % Restrict_Residuals(level, a, u, f, ia, ja, iw, ifg)
            call Amg % Solve_On_Coarsest_Level(m, ifac,          &
                                               a, u, f, ia, ja,  &
                                               iw, icg)
            cycle upward
          endif
        end if

        if(level .eq. l) then
          call Amg % Normalize_U(l, u)
          return
        end if

        !------------------------------------------!
        !   Grid switching corresponding to igam   !
        !------------------------------------------!
        if(igam .lt. 3) then
          if(igam .ge. 0) cycle upward
          if(level .eq. l+1 .and. ng(level) .lt. abs(igam)) then
            moredown = .true.
            exit upward
          end if
          cycle upward
        end if
        if(ng(level) .lt. 2) then
          moredown = .true.
          exit upward
        end if

        if(igam .eq. 4)  ng(level) = 0
      end do upward

    if(.not. moredown) exit downward
  end do downward

  end subroutine
