!==============================================================================!
  subroutine Solve(Amg, madapt, levels)
!------------------------------------------------------------------------------!
!   Solution phase of Amg1r5
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Arguments]----------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: madapt
  integer                 :: levels
!-----------------------------------[Locals]-----------------------------------!
  real    :: ama, cfac, epsi, epsil, fac, fmax
  real    :: rescg, resold, umax
  integer :: digit(AMG_MAX_LEVELS)
  integer :: i, iter, m, mfirst, msel, n
  integer :: n_digits
  real,    contiguous, pointer :: a(:), u(:), f(:)
  integer, contiguous, pointer :: ia(:)
!------------------------------------[Save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ! Number of unknowns at the finest level
  n  =  Amg % lev(1) % n
  a  => Amg % lev(1) % a
  ia => Amg % lev(1) % ia
  u  => Amg % lev(1) % u
  f  => Amg % lev(1) % f

  m = levels
  Amg % ncyc0 = 0
  if(Amg % eps .ne. 0.0) then
    epsi = Amg % eps
  else
    epsi = 1.0e-12
  end if

  !----------------------!
  !   Decompose madapt   !
  !----------------------!
  if (madapt.ne.0) then
    call Amg % Get_Integer_Digits(madapt, 2, n_digits, digit)
    msel = digit(1)
    if (msel.eq.2) then
      if (digit(2).ne.0) then
        fac = real(digit(2))
        do i = 1, 100
          fac = fac/10.0
          if(fac .le. 1.0) exit
        end do
      else
        fac = 0.7
      end if
    end if
  else
    msel = 2
    fac = 0.7
  end if

  !-------------------------------------------------!
  !   Default values proposed by Ruge and Stueben   !
  !-------------------------------------------------!
  ! Amg % cycle % type           = AMG_V_CYCLE
  ! Amg % cycle % cg_usage       = AMG_NO_CG_STEPS
  ! Amg % cycle % stop_criterion = AMG_PERFORM_ALL_CYCLES
  ! Amg % cycle % max_cycles     = 10

  !----------------------------------------------------------------!
  !   Set epsi according to convergence criterion given by iconv   !
  !----------------------------------------------------------------!
  if(Amg % cycle % stop_criterion .ne. AMG_STOP_IF_RES_LT_EPS_F) then
    if(Amg % cycle % stop_criterion .eq. AMG_STOP_IF_RES_LT_EPS_UD) then
    ama = 0.0
    do i = 1, n
      ama = max(ama, a(ia(i)))
    end do
    epsi = epsi*ama
    end if
  else
    fmax = 0.0
    do i = 1, n
      fmax = max(fmax, abs(f(i)))
    end do
    epsi = epsi * fmax
  end if

  !--------------------------------------------------------------!
  !   In the case coarse grid solution is not obtained with GS   !
  !   sweeps do not use coarsest grid with less than 10 points   !
  !--------------------------------------------------------------!
  if(Amg % coarse_solver .ne. AMG_SOLVER_GS) then
    do i = m, 1, -1
      if(Amg % imax(i) - Amg % imin(i) .ge. 9) exit
    end do
    m = i
    levels = i
  end if

  !-------------!
  !   Cycling   !
  !-------------!
  if(Amg % iout .ne. 0) then
    call Amg % Calculate_Residual(1, Amg % res0)
    if(Amg % iout .ge. 3) then
      write (6, 9005)
      write (6, 9000) Amg % res0
    end if
    resold = Amg % res0
  end if

  do iter = 1, Amg % cycle % max_cycles
    call Amg % Backup_U(1)
    call Amg % One_Cycle(1, m, iter, msel, fac, levels)
    call Amg % Cg_Step(1, iter)
    if(Amg % ierr.gt.0) return
    if(iter .eq. 1) then
      mfirst = m
      if(Amg % iout .ge. 3) write(6, 9040) m
    else if(Amg % iout .ge. 3.and.m.ne.mfirst) then
      mfirst = m
      write(6, 9040) m
    end if
    if(Amg % iout .ge. 3 .or.  &
       Amg % cycle % stop_criterion .ne. AMG_PERFORM_ALL_CYCLES) then
      call Amg % Calculate_Residual(1, Amg % res)
    end if
    Amg % ncyc0 = iter
    if(Amg % iout .ge. 3) then
      cfac = Amg % res / (resold + TINY)
      resold = Amg % res
      if(m .ne. 1) then
        call Amg % Calculate_Residual(2, rescg)
        write(6, 9010) iter, rescg, Amg % res, cfac
      else
        write(6, 9020) iter, Amg % res, cfac
      end if
    end if
    if(Amg % cycle % stop_criterion .ne. AMG_PERFORM_ALL_CYCLES) then
      epsil = epsi
      if(Amg % cycle % stop_criterion .eq. AMG_STOP_IF_RES_LT_EPS_UD) then
        umax = 0.0
        do i = 1, n
          umax = max(umax, abs(u(i)))
        end do
        epsil = epsi * umax
      end if
    end if
    if(Amg % res .lt. epsil) exit
  end do
  if(Amg % iout .lt. 3 .and. Amg % iout .ne. 0) then
    call Amg % Calculate_Residual(1, Amg % res)
  end if

9000  format (' cycle  0:',3x,'res=',d9.3)
9005  format (/' ************* cycling..... *************'/)
9010  format ( ' cycle ',i2,':   rescg=',d9.3,'   res=',d9.3,  &
               '   cfac=',d9.3)
9020  format ( ' cycle ',i2,':   res=',d9.3,'   cfac=',d9.3)
9040  format (/' cycling between grids 1 and',i3,':'/)

      end subroutine
