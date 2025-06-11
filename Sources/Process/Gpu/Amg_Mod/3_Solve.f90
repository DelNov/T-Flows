!==============================================================================!
  subroutine Solve(Amg, madapt, ncyc,  &
                   a, u, f, ia, ja,    &  ! holds linear system
                   iw, icg, ifg,       &
                   levels)
!------------------------------------------------------------------------------!
!   Solution phase of Amg1r5
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: madapt, ncyc
  real            :: a(:), u(:), f(:)
  integer         :: ia(:), ja(:)
  integer         :: iw(:), icg(:), ifg(:)
  integer         :: levels
!-----------------------------------[locals]-----------------------------------!
  real    :: ama, cfac, epsi, epsil, fac, fmax
  real    :: rescg, resold, umax
  integer :: digit(AMG_MAX_LEVELS)
  integer :: i, icgr, iconv, igam, iter, l, m, mfirst, msel, n
  integer :: ncycle, n_digits
  integer :: ndu
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ndu = size(u, 1)

  !-------------------------------!
  !   Test of available storage   !
  !-------------------------------!
  if(ndu .lt. Amg % imax(levels)) then
    write(6, '(a)')  &
      ' *** error in Solve: ndu too small ***'
    Amg % ierr = AMG_ERR_DIM_U_TOO_SMALL
    return
  end if
  if(ndu .lt. Amg % imax(levels)) then
    write(6, '(a)')  &
      ' *** error in Solve: ndu too small ***'
    Amg % ierr = AMG_ERR_DIM_F_TOO_SMALL
    return
  end if

  m = levels
  Amg % ncyc0 = 0
  do n = 11, 20
    Amg % time(n) = 0.0
  end do
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

  !--------------------!
  !   Decompose ncyc   !
  !--------------------!
  if(ncyc .ne. 0) then
    call Amg % Get_Integer_Digits(abs(ncyc), 4, n_digits, digit)
    igam   = sign(digit(1), ncyc)
    icgr   = digit(2)
    iconv  = digit(3)
    ncycle = digit(4)
    if(ncycle .eq. 0) return
  else
    igam   = 1
    icgr   = 0
    iconv  = 1
    ncycle = 10
  end if

  !----------------------------------------------------------------!
  !   Set epsi according to convergence criterion given by iconv   !
  !----------------------------------------------------------------!
  if(iconv.ne.3) then
    if(iconv.eq.4) then
    ama = 0.0
    do i = 1, Amg % imax(1)
      ama = max(ama,a(ia(i)))
    end do
    epsi = epsi*ama
    end if
  else
    fmax = 0.0
    do i = 1, Amg % imax(1)
      fmax = max(fmax,abs(f(i)))
    end do
    epsi = epsi*fmax
  end if

  !------------------------------------!
  !   Decompose Amg % def_relax_down   !
  !------------------------------------!
  if(Amg % def_relax_down .ne. 0) then
    call Amg % Get_Integer_Digits(Amg % def_relax_down,  &
                                  9,                     &
                                  n_digits,              &
                                  Amg % type_relax_down)
    Amg % n_relax_down = Amg % type_relax_down(2)
    Amg % nrdlen       = n_digits-2
    do i = 1, Amg % nrdlen  ! seems to shift the info two places up
      Amg % type_relax_down(i) = Amg % type_relax_down(i+2)
    end do
  else
    Amg % n_relax_down       = 6
    Amg % nrdlen             = 2
    Amg % type_relax_down(1) = 3
    Amg % type_relax_down(2) = 1
  end if

  !----------------------------------!
  !   Decompose Amg % def_relax_up   !
  !----------------------------------!
  if(Amg % def_relax_up .ne. 0) then
    call Amg % Get_Integer_Digits(Amg % def_relax_up,  &
                                  9,                   &
                                  n_digits,            &
                                  Amg % type_relax_up)
    Amg % n_relax_up = Amg % type_relax_up(2)
    Amg % nrulen     = n_digits-2
    do i = 1, Amg % nrulen  ! seems to shift the info two places up
      Amg % type_relax_up(i) = Amg % type_relax_up(i+2)
    end do
  else
    Amg % n_relax_up       = 6
    Amg % nrulen           = 2
    Amg % type_relax_up(1) = 3
    Amg % type_relax_up(2) = 1
  end if

  !---------------------------------------!
  !   Decompose Amg % def_coarse_solver   !
  !---------------------------------------!
  if(Amg % def_coarse_solver .ne. 0) then
    call Amg % Get_Integer_Digits(Amg % def_coarse_solver, 2, n_digits, digit)
    Amg % coarse_solver  = digit(1)
    Amg % n_relax_coarse = digit(2)

    !--------------------------------------------------------------!
    !   In the case coarse grid solution is not obtained with GS   !
    !   sweeps do not use coarsest grid with less than 10 points   !
    !--------------------------------------------------------------!
    if(Amg % coarse_solver .ne. AMG_SOLVER_GS) then
      do i = m, 1, -1
        l = i
        if(Amg % imax(i) - Amg % imin(i).ge.9) exit
      end do
      m = i
      levels = i
    end if
  else
    Amg % coarse_solver  = 1
    Amg % n_relax_coarse = 0
  end if

  !-------------!
  !   Cycling   !
  !-------------!
  if(Amg % iout .ne. 0) then
    call Amg % Calculate_Residual(1, Amg % res0, a, u, f, ia, ja, iw)
    if(Amg % iout .ge. 3) then
      write (6, 9005)
      write (6, 9000) Amg % res0
    end if
    resold = Amg % res0
  end if

  do iter = 1, ncycle
    call Amg % Backup_U(1, icgr, u, m)
    call Amg % One_Cycle(1, igam,              &
                         a, u, f, ia, ja,      &
                         iw,  ifg, icg,        &
                         m, iter, msel, fac, levels)
    call Amg % Cg_Step(1, icgr, iter, a, u, f, ia, ja, iw, m)
    if(Amg % ierr.gt.0) return
    if(iter .eq. 1) then
      mfirst = m
      if(Amg % iout .ge. 3) write(6, 9040) m
    else if(Amg % iout .ge. 3.and.m.ne.mfirst) then
      mfirst = m
      write(6, 9040) m
    end if
    if(Amg % iout .ge. 3 .or. iconv .ne. 1) then
      call Amg % Calculate_Residual(1, Amg % res, a, u, f, ia, ja, iw)
    end if
    Amg % ncyc0 = iter
    if(Amg % iout .ge. 3) then
      cfac = Amg % res / (resold + TINY)
      resold = Amg % res
      if(m .ne. 1) then
        call Amg % Calculate_Residual(2, rescg, a, u, f, ia, ja, iw)
        write(6, 9010) iter, rescg, Amg % res, cfac
      else
        write(6, 9020) iter, Amg % res, cfac
      end if
    end if
    if(iconv .ne. 1) then
      epsil = epsi
      if(iconv.eq.4) then
        umax = 0.0
        do i = Amg % imin(1), Amg % imax(1)
          umax = max(umax,abs(u(i)))
        end do
        epsil = epsi*umax
      end if
    end if
    if(Amg % res .lt. epsil) exit
  end do
  if(Amg % iout .lt. 3 .and. Amg % iout .ne. 0) then
    call Amg % Calculate_Residual(1, Amg % res, a, u, f, ia, ja, iw)
  end if

9000  format (' cycle  0:',3x,'res=',d9.3)
9005  format (/' ************* cycling..... *************'/)
9010  format ( ' cycle ',i2,':   rescg=',d9.3,'   res=',d9.3,  &
               '   cfac=',d9.3)
9020  format ( ' cycle ',i2,':   res=',d9.3,'   cfac=',d9.3)
9040  format (/' cycling between grids 1 and',i3,':'/)

      end subroutine
