!==============================================================================!
  subroutine Solve(Amg, madapt, ncyc, iout,   &
                   a, u, f, ia, ja,           &
                   iw, icg, ifg,              &
                   ncyc0, levels)
!------------------------------------------------------------------------------!
!   Solution phase of Amg1r5
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type)  :: Amg
  integer          :: madapt, ncyc, iout
  double precision :: a(:), u(:), f(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:), icg(:), ifg(:)
  integer          :: ncyc0, levels
!-----------------------------------[locals]-----------------------------------!
  double precision :: ama, cfac, epsi, epsil, fac, fmax
  double precision :: rescg, resold, umax
  integer          :: digit(AMG_MAX_LEVELS)
  integer          :: i, icgr, iconv, igam, iter, l, m, mfirst, msel, n
  integer          :: ncycle, n_digits
  integer          :: ndu
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
  endif
  if(ndu .lt. Amg % imax(levels)) then
    write(6, '(a)')  &
      ' *** error in Solve: ndu too small ***'
    Amg % ierr = AMG_ERR_DIM_F_TOO_SMALL
    return
  endif

  m = levels
  ncyc0 = 0
  do n = 11, 20
    Amg % time(n) = 0.0
  end do
  if(Amg % eps .ne. 0.d0) then
    epsi = Amg % eps
  else
    epsi = 1.d-12
  endif

  !----------------------!
  !   Decompose madapt   !
  !----------------------!
  if (madapt.ne.0) then
    call Amg % Get_Integer_Digits(madapt, 2, n_digits, digit)
    msel = digit(1)
    if (msel.eq.2) then
      if (digit(2).ne.0) then
        fac = dble(digit(2))
        do i = 1, 100
          fac = fac/10.d0
          if(fac .le. 1.d0) exit
        end do
      else
        fac = 0.7d0
      endif
    endif
  else
    msel = 2
    fac = 0.7d0
  endif

  !--------------------!
  !   Decompose ncyc   !
  !--------------------!
  if (ncyc.ne.0) then
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
  endif

  !----------------------------------------------------------------!
  !   Set epsi according to convergence criterion given by iconv   !
  !----------------------------------------------------------------!
  if (iconv.ne.3) then
    if (iconv.eq.4) then
    ama = 0.d0
    do i = 1, Amg % imax(1)
      ama = max(ama,a(ia(i)))
    end do
    epsi = epsi*ama
    endif
  else
    fmax = 0.d0
    do i = 1, Amg % imax(1)
      fmax = max(fmax,abs(f(i)))
    end do
    epsi = epsi*fmax
  endif

  !-------------------------!
  !   Decompose Amg % nrd   !
  !-------------------------!
  if(Amg % nrd .ne. 0) then
    call Amg % Get_Integer_Digits(Amg % nrd, 9, n_digits, Amg % nrdtyp)
    Amg % nrdx   = Amg % nrdtyp(2)
    Amg % nrdlen = n_digits-2
    do i = 1, Amg % nrdlen
      Amg % nrdtyp(i) = Amg % nrdtyp(i+2)
    end do
  else
    Amg % nrdx      = 1
    Amg % nrdlen    = 2
    Amg % nrdtyp(1) = 3
    Amg % nrdtyp(2) = 1
  endif

  !-------------------------!
  !   Decompose Amg % nru   !
  !-------------------------!
  if(Amg % nru .ne. 0) then
    call Amg % Get_Integer_Digits(Amg % nru, 9, n_digits, Amg % nrutyp)
    Amg % nrux   = Amg % nrutyp(2)
    Amg % nrulen = n_digits-2
    do i = 1, Amg % nrulen
      Amg % nrutyp(i) = Amg % nrutyp(i+2)
    end do
  else
    Amg % nrux      = 1
    Amg % nrulen    = 2
    Amg % nrutyp(1) = 3
    Amg % nrutyp(2) = 1
  endif

  !----------------------------!
  !   Decompose Amg % nsolco   !
  !----------------------------!
  if(Amg % nsolco .ne. 0) then
    call Amg % Get_Integer_Digits(Amg % nsolco, 2, n_digits, digit)
    Amg % nsc  = digit(1)
    Amg % nrcx = digit(2)

    !--------------------------------------------------!
    !   In case of yale-smp coarse grid solution, do   !
    !   not use coarsest grid with less than 10 pnts   !
    !--------------------------------------------------!
    if(Amg % nsc .eq. 2) then
      do i = m, 1, -1
        l = i
        if(Amg % imax(i) - Amg % imin(i).ge.9) exit
      end do
      m = i
      levels = i
    endif
  else
    Amg % nsc  = 1
    Amg % nrcx = 0
  endif

  !-------------!
  !   Cycling   !
  !-------------!
  if (iout.ne.0) then
    call Amg % Calculate_Residual(1, Amg % res0, a, u, f, ia, ja, iw)
    if (iout.eq.3) then
      write (6, 9005)
      write (6, 9000) Amg % res0
    endif
    resold = Amg % res0
  endif

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
      if (iout.eq.3) write(6, 9040) m
    elseif (iout.eq.3.and.m.ne.mfirst) then
      mfirst = m
      write(6, 9040) m
    endif
    if(iout .eq. 3 .or. iconv .ne. 1) then
      call Amg % Calculate_Residual(1, Amg % res, a, u, f, ia, ja, iw)
    end if
    ncyc0 = iter
    if(iout .eq. 3) then
      cfac = Amg % res / (resold+1.0d-40)
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
        umax = 0.d0
        do i = Amg % imin(1), Amg % imax(1)
          umax = max(umax,abs(u(i)))
        end do
        epsil = epsi*umax
      end if
    end if
    if(Amg % res .lt. epsil) exit
  end do
  if(iout .ne. 3 .and. iout .ne. 0) then
    call Amg % Calculate_Residual(1, Amg % res, a, u, f, ia, ja, iw)
  end if

9000  format (' cycle  0:',3x,'res=',d9.3)
9005  format (/' ************* cycling..... *************'/)
9010  format ( ' cycle ',i2,':   rescg=',d9.3,'   res=',d9.3,  &
               '   cfac=',d9.3)
9020  format ( ' cycle ',i2,':   res=',d9.3,'   cfac=',d9.3)
9040  format (/' cycling between grids 1 and',i3,':'/)

      end subroutine
