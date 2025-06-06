!==============================================================================!
  subroutine solve(amg, madapt, ncyc, iout,   &
                   a, u, f, ia, ja,           &
                   iw, icg, ifg,              &
                   ncyc0, levels)
!------------------------------------------------------------------------------!
!   Solution phase of amg1r5
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
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
  if(ndu .lt. amg % imax(levels)) then
    write(6, '(a)')  &
      ' *** error in solve: ndu too small ***'
    amg % ierr = AMG_ERR_DIM_U_TOO_SMALL
    return
  endif
  if(ndu .lt. amg % imax(levels)) then
    write(6, '(a)')  &
      ' *** error in solve: ndu too small ***'
    amg % ierr = AMG_ERR_DIM_F_TOO_SMALL
    return
  endif

  m = levels
  ncyc0 = 0
  do n = 11, 20
    amg % time(n) = 0.0
  end do
  if(amg % eps .ne. 0.d0) then
    epsi = amg % eps
  else
    epsi = 1.d-12
  endif

  !----------------------!
  !   Decompose madapt   !
  !----------------------!
  if (madapt.ne.0) then
    call amg % get_integer_digits(madapt, 2, n_digits, digit)
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
    call amg % get_integer_digits(abs(ncyc), 4, n_digits, digit)
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
    do i = 1, amg % imax(1)
      ama = max(ama,a(ia(i)))
    end do
    epsi = epsi*ama
    endif
  else
    fmax = 0.d0
    do i = 1, amg % imax(1)
      fmax = max(fmax,abs(f(i)))
    end do
    epsi = epsi*fmax
  endif

  !-------------------------!
  !   Decompose amg % nrd   !
  !-------------------------!
  if(amg % nrd .ne. 0) then
    call amg % get_integer_digits(amg % nrd, 9, n_digits, amg % nrdtyp)
    amg % nrdx   = amg % nrdtyp(2)
    amg % nrdlen = n_digits-2
    do i = 1, amg % nrdlen
      amg % nrdtyp(i) = amg % nrdtyp(i+2)
    end do
  else
    amg % nrdx      = 1
    amg % nrdlen    = 2
    amg % nrdtyp(1) = 3
    amg % nrdtyp(2) = 1
  endif

  !-------------------------!
  !   Decompose amg % nru   !
  !-------------------------!
  if(amg % nru .ne. 0) then
    call amg % get_integer_digits(amg % nru, 9, n_digits, amg % nrutyp)
    amg % nrux   = amg % nrutyp(2)
    amg % nrulen = n_digits-2
    do i = 1, amg % nrulen
      amg % nrutyp(i) = amg % nrutyp(i+2)
    end do
  else
    amg % nrux      = 1
    amg % nrulen    = 2
    amg % nrutyp(1) = 3
    amg % nrutyp(2) = 1
  endif

  !----------------------------!
  !   Decompose amg % nsolco   !
  !----------------------------!
  if(amg % nsolco .ne. 0) then
    call amg % get_integer_digits(amg % nsolco, 2, n_digits, digit)
    amg % nsc  = digit(1)
    amg % nrcx = digit(2)

    !--------------------------------------------------!
    !   In case of yale-smp coarse grid solution, do   !
    !   not use coarsest grid with less than 10 pnts   !
    !--------------------------------------------------!
    if(amg % nsc .eq. 2) then
      do i = m, 1, -1
        l = i
        if(amg % imax(i) - amg % imin(i).ge.9) exit
      end do
      m = i
      levels = i
    endif
  else
    amg % nsc  = 1
    amg % nrcx = 0
  endif

  !-------------!
  !   Cycling   !
  !-------------!
  if (iout.ne.0) then
    call amg % calculate_residual(1, amg % res0, a, u, f, ia, ja, iw)
    if (iout.eq.3) then
      write (6, 9005)
      write (6, 9000) amg % res0
    endif
    resold = amg % res0
  endif

  do iter = 1, ncycle
    call amg % backup_u(1, icgr, u, m)
    call amg % one_cycle(1, igam,              &
                         a, u, f, ia, ja,      &
                         iw,  ifg, icg,        &
                         m, iter, msel, fac, levels)
    call amg % cg_step(1, icgr, iter, a, u, f, ia, ja, iw, m)
    if(amg % ierr.gt.0) return
    if(iter .eq. 1) then
      mfirst = m
      if (iout.eq.3) write(6, 9040) m
    elseif (iout.eq.3.and.m.ne.mfirst) then
      mfirst = m
      write(6, 9040) m
    endif
    if(iout .eq. 3 .or. iconv .ne. 1) then
      call amg % calculate_residual(1, amg % res, a, u, f, ia, ja, iw)
    end if
    ncyc0 = iter
    if(iout .eq. 3) then
      cfac = amg % res / (resold+1.0d-40)
      resold = amg % res
      if(m .ne. 1) then
        call amg % calculate_residual(2, rescg, a, u, f, ia, ja, iw)
        write(6, 9010) iter, rescg, amg % res, cfac
      else
        write(6, 9020) iter, amg % res, cfac
      end if
    end if
    if(iconv .ne. 1) then
      epsil = epsi
      if(iconv.eq.4) then
        umax = 0.d0
        do i = amg % imin(1), amg % imax(1)
          umax = max(umax,abs(u(i)))
        end do
        epsil = epsi*umax
      end if
    end if
    if(amg % res .lt. epsil) exit
  end do
  if(iout .ne. 3 .and. iout .ne. 0) then
    call amg % calculate_residual(1, amg % res, a, u, f, ia, ja, iw)
  end if

9000  format (' cycle  0:',3x,'res=',d9.3)
9005  format (/' ************* cycling..... *************'/)
9010  format ( ' cycle ',i2,':   rescg=',d9.3,'   res=',d9.3,  &
               '   cfac=',d9.3)
9020  format ( ' cycle ',i2,':   res=',d9.3,'   cfac=',d9.3)
9040  format (/' cycling between grids 1 and',i3,':'/)

      end subroutine
