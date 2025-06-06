!==============================================================================!
  subroutine interpolation_weights(amg, level, ichk, mmax,  &  !  4
                                   a, ia, ja,               &  !  7
                                   iw, ifg, icg,            &  ! from "kwork"
                                   ncolor,                  &  ! from "iwork"
                                   ncolx, iajas)
!------------------------------------------------------------------------------!
!     Set up final coarser grid k+1 and interpolation formula from
!     grid k+1 to grid k. This is the version as described in ruge/
!     stueben (bristol). interpolation_weights assumes the grid to be
!     pre-colored by subroutine pre_color.
!
!     On exit, a, ja and iw are set to contain the interpolation
!     weights and corresponding pointers as required in the solution
!     phase of amg1r5. Also, icg(i) (imin(k)<=i<=imax(k)) are set to
!     their final values, except for those i with icg(i)<0.
!
!     Comments on input
!
!     - ja(ia(i)) ----- i=imin(k),...,imax(k)
!
!     Assumed to point to the last strong connection of point i (or to
!     ia(i) if there is no such connection). on exit, reset to original
!     values (i.e. ja(ia(i))=i).
!
!     - icg(i) -------- i=imin(k),...,imax(k)
!
!     Assumed to contain information on pre-coloring:
!
!       icg(i)>0: i is c-point
!       icg(i)=0: i is forced f-point (i.e. i has no strong connnection)
!       icg(i)<0: i is f-point with (at least) one strong connection.
!
!     Comments on work space used
!
!     - ifg(i) -------- i=imin(k),...,imax(k)
!
!     Is used for several purposes. in particular, to distinguish
!     interpolatory and non-interpolatory points.
!
!     - ncolor(i) ----- i=1,...,#points on grid k
!
!     Is set to f-point-colors to be used later in subroutine opdfn.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: level
  integer          :: ichck
  double precision :: a(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:), ifg(:), icg(:), ncolor(:)
  integer          :: ncolx, iajas
!-----------------------------------[locals]-----------------------------------!
  integer :: i, ic, icgp, ichk, ihi, ilo, ip, is
  integer :: ii, ijas, iblck, iblck1
  integer :: j, jj, jhi, jlo, j1, j2, jw0, jwx, jwpos, jjhi, jjlo
  integer :: mmax, nc, ncondc, ncount, npts, nptsc, ndaja
  double precision :: ewt2i, s, si, scale, ww
  integer :: nda, ndu, ndicg
  logical :: skip
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  nda   = size(a,   1)
  ndu   = size(ia,  1)
  ndicg = size(icg, 1)

  call amg % timer_start()
  ndaja = nda
  ncount = 0
  if (level.eq.1) then
    amg % iminw(1) = 1
    iw(1) = ia(amg % imax(1)+1)
  endif

  ilo = amg % imin(level)
  ihi = amg % imax(level)
  do i = ilo, ihi
    ifg(i) = 0
  end do

  !---------------------------------------------------------------------!
  !                                                                     !
  !   Sweep over f-points i which have at least one strong connection   !
  !                                                                     !
  !---------------------------------------------------------------------!
  iblck = amg % iminw(level)
  outer: do i = ilo, ihi
    if(icg(i) .lt. 0) then
      jlo = ia(i)+1
      jhi = ja(ia(i))
      ewt2i = amg % ewt2 / a(jlo)
      ncondc = 0

      !----------------------------------------------------------!
      !   Initialize "block" of interpolation weights of point   !
      !----------------------------------------------------------!
      inner: do
        jw0 = iw(iblck)
        jwx = jw0
        if(jwx .gt. ndaja) then
          call error_pwint(nda)
          return
        end if
        do j = jlo, jhi
          ii = ja(j)
          if(icg(ii) .le. 0) cycle
          a(jwx)  = a(j)
          ja(jwx) = ii
          ifg(ii) = jwx
          jwx = jwx+1
          if(jwx .gt. ndaja) then
            call error_pwint(nda)
            return
          end if
        end do
        a(jwx)  = a(ia(i))
        ja(jwx) = i
        amg % mda  = max(amg % mda, jwx + iajas)
        ifg(i) = jwx

        !---------------------------------------------------------------------!
        !   Sweep over strongly connected f-points. these must be "covered"   !
        !   by a total weight defined by ewt2. if an f-point has no strong    !
        !   connections, regard it to be covered, but do not distribute       !
        !   the corresponding weight (error at such a point can be assumed    !
        !   to be very small!).                                               !
        !---------------------------------------------------------------------!
        do j = jlo, jhi
          ii = ja(j)
          if(icg(ii) .lt. 0) then

            !-------------------------------------------------------!
            !   Compute dependence on set of interpolation points   !
            !-------------------------------------------------------!
            s  = 0.0d0
            si = 0.0d0
            jjlo = ia(ii)+1
            jjhi = ia(ii+1)-1
            do jj = jjlo, jjhi
              if(ifg(ja(jj)) .lt. jw0) cycle
              if(ja(jj) .eq. i) si = a(jj)
              s = s + a(jj)
            end do
            skip = .false.
            if(ichk .ne. 2) then
              if(s .eq. 0.0d0) then
                a(jwx) = a(jwx)+a(j)
                cycle
              else
                skip = .true.
              end if
            end if
            if(.not. skip) then

              !-----------------------------------------------------!
              !   Check dependence on set of interpolation points   !
              !-----------------------------------------------------!
              if(s-si .gt. ewt2i*a(j)*a(jjlo)) then

                ! Dependence too small: if there is not yet a conditional
                ! c-point, make ii such a point and restart the
                ! process for defining interpolation weights for point i.
                ! otherwise make i itself a c-point and leave ii an f-point.
                if(ncondc .eq. 0) then
                  ncount = ncount+1
                  ncondc = 1
                  ip = ii
                  icgp = icg(ii)
                  icg(ii) = 1
                  cycle inner
                else
                  icg(i) = 1
                  icg(ip) = icgp
                  do jj = jw0, jwx
                    ifg(ja(jj)) = 0
                  end do
                  cycle outer
                end if
              end if
            end if

            !---------------------------------------!
            !   Distribute the weight of point ii   !
            !---------------------------------------!
            ww = a(j)/s
            do jj = jjlo, jjhi
              if(ifg(ja(jj)).ge.jw0)  &
                a(ifg(ja(jj))) = a(ifg(ja(jj)))+a(jj)*ww
            end do
          end if

        end do
        exit inner
      end do inner

      !--------------------------------------------------------------!
      !   All necessary points are covered. now distribute weights   !
      !   from weak connections of point i (analogous as above)      !
      !--------------------------------------------------------------!
      do j = jhi + 1, ia(i+1) - 1
        ii = ja(j)
        if(icg(ii) .eq. 0) cycle
        s = 0.0d0
        jjlo = ia(ii)+1
        jjhi = ia(ii+1)-1
        do jj = jjlo, jjhi
          if (ifg(ja(jj)).ge.jw0) s = s+a(jj)
        end do
        if (s.eq.0.0d0) then
          a(jwx) = a(jwx)+a(j)
        else
          ww = a(j)/s
          do jj = jjlo, jjhi
            if (ifg(ja(jj)).ge.jw0)  &
               a(ifg(ja(jj))) = a(ifg(ja(jj)))+a(jj)*ww
          end do
        endif
      end do

      icg(i) = -iblck
      iblck = iblck+1
      iw(iblck) = jwx+1
    end if
  end do outer

  !-----------------------------------------------------------!
  !   Set icg; reset ja(ia(i)); check size of coarsest grid   !
  !-----------------------------------------------------------!
  ic = ihi
  amg % imin(level+1) = ic+1
  do i = ilo, ihi
    ja(ia(i)) = i
    if(icg(i) .le. 0) cycle
    ic = ic+1
    icg(i) = ic
    if(ic .lt. ndicg) icg(ic) = 0
  end do
  amg % imax(level+1) = ic

  npts = ihi-ilo+1
  nptsc = amg % imax(level+1) - amg % imin(level+1)+1
  if(nptsc.eq.1)                                             mmax = level+1
  if(nptsc.eq.1 .and. amg % irow0 .eq. AMG_SINGULAR_MATRIX)  mmax = level
  if(nptsc.eq.npts.or.nptsc.eq.0)                            mmax = level
  if(level .lt. mmax) then
    if(ic .ge. ndu) then
      write(6, '(a)')  &
        ' *** error in interpolation_weights: ndu too small ***'
      amg % ierr = AMG_ERR_DIM_IA_TOO_SMALL
      return
    end if
    if(ic .ge. ndicg) then
      write(6, '(a)')  &
        ' *** error in interpolation_weights: ndw too small ***'
      amg % ierr = AMG_ERR_DIM_ICG_TOO_SMALL
      return
    end if
  end if

  !------------------!
  !   Re-arrange a   !
  !------------------!
  iblck1 = amg % iminw(level)
  jwpos = iw(iblck1)
  do i=ilo,ihi
    if(icg(i) .lt. 0) then
      iblck = -icg(i)
      j1 = iw(iblck)
      j2 = iw(iblck+1)-1
      if (j2.le.j1) then
        icg(i) = 0
      else
        icg(i) = -iblck1
        iw(iblck1) = jwpos
        scale = -1.0d0/a(j2)
        do j = j1, j2 - 1
          a(jwpos)  = a(j)*scale
          ja(jwpos) = icg(ja(j))
          jwpos = jwpos+1
        end do
        iblck1 = iblck1+1
      end if
    end if
  end do
  amg % imaxw(level) = iblck1-1
  iw(iblck1) = jwpos

  !--------------------------------------------------------------!
  !   Store type of points (i.e. c, f or ff)  on vector ncolor   !
  !--------------------------------------------------------------!
  is = 1-ilo
  do i = ilo, ihi
    nc = icg(i)
    if (nc.lt.0) then
      ncolor(i+is) = 1
    elseif (nc.gt.0) then
      ncolor(i+is) = 2
    else
      ncolor(i+is) = 3
    endif
  end do
  ncolx = 1

  !----------!
  !   Exit   !
  !----------!
# ifdef VERBOSE
    write(6, '(a,i3,a,i4)')   &
      ' interpolation operator no.', level,  &
      ' completed. c-points added in interpolation_weights:', ncount
# endif

  call amg % timer_stop(4)

  end subroutine

!==============================================================================!
  subroutine error_pwint(nda)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  integer nda
!==============================================================================!

  write(6, '(a)')  &
    ' *** error in interpolation_weights: nda too small ***'
  amg % ierr = AMG_ERR_DIM_JA_TOO_SMALL

  end subroutine
