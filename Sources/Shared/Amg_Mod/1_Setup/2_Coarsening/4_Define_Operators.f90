!==============================================================================!
  subroutine Define_Operators(Amg, level, mmax,  &
                              a, ia, ja,         &   !  linear system
                              iw, icg, ifg,      &   !  work arrays
                              ncolor, ncolx, iajas)
!------------------------------------------------------------------------------!
!     This subroutine constructs the cg-operator a(level) (level>1).
!
!     one row is constructed at a time, and the row structure
!     (i.e., which connection goes where) must be determined
!     at the same time.  in order to avoid searches through
!     the current row to determine if a position for a connection
!     has already been defined, icg (for level level) and ifg (for level
!     level-1) are used auxiliarily. furthermore, to speed up computation,
!     the transpose of interpolation is temporarily stored on a (at
!     the end of the available space). ifg (for level level) is used
!     as pointer to the corresponding rows. more details: see below.
!
!     Comments on input
!
!     - imin(level), imax(level)
!
!     first/last number of grid points on grid level.
!
!     - iminw(level-1), imaxw(level-1)
!
!     first/last "block" number of interpolation weights contributing
!     in interpolation to grid level-1.
!
!     - iw(ibl) ----- ibl=iminw(level-1),imaxw(level-1)+1
!     - a(j) -------- j=jlo,jhi      (jlo=iw(iminw(level-1)))
!     - ja(j) ------- j=jlo,jhi      (jhi=iw(imaxw(level-1)+1)-1)
!
!     are assumed to contain the weights of interpolation along with
!     the necessary pointers. iw(ibl) points to the first entry of
!     "block" number ibl in a (cf. below).
!
!     - icg(if) ----- if=imin(level-1),imax(level-1)
!
!     are assumed to be set to their final values:
!
!       icg(if)>0: if is c-point of grid level-1 and icg(if) points just
!                  to the corresponding point on grid level;
!       icg(if)=0: if is f-point of grid level-1 without any contribution
!                  in interpolation from grid level;
!       icg(if)<0: if is f-point of grid level-1; ibl=-icg(if) points just
!                  to the "block" of interpolation weights. that means,
!                  jw(j) (iw(ibl)<=j<=iw(ibl+1)-1) points to the points
!                  (on grid level) which contribute in interpolation to if.
!                  the corresponding weights are stored in a(j).
!
!     - ifg(if) ----- if=imin(level-1),imax(level-1)
!
!     work space (see below). ifg(if) is assumed to be > -imin(level).
!
!     - ifg(ic) ----- ic=imin(level),imax(level)+1
!
!     work space (see below). the contents of ifg(ic) is arbitrary.
!
!     - icg(ic) ----- ic=imin(level),imax(level)
!
!     work space (see below). icg(ic) is assumed to be zero.
!
!     Comments on work space used
!
!     - a(j) -------- j=jjlo,jjhi
!     - ja(j) ------- j=jjlo,jjhi
!
!     (jjlo=min(nda,nda)-iw(imaxw(level-1)+1)+iw(iminw(level-1)+1),
!      jjhi=min(nda,nda), i.e. the space of length of interpolation
!      (level-1) at the end of a and ja, respectively.)
!
!     is used to store the transpose of interpolation w(level-1). this is
!     done in order to speed up the computation of the operator a(level).
!     if during assemblage of a(level) this work space is required by a(level),
!     processing con tinues using the not yet redefined part of the work
!     space and recalculating the lost entries each time they are
!     needed. thus the calculation slows down.
!
!     - ifg(ic) ----- ic=imin(level),imax(level)+1
!
!     is used as pointer for the transpose of interpolation: the coarse-
!     grid point ic contributes in interpolation to the fine-grid points
!     if=ja(j) (ifg(ic)<=j<=ifg(ic+1)-1) with weight a(j). the
!     contribution to itself (by the weight 1.0) is not contained.
!
!     - icg(ic) ----- ic=imin(level),imax(level)
!
!     is used for several purposes. in assembling the cg-operator a(level),
!     it serves as pointer to positions in a(level) which have already been
!     defined: if the current row corresponds to point ic1, and a
!     connection to another cg-point ic2 has just been found, then if
!     icg(ic2)<ia(ic1), the corresponding entry in row ic1 of a(level)
!     has not yet been defined. otherwise, icg(ic2) points to the
!     location for that entry. (also see ifg below.)
!
!     - ifg(if) ----- if=imin(level-1),imax(level-1)
!
!     In assembling the cg-operator a(level), this vector contains informat-
!     ion on whether the existence of entries in a(level) has to be checked
!     or not.
!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: level, mmax
  real            :: a(:)
  integer         :: ia(:), ja(:)
  integer         :: iw(:), icg(:), ifg(:), ncolor(:)
  integer         :: ncolx
  integer         :: iajas
!-----------------------------------[locals]-----------------------------------!
  real    :: wjf1, ww
  integer :: i, iadr, ialow, iadrs, ibl, ic, ic1, ic2, ic3, icol,   &
             if, if1, if2, ihi, ilo, ilo1, ist, istti, j, jb, jc3,  &
             jf1, jf2, jf3, jpos, k1, ndaja, npts
  integer :: nda
  logical :: found
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  nda = size(a, 1)

  !------------------------------------------------------!
  !   Extend a, ja to store transpose of interpolation   !
  !------------------------------------------------------!
  call Amg % timer_start()
  ndaja = nda
  Amg % iminw(level) = Amg % imaxw(level-1)+1
  if(level .le. mmax) then
    jpos = iw(Amg % iminw(level))
    do j = iw(Amg % iminw(level-1)), iw(Amg % imaxw(level-1)+1)-1
      icg(ja(j)) = icg(ja(j))+1
    end do

    Amg % mda = max(Amg % mda, jpos + jpos - iw(Amg % iminw(level-1))+iajas)
    ifg(Amg % imin(level)) = ndaja - jpos + iw(Amg % iminw(level-1))+1
    if(ifg(Amg % imin(level)).le.jpos) then
      call error_9900(nda)
      return
    end if
    do ic = Amg % imin(level), Amg % imax(level)
      ifg(ic+1) = ifg(ic)+icg(ic)
      icg(ic) = ifg(ic)
    end do

    ! Watch out: if is the loop variable
    do if = Amg % imin(level-1), Amg % imax(level-1)
      if(icg(if).lt.0) then
        ibl = -icg(if)
        do j = iw(ibl), iw(ibl+1) - 1
          ic = ja(j)
          a(icg(ic))=a(j)
          ja(icg(ic))=if
          icg(ic)=icg(ic)+1
        end do
      end if
    end do

    do ic = Amg % imin(level), Amg % imax(level)
      icg(ic) = 0
    end do

    !---------------------------------------------------------------!
    !   Sweep over all cg-points ic to assemble rows of cg matrix   !
    !---------------------------------------------------------------!
    istti = 1
    do if = Amg % imin(level-1), Amg % imax(level-1)
      if (icg(if).le.0) cycle
      ic = icg(if)
      if(jpos .gt. ndaja) then
        call error_9900(nda)
        return
      end if
      ialow = jpos
      icg(ic) = jpos
      a(jpos)  = a(ia(if))
      ja(jpos) = ic
      jpos = jpos+1

      !-------------------------------------------------!
      !                                                 !
      !   Search for c-c-c-c and c-c-f-c connections    !
      !                                                 !
      !-------------------------------------------------!
      do jf1 = ia(if)+1, ia(if+1)-1
        if1 = ja(jf1)
        ic1 = icg(if1)

        !----------------------------------------------------!
        !   if1 is f-point: sweep over c-c-f-c connections   !
        !----------------------------------------------------!
        if(ic1.lt.0) then
          ifg(if1) = -ic
          do jf2 = iw(-ic1), iw(-ic1+1)-1
            ic2 = ja(jf2)
            if(icg(ic2) .ge. ialow) then
              a(icg(ic2)) = a(icg(ic2))+a(jf1)*a(jf2)
            else
              if(jpos .gt. ndaja) then
                call error_9900(nda)
                return
              end if
              icg(ic2) = jpos
              a(jpos)  = a(jf1)*a(jf2)
              ja(jpos) = ic2
              jpos = jpos+1
            end if
          end do
        end if

        !----------------------------------------!
        !   if1 is c-point: c-c-c-c connection   !
        !----------------------------------------!
        if(ic1.gt.0) then
          if(icg(ic1) .lt. ialow) then
            if(jpos .gt. ndaja) then
              call error_9900(nda)
              return
            end if
            icg(ic1) = jpos
            a(jpos)  = a(jf1)
            ja(jpos) = ic1
            jpos = jpos+1
          else
            a(icg(ic1)) = a(icg(ic1))+a(jf1)
          end if
        end if
      end do
      ist = Amg % imin(level-1)-1

      !------------------------------------------------!
      !                                                !
      !   Search for c-f-c-c and c-f-f-c connections   !
      !                                                !
      !------------------------------------------------!
      do jf1 = ifg(ic), ifg(ic+1) - 1
        if (jf1.ge.jpos) then
          if1 = ja(jf1)
          wjf1 = a(jf1)
          ist = if1
        else
          istti = 0
          found = .false.
          do jf2 = ist + 1, Amg % imax(level-1)
            if(icg(jf2) .lt. 0) then
              ibl = -icg(jf2)
              do jb = iw(ibl), iw(ibl+1) - 1
                jc3 = ja(jb)
                if(jc3 .eq. ic) then
                  if1 = jf2
                  wjf1 = a(jb)
                  found = .true.
                  exit
                end if
              end do
              if(found) exit
            end if
          end do

          !----------------!
          !   Error exit   !
          !----------------!
          if(.not. found) then
            write(6, '(a,a,i3)')                                    &
              ' *** error in opdfn: interpolation entry missing ',  &
              'on grid ', level-1
            Amg % ierr = AMG_ERR_INTERP_MISSING
            return
          end if
          ist = if1
        endif
        do jf2 = ia(if1), ia(if1+1) - 1
          if2 = ja(jf2)
          ic2 = icg(if2)
          ww = wjf1*a(jf2)
          if(ic2 .ne. 0) then
            if(ic2 .lt. 0) then

              !----------------------------------------------------!
              !   if2 is f-point: sweep over c-f-f-c connections   !
              !----------------------------------------------------!
              if(ifg(if2) .ne. -ic) then
                ifg(if2) = -ic
                do jf3 = iw(-ic2), iw(-ic2+1) - 1
                  ic3 = ja(jf3)
                  if(icg(ic3) .ge. ialow) then
                    a(icg(ic3)) = a(icg(ic3))+ww*a(jf3)
                  else
                    if(jpos .gt. ndaja) then
                      call error_9900(nda)
                      return
                    end if
                    icg(ic3) = jpos
                    a(jpos) = ww*a(jf3)
                    ja(jpos) = ic3
                    jpos = jpos+1
                 end if
                end do

             !-------------------------------------------------------------!
             !   if2 has been encountered before; do not check positions   !
             !-------------------------------------------------------------!
              else
                do jf3 = iw(-ic2), iw(-ic2+1) - 1
                  iadrs = icg(ja(jf3))
                  a(iadrs) = a(iadrs)+ww*a(jf3)
                end do
              end if
            end if

            !----------------------------------------!
            !   if2 is c-point: c-f-c-c connection   !
            !----------------------------------------!
            if(ic2 .gt. 0) then
              if(icg(ic2) .ge. ialow) then
                a(icg(ic2)) = a(icg(ic2))+ww
              else
                if(jpos .gt. ndaja) then
                  call error_9900(nda)
                  return
                end if
                icg(ic2) = jpos
                a(jpos) = ww
                ja(jpos) = ic2
                jpos = jpos+1
              end if
            end if
          end if
        end do
      end do
      ia(ic+1) = jpos
    end do
    Amg % mda = max(Amg % mda, ia(Amg % imax(level)+1)-1 + iajas)

    !----------------------------------!
    !   Warning for storage shortage   !
    !----------------------------------!
    if (istti.ne.1) then
      k1 = level-1
      write(6, '(a,a,i2,a,a)')                                     &
        ' --- warng: unable to store transpose of interpolation',  &
        ' on grid ',k1,' during execution of opdfn, because nda',  &
        ' or nda too small. setup computation is slowing down.'
      Amg % ierr = AMG_WARN_YALE_STORAGE_A
    endif
#   ifdef VERBOSE
      write(6, '(a,i3,a)')  &
        ' coarse grid operator no.', level, ' completed'
#   endif

  !-------------------------------------------------------!
  !   Set up linked list for relaxation on grid level-1   !
  !-------------------------------------------------------!
  end if
  ia(Amg % imin(level)) = iw(Amg % iminw(level))
  ilo = Amg % imin(level-1)
  ihi = Amg % imax(level-1)
  ilo1 = ilo-1
  npts = ihi-ilo1
  ist = AMG_BIG_INTEGER
  do icol = 1, ncolx
    do i = npts, 1, -1
      if(ncolor(i) .ne. icol) cycle
      icg(i+ilo1) = -ist
      ist=i+ilo1
    end do
  end do
  Amg % nstcol(level-1) = ist

  !---------------------------!
  !   Exit / error messages   !
  !---------------------------!

  call Amg % timer_stop(6)

  end subroutine

!==============================================================================!
  subroutine error_9900(nda)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  integer :: nda
!==============================================================================!

  write(6, '(a)')  &
    ' *** error in opdfn: nda too small ***'
  Amg % ierr = AMG_ERR_DIM_JA_TOO_SMALL

  end subroutine
