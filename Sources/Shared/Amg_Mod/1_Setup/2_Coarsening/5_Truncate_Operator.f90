!==============================================================================!
  subroutine Truncate_Operator(Amg, level, ntr_local,  &
                               a, ia, ja)           ! defining system
!------------------------------------------------------------------------------!
!   Truncates operator on grid level corresponding to the value of ntr.
!   ntr has to be 0 or 1:
!
!   =0:    pairs of zeroes are removed from coarse grid operators;
!   =1:    coarse grid operators remain unchanged.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: level, ntr_local
  real            :: a(:)
  integer         :: ia(:), ja(:)
!-----------------------------------[locals]-----------------------------------!
  real    :: at
  integer :: i, i1, imn, imx, j, j1, j2, j3, jt, jpos, nna
  logical :: found
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(ntr_local .eq. 1) return

  call Amg % timer_start()

  imn = Amg % imin(level)
  imx = Amg % imax(level)
  nna = ia(imx+1)-ia(imn)
  jpos = ia(imn)

  do i = imn, imx
    j1 = ia(i)+1
    j2 = ia(i+1)-1
    a(jpos)  = a(ia(i))
    ja(jpos) = i
    ia(i) = jpos
    jpos = jpos+1
    do j = j1, j2
      i1 = ja(j)
      if(i1 .ge. 0) then
        if(i1.ge.i .and.a(j).eq.0.0) then
          found = .false.
          do j3 = ia(i1), ia(i1+1) - 1
            if(ja(j3) .eq. i) then
              jt = j3
              at = a(j3)
              found = .true.
              exit
            end if
          end do

          ! Error message
          if(.not. found) then
            write(6, '(a,a,i3,a)')                               &
              ' *** error in trunc: transpose a-entry missing',  &
              ' on grid', level, ' ***'
            Amg % ierr = AMG_ERR_TRANSPOSE_MISSING
            return
          end if

          if(at.eq.0.0) then
            ja(jt) = -ja(jt)
            ! I am not happy with this cycle here, but it
            ! lesser evil than GOTO which was here before
            cycle
          end if
        end if
        a(jpos)  = a(j)
        ja(jpos) = ja(j)
        jpos = jpos+1
      end if
    end do
  end do
  ia(imx+1) = jpos

  call Amg % timer_stop(7)

  end subroutine
