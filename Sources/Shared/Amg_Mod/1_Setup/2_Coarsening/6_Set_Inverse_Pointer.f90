!==============================================================================!
  subroutine set_inverse_pointer(amg, icg, ifg, levels)
!------------------------------------------------------------------------------!
!   Set "inverse" pointer ifg
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: icg(:), ifg(:)
  integer          :: levels
!-----------------------------------[locals]-----------------------------------!
  integer :: i, ib, ist, level
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  call amg % timer_start()

  do i = amg % imin(1), amg % imax(levels-1)
    if(icg(i) .gt. 0) ifg(icg(i)) = i
  end do

  ib = 1
  do level = 1, levels - 1
    ist = amg % nstcol(level)
    if(ist .lt. AMG_BIG_INTEGER) then
      do
        ifg(ib) = ist
        ist = -icg(ist)
        ib = ib + 1
        if(ist .ge. AMG_BIG_INTEGER) exit
      end do
    end if
  end do

  call amg % timer_stop(8)

  end subroutine
