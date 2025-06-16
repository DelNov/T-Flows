!==============================================================================!
  subroutine Set_Inverse_Pointer(Amg, icg, ifg, levels)
!------------------------------------------------------------------------------!
!   Set "inverse" pointer ifg
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type)  :: Amg
  integer          :: icg(:), ifg(:)
  integer          :: levels
!-----------------------------------[locals]-----------------------------------!
  integer :: i, ib, ist, level
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  do i = Amg % imin(1), Amg % imax(levels-1)
    if(icg(i) .gt. 0) ifg(icg(i)) = i
  end do

  ib = 1
  do level = 1, levels - 1
    ist = Amg % start_of_color(level)
    if(ist .lt. AMG_BIG_INTEGER) then
      do
        ifg(ib) = ist
        ist = -icg(ist)
        ib = ib + 1
        if(ist .ge. AMG_BIG_INTEGER) exit
      end do
    end if
  end do

  end subroutine
