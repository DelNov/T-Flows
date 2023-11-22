!==============================================================================!
  subroutine Print_Regions_List(Grid)
!------------------------------------------------------------------------------!
!   Prints a list of regions (boundary conditions) in a grid.                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer                   :: j, ll, lc, ln
  character(:), allocatable :: reg_list
!==============================================================================!

  !----------------------------------------------------------!
  !   Work out the list length and allocate memeory for it   !
  !----------------------------------------------------------!
  ln = 0
  do j = 1, Grid % n_bnd_regions
    ln = ln + len_trim(Grid % region % name(j)) + 10
  end do
  Assert(ln > 0)
  allocate(character(ln) :: reg_list)
  reg_list(:) = ' '                     ! clear the list

  !--------------------------------------!
  !   Form the boundary condition list   !
  !--------------------------------------!
  ll = 1
  do j = 1, Grid % n_bnd_regions
    if(j < Grid % n_bnd_regions) then
      lc = len_trim(Grid % region % name(j)) + 10
      write(reg_list(ll:ll+lc), '(i4,a2,a,a4)')  &
            j,                                   &
            '. ',                                &
            trim(Grid % region % name(j)),       &
            ' \n '
    else
      lc = len_trim(Grid % region % name(j)) +  6
      write(reg_list(ll:ll+lc), '(i4,a2,a,a4)')  &
            j,                                   &
            '. ',                                &
            trim(Grid % region % name(j))
    end if
    ll = ll + lc
  end do

  !-------------------------------!
  !   Call the Message function   !
  !-------------------------------!
  call Message % Framed(20,                                   &
    "Grid currently has the following boundary conditions:",  &
    reg_list(1:ll))

  end subroutine
