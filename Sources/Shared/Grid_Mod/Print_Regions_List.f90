!==============================================================================!
  subroutine Print_Regions_List(Grid)
!------------------------------------------------------------------------------!
!   Prints a list of regions (boundary conditions) in a grid.                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer                 :: j, ll, lc
  character(MAX_TOKENS*2) :: reg_list
!==============================================================================!

  ! Clear the list
  reg_list(:) = ' '

  !--------------------------------------!
  !   Form the boundary condition list   !
  !--------------------------------------!
  ll = 1
  do j = 1, Grid % n_bnd_regions
    if(j < Grid % n_bnd_regions) then
      lc = len_trim(Grid % region % name(j)) + 8
      write(reg_list(ll:ll+lc), '(i2,a2,a,a4)')  &
            j,                                   &
            '. ',                                &
            trim(Grid % region % name(j)),       &
            ' \n '
    else
      lc = len_trim(Grid % region % name(j)) + 4
      write(reg_list(ll:ll+lc), '(i2,a2,a,a4)')  &
            j,                                   &
            '. ',                                &
            trim(Grid % region % name(j))
    end if
    ll = ll + lc + 1
  end do

  !-------------------------------!
  !   Call the Message function   !
  !-------------------------------!
  call Message % Framed(20,                                   &
    "Grid currently has the following boundary conditions:",  &
    reg_list(1:ll))

  end subroutine
