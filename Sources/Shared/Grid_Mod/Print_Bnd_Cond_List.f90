!==============================================================================!
  subroutine Print_Bnd_Cond_List(Grid)
!------------------------------------------------------------------------------!
!   Prints a list of boundary conditions in a grid.                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer                 :: j, ll, lc
  character(MAX_TOKENS*2) :: bnd_cond_list
!==============================================================================!

  ! Clear the list
  bnd_cond_list(:) = ' '

  !--------------------------------------!
  !   Form the boundary condition list   !
  !--------------------------------------!
  ll = 1
  do j = 1, Grid % n_bnd_cond
    if(j < Grid % n_bnd_cond) then
      lc = len_trim(Grid % bnd_cond % name(j)) + 8
      write(bnd_cond_list(ll:ll+lc), '(i2,a2,a,a4)')  &
            j,                                        &
            '. ',                                     &
            trim(Grid % bnd_cond % name(j)),          &
            ' \n '
    else
      lc = len_trim(Grid % bnd_cond % name(j)) + 4
      write(bnd_cond_list(ll:ll+lc), '(i2,a2,a,a4)')  &
            j,                                        &
            '. ',                                     &
            trim(Grid % bnd_cond % name(j))
    end if
    ll = ll + lc + 1
  end do

  !-------------------------------!
  !   Call the Message function   !
  !-------------------------------!
  call Message % Print_Framed_Text(20,                        &
    "Grid currently has the following boundary conditions:",  &
    bnd_cond_list(1:ll))

  end subroutine
