!==============================================================================!
  subroutine Print_Regions_List(Grid)
!------------------------------------------------------------------------------!
!>  This subroutine compiles and displays a list of boundary condition regions
!>  defined in the computational grid. It's typically used during pre-processing
!>  stages such as in the Generate and Convert programs to verify the setup of
!>  boundary conditions within the grid.
!------------------------------------------------------------------------------!
!   Process                                                                    !
!                                                                              !
!   * List Length Calculation:                                                 !
!     - Determines the total length needed for the regions list by summing the !
!       lengths of the names of all boundary condition regions.                !
!   * Memory Allocation for List:                                              !
!     - Allocates memory for a character array based on the calculated length. !
!   * Forming the List:                                                        !
!     - Constructs the list by concatenating the names of the boundary         !
!       conditions, including their sequential numbers and formatting.         !
!   * Displaying the List:                                                     !
!     - Utilizes the Message function to display the compiled list of boundary !
!       conditions in a framed format for easy reading and verification.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! computational grid
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
