!==============================================================================!
  subroutine Update_Grid_On_Host(Grid)
!------------------------------------------------------------------------------!
!>  Copy all the grid variables (so far only the wall distance) you need
!>  for post-processing back to CPU
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! field to transfer to device
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Grid)
# endif
!==============================================================================!

# if T_FLOWS_GPU == 0
    return
# else
    O_Print '(a)', ' # Copying grid (wall distance) back to the host'
# endif

  call Gpu % Vector_Update_Host(Grid % wall_dist)

  end subroutine

