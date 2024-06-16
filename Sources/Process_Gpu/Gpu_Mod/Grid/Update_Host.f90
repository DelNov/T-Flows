!==============================================================================!
  subroutine Grid_Update_Host(Gpu, Grid, Turb)
!------------------------------------------------------------------------------!
!>  Copy all the grid variables (so far only the wall distance) you need
!>  for post-processing back to CPU
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type) :: Gpu   !! parent class
  type(Grid_Type) :: Grid  !! field to transfer to device
  type(Turb_Type) :: Turb  !! to check if wall distance should be updated
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(Grid)
    Unused(Turb)
# endif
!==============================================================================!

# if T_FLOWS_GPU == 0
    return
# else
    O_Print '(a)', ' # Copying grid (wall distance) back to the host'
# endif

  if(Turb % model .ne. NO_TURBULENCE_MODEL) then
    call Gpu % Vector_Update_Host(Grid % wall_dist)
  end if

  end subroutine

