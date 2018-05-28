!==============================================================================!
  subroutine Cgns_Mod_Get_Arrays_Dimensions(idx, nn_or_nc)
!------------------------------------------------------------------------------!
!   Fetches correct dimensions for arrays in CGNS lib dependent functions      !
!------------------------------------------------------------------------------!
!   Arrays structure in CGNS parallel functions are strictly followings:       !
!                                                                              !
!   Processor:    |        p_1        |               p_2               | ...  !
!   x,y,z:        |      (1 : nn_1)   |       nn_1 + 1 : nn_1 + nn_2    | ...  !
!                                                                              !
!------------------------------------------------------------------------------!
!                                                                              !
!   Cell type:    |      HEXA_8      |     PENTA_6      |       PYRA_5     |...!
!   Connections:  |-p1-|-p2-|...|-pN-|-p1-|-p2-|...|-pN-|-p1-|-p2-|...|-pN-|...!
!                                                                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: idx      !(out)
  integer :: nn_or_nc !(in )
!-----------------------------------[Locals]-----------------------------------!
  integer :: tmp(n_proc)
!------------------------------------------------------------------------------!

  ! single-processor case
  if(this_proc .eq. 0) then
    idx = 1
    return
  end if

  ! multi-processor case
  tmp = 0
  tmp(this_proc) = nn_or_nc

  call Comm_Mod_Global_Sum_Int_Array(tmp, n_proc)

  idx = 1

  if (this_proc > 1) idx = sum(tmp(1:this_proc-1)) + 1

  end subroutine
