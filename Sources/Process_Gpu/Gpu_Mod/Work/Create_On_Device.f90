!==============================================================================!
  subroutine Work_Create_On_Device(Gpu, Work)
!------------------------------------------------------------------------------!
!>  Copy all the working arrays you need in your simulation to GPU.  Clearly
!>  enough, this should be called after the Work_Mod has allocated its arrays.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type) :: Gpu   !! parent class
  type(Work_Type) :: Work  !! to check if wall distance should be coppied
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(Work)
# endif
!==============================================================================!

# if T_FLOWS_GPU == 1
    O_Print '(a)', ' #---------------------------------------'
    O_Print '(a)', ' # Creating working arrays to the device'
    O_Print '(a)', ' #---------------------------------------'
# endif

  O_Print '(a)', ' # Working real arrays'
  do i = 1, Work % Allocated_Real_Cell_Arrays()
    Assert(allocated(Work % r_cell(i) % array))
    call Gpu % Vector_Real_Create_On_Device(Work % r_cell(i) % array)
  end do

  do i = 1, Work % Allocated_Real_Face_Arrays()
    Assert(allocated(Work % r_face(i) % array))
    call Gpu % Vector_Real_Create_On_Device(Work % r_face(i) % array)
  end do

  do i = 1, Work % Allocated_Real_Node_Arrays()
    Assert(allocated(Work % r_node(i) % array))
    call Gpu % Vector_Real_Create_On_Device(Work % r_node(i) % array)
  end do

  O_Print '(a)', ' # Working integer arrays'
  do i = 1, Work % Allocated_Int_Cell_Arrays()
    Assert(allocated(Work % i_cell(i) % array))
    call Gpu % Vector_Int_Create_On_Device(Work % i_cell(i) % array)
  end do

  do i = 1, Work % Allocated_Int_Face_Arrays()
    Assert(allocated(Work % i_face(i) % array))
    call Gpu % Vector_Int_Create_On_Device(Work % i_face(i) % array)
  end do

  do i = 1, Work % Allocated_Int_Node_Arrays()
    Assert(allocated(Work % i_node(i) % array))
    call Gpu % Vector_Int_Create_On_Device(Work % i_node(i) % array)
  end do

  end subroutine

