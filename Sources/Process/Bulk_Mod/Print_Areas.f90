!==============================================================================!
  subroutine Bulk_Mod_Print_Areas(bulk)
!------------------------------------------------------------------------------!
!>  Prints the computed cross-sectional areas of monitoring planes in the
!>  Bulk_Mod module. This subroutine facilitates the verification of the
!>  placement of the point defining monitoring planes.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Bulk_Type) :: bulk  !! bulk flow properties
!==============================================================================!

  if(First_Proc()) then
    write(*,'(a7,es12.5)') ' # Ax :', bulk % area_x
    write(*,'(a7,es12.5)') ' # Ay :', bulk % area_y
    write(*,'(a7,es12.5)') ' # Az :', bulk % area_z
  end if

  end subroutine
