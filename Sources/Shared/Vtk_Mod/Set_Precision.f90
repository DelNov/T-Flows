!==============================================================================!
  subroutine Vtk_Mod_Set_Precision()
!------------------------------------------------------------------------------!
!>  Sets the strings describing precision of VTK files.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  if(sizeof(1) .eq. SP) then
    intp = '"Int32"'
  else
    intp = '"Int64"'
  end if

  if(sizeof(1.0) .eq. SP) then
    floatp = '"Float32"'
  else
    floatp = '"Float64"'
  end if

  end subroutine
