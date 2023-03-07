!==============================================================================!
  character(SL) function Bnd_Cond_Name_At_Face(Grid, face)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition name at a face.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
  integer, intent(in) :: face
!==============================================================================!

  Assert(Grid % region % at_face(face) .ge. 0)
  Assert(Grid % region % at_face(face) .le. Grid % per_z_reg)

  if(Grid % region % at_face(face) .eq. 0) then
    Bnd_Cond_Name_At_Face = 'UNDEFINED'
  else
    Bnd_Cond_Name_At_Face = Grid % region % name(Grid % region % at_face(face))
  end if

  end function

