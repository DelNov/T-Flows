!==============================================================================!
  character(SL) function Bnd_Cond_Name_At_Face(Grid, face)
!------------------------------------------------------------------------------!
!>  Provides a shortcut to obtain boundary condition name at boundary face.
!>  It is a bit of the relict of the past, it would be better to browse
!>  through bundary regions first, check the name of the region, and then
!>  browse through all the faces in that region.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid  !! grid under consideration
  integer, intent(in) :: face  !! face number
!==============================================================================!

  Assert(Grid % region % at_face(face) .ge. 0)
  Assert(Grid % region % at_face(face) .le. Grid % per_z_reg)

  if(Grid % region % at_face(face) .eq. 0) then
    Bnd_Cond_Name_At_Face = 'UNDEFINED'
  else
    Bnd_Cond_Name_At_Face = Grid % region % name(Grid % region % at_face(face))
  end if

  end function

