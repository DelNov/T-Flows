!==============================================================================!
  subroutine Control_Mod_Nodal_Curvature(nodal_curvature, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: nodal_curvature
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('NODAL_CURVATURE', 'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    nodal_curvature = .true.

  else if( val .eq. 'NO' ) then
    nodal_curvature = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for nodal curvature: ',   &
                trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
