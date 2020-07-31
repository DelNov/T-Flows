!==============================================================================!
  subroutine Control_Mod_Skewness_Correction_Vof(skew_corr, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: skew_corr
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('SKEWNESS_CORRECTION_VOF',   &
                                  'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    skew_corr = .true.

  else if( val .eq. 'NO' ) then
    skew_corr = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for temporal skewure correction: ',   &
                trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
