!==============================================================================!
  integer function Numerics_Mod_Advection_Scheme_Code(name)
!------------------------------------------------------------------------------!
!   Decode the string value on advection scheme from control file              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL) :: name
!==============================================================================!

  select case(name)

    case('UPWIND')
      Numerics_Mod_Advection_Scheme_Code = UPWIND
    case('CENTRAL')
      Numerics_Mod_Advection_Scheme_Code = CENTRAL
    case('LUDS')
      Numerics_Mod_Advection_Scheme_Code = LUDS
    case('QUICK')
      Numerics_Mod_Advection_Scheme_Code = QUICK
    case('SMART')
      Numerics_Mod_Advection_Scheme_Code = SMART
    case('GAMMA')
      Numerics_Mod_Advection_Scheme_Code = GAMMA
    case('MINMOD')
      Numerics_Mod_Advection_Scheme_Code = MINMOD
    case('BLENDED')
      Numerics_Mod_Advection_Scheme_Code = BLENDED
    case('SUPERBEE')
      Numerics_Mod_Advection_Scheme_Code = SUPERBEE
    case('AVL_SMART')
      Numerics_Mod_Advection_Scheme_Code = AVL_SMART
    case('CICSAM')
      Numerics_Mod_Advection_Scheme_Code = CICSAM
    case('STACS')
      Numerics_Mod_Advection_Scheme_Code = STACS

    case default
      call Message % Error(64, 'Unknown advection scheme: '//trim(name)    //  &
                               '. This error is critical, exiting! '       //  &
                               'Check the file: Documents/all_control'     //  &
                               '_keywords to see which advection schemes ' //  &
                               'are currently available in the code.',         &
                               file=__FILE__, line=__LINE__, one_proc=.true.)
  end select

  end function
