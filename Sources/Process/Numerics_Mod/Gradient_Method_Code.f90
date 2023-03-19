!==============================================================================!
  integer function Numerics_Mod_Gradient_Method_Code(name)
!------------------------------------------------------------------------------!
!   Decode the string value on advection scheme from control file              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL) :: name
!==============================================================================!

  select case(name)

    case('LEAST_SQUARES')
      Numerics_Mod_Gradient_Method_Code = LEAST_SQUARES
    case('GAUSS_THEOREM')
      Numerics_Mod_Gradient_Method_Code = GAUSS_THEOREM

    case default
      call Message % Error(80, 'Unknown gradient calculation method: '   //  &
                                trim(name)//'. This error is critical, ' //  &
                               ' exiting! Check the file: '              //  &
                               'Documents/all_control_keywords to see '  //  &
                               'which gradient computation methods are ' //  &
                               'currently available in the code.',           &
                               file=__FILE__, line=__LINE__, one_proc=.true.)
  end select

  end function
