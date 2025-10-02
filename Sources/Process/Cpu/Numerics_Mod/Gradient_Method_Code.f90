!==============================================================================!
  integer function Numerics_Mod_Gradient_Method_Code(name)
!------------------------------------------------------------------------------!
!>  This function in Numerics_Mod module converts a textual representation of
!>  the method for gradient calculation into its corresponding numerical code.
!>  The function ensures that the user input from the control file is accurately
!>  translated to a predefined numerical value that Process can act upon.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(in) :: name  !! textual representation of the
                                     !! method for gradient calculation
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
