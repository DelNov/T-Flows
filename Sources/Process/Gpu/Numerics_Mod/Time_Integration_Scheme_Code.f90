!==============================================================================!
  integer function Numerics_Mod_Time_Integration_Scheme_Code(name)
!------------------------------------------------------------------------------!
!>  This function in Numerics_Mod module converts a textual representation of
!>  a time integration scheme into its corresponding numerical code.  The
!>  function ensures that the user input from the control file is accurately
!>  translated to a predefined numerical value that Process can act upon.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(in) :: name  !! textual representation of
                                     !! the time-integration scheme
!==============================================================================!

  select case(name)

    case('LINEAR')
      Numerics_Mod_Time_Integration_Scheme_Code = LINEAR
    case('PARABOLIC')
      Numerics_Mod_Time_Integration_Scheme_Code = PARABOLIC

    case default
      call Message % Error(60,                                    &
               'Unknown time-integration scheme: '//trim(name)//  &
               '. \n Exiting!', one_proc=.true.)
  end select

  end function
