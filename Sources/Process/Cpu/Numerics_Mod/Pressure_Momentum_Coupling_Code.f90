!==============================================================================!
  integer function Numerics_Mod_Pressure_Momentum_Coupling_Code(name)
!------------------------------------------------------------------------------!
!>  This function in Numerics_Mod module converts a textual representation of
!>  a pressure-momentum scheme into its corresponding numerical code.  The
!>  function ensures that the user input from the control file is accurately
!>  translated to a predefined numerical value that Process can act upon.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(in) :: name  !! textual representation of
                                     !! the pressure-velocity scheme
!==============================================================================!

  select case(name)

    case('SIMPLE')
      Numerics_Mod_Pressure_Momentum_Coupling_Code = SIMPLE

    case('PISO')
      Numerics_Mod_Pressure_Momentum_Coupling_Code = PISO

    case default
      call Message % Error(72,                                                 &
               'Unknown pressure-momentum coupling algorithm: '//trim(name)//  &
               '. \n Exiting!', one_proc=.true.)
  end select

  end function
