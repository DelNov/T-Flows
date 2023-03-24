!==============================================================================!
  integer function Numerics_Mod_Pressure_Momentum_Coupling_Code(name)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL) :: name
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
