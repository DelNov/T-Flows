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
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown pressure-momentum coupling algorithm: ',  &
                 trim(name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End
      stop

  end select

  end function
