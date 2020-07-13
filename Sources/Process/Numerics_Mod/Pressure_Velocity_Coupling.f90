!==============================================================================!
  integer function Numerics_Mod_Pressure_Velocity_Coupling(algorithm_name)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: algorithm_name
!==============================================================================!

  select case(algorithm_name)

    case('SIMPLE')
      Numerics_Mod_Pressure_Velocity_Coupling = SIMPLE
    case('PISO')
      Numerics_Mod_Pressure_Velocity_Coupling = PISO

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown coupling algorithm: ',  &
                 trim(algorithm_name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End

  end select

  end function
