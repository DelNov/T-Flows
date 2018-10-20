!==============================================================================!
  subroutine Control_Mod_Pressure_Momentum_Coupling(verbose)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc, Comm_Mod_End
  use Numerics_Mod, only: pressure_momentum_coupling, SIMPLE, PROJECTION
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  ! Set a default value
  pressure_momentum_coupling = SIMPLE

  call Control_Mod_Read_Char_Item('PRESSURE_MOMENTUM_COUPLING', 'simple',  &
                                   val, verbose)
  call To_Upper_Case(val)

  !----------------------------!
  !   Select coupling method   !
  !----------------------------!
  select case(val)

    case('SIMPLE')
      pressure_momentum_coupling = SIMPLE
    case('PROJECTION')
      pressure_momentum_coupling = PROJECTION

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown pressure-momentum coupling: ', trim(val)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End

  end select

  end subroutine
