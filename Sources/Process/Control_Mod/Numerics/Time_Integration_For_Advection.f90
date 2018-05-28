!==============================================================================!
  subroutine Control_Mod_Time_Integration_For_Advection(scheme, verbose)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Numerics_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: scheme
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('TIME_INTEGRATION_FOR_ADVECTION',  &
                                  'fully_implicit',                  &
                                   val,                              &
                                   verbose)
  call To_Upper_Case(val)

  select case(val)

    case('FULLY_IMPLICIT')                 
      scheme = FULLY_IMPLICIT
    case('ADAMS_BASHFORTH')              
      scheme = ADAMS_BASHFORTH
    case('CRANK_NICOLSON')          
      scheme = CRANK_NICOLSON

    case default
      print *, '# Unknown time-integration scheme for advection: ', trim(val)
      print *, '# Exiting!'
      stop 

  end select

  end subroutine
