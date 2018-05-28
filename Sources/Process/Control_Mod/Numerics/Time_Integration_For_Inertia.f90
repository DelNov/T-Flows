!==============================================================================!
  subroutine Control_Mod_Time_Integration_For_Inertia(scheme, verbose)
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

  call Control_Mod_Read_Char_Item('TIME_INTEGRATION_FOR_INERTIA', 'linear',  &
                                   val, verbose)
  call To_Upper_Case(val)

  select case(val)

    case('LINEAR')                 
      scheme = LINEAR
    case('PARABOLIC')              
      scheme = PARABOLIC

    case default
      print *, '# Unknown time-integration scheme for inertia: ', trim(val)
      print *, '# Exiting!'
      stop 

  end select

  end subroutine
