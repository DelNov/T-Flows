!==============================================================================!
  subroutine Control_Mod_Time_Integration_Scheme(scheme, verbose)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc, Comm_Mod_End
  use Numerics_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: scheme
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('TIME_INTEGRATION_SCHEME', 'linear',  &
                                   val, verbose)
  call To_Upper_Case(val)

  select case(val)

    case('LINEAR')
      scheme = LINEAR
    case('PARABOLIC')
      scheme = PARABOLIC

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown time-integration scheme for inertia: ',  &
                 trim(val)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End

  end select

  end subroutine
