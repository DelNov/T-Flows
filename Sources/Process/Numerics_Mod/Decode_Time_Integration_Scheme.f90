!==============================================================================!
  subroutine Numerics_Mod_Decode_Time_Integration_Scheme(scheme_name,  &
                                                         scheme_code)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc, Comm_Mod_End
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: scheme_code
  character(len=80) :: scheme_name
!==============================================================================!

  select case(scheme_name)

    case('LINEAR')
      scheme_code = LINEAR
    case('PARABOLIC')
      scheme_code = PARABOLIC

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown time-integration scheme for inertia: ',  &
                 trim(scheme_name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End

  end select

  end subroutine
