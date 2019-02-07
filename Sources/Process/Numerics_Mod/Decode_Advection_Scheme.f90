!==============================================================================!
  subroutine Numerics_Mod_Decode_Advection_Scheme(scheme_name, scheme_code)
!------------------------------------------------------------------------------!
!   Decode the string value on advection scheme_code from control file              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc, Comm_Mod_End
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: scheme_name
  integer           :: scheme_code
!==============================================================================!

  select case(scheme_name)

    case('UPWIND')
      scheme_code = UPWIND
    case('CENTRAL')
      scheme_code = CENTRAL
    case('LUDS')
      scheme_code = LUDS
    case('QUICK')
      scheme_code = QUICK
    case('SMART')
      scheme_code = SMART
    case('GAMMA')
      scheme_code = GAMMA
    case('MINMOD')
      scheme_code = MINMOD
    case('BLENDED')
      scheme_code = BLENDED
    case('SUPERBEE')
      scheme_code = SUPERBEE
    case('AVL_SMART')
      scheme_code = AVL_SMART

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Advection scheme_code for energy not found.  Exiting!'
      end if
      call Comm_Mod_End

  end select

  end subroutine
