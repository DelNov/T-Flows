!==============================================================================!
  subroutine Control_Mod_Advection_Scheme_For_User_Scalars(scheme, verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model from the control file.                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod,     only: this_proc, Comm_Mod_End
  use Numerics_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: scheme
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('ADVECTION_SCHEME_FOR_USER_SCALARS',  &
                                  'upwind', val, verbose)
  call To_Upper_Case(val)

  select case(val)

    case('UPWIND')
      scheme = UPWIND
    case('CENTRAL')
      scheme = CENTRAL
    case('LUDS')
      scheme = LUDS
    case('QUICK')
      scheme = QUICK
    case('SMART')
      scheme = SMART
    case('GAMMA')
      scheme = GAMMA
    case('MINMOD')
      scheme = MINMOD
    case('BLENDED')
      scheme = BLENDED
    case('SUPERBEE')
      scheme = SUPERBEE
    case('AVL_SMART')
      scheme = AVL_SMART

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Advection scheme for user scalar not found.' //  &
                 '  Exiting!'
      end if
      call Comm_Mod_End

  end select

  end subroutine
