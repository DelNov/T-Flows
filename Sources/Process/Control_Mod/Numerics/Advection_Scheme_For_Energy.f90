!==============================================================================!
  subroutine Control_Mod_Advection_Scheme_For_Energy(scheme, verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model from the control file.                            !
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

  call Control_Mod_Read_Char_Item('ADVECTION_SCHEME_FOR_ENERGY', 'upwind',  &
                                   val, verbose)
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
      print *, '# Exiting!'
      stop 

  end select

  end subroutine
