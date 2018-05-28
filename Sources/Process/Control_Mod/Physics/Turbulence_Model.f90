!==============================================================================!
  subroutine Control_Mod_Turbulence_Model(verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model from the control file.                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Turbulence_Mod, only: turbulence_model,      &
                            NONE,                  &
                            K_EPS,                 &
                            K_EPS_ZETA_F,          &
                            HYBRID_K_EPS_ZETA_F,   &
                            LES,                   &
                            DNS,                   &
                            DES_SPALART,           &
                            SPALART_ALLMARAS,      &
                            HANJALIC_JAKIRLIC,     &
                            REYNOLDS_STRESS
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('TURBULENCE_MODEL', 'none',  &
                                   val, verbose)
  call To_Upper_Case(val)

  select case(val)

    case('NONE')
      turbulence_model = NONE
    case('K_EPS')                 
      turbulence_model = K_EPS
    case('K_EPS_ZETA_F')          
      turbulence_model = K_EPS_ZETA_F
    case('HYBRID_K_EPS_ZETA_F')   
      turbulence_model = HYBRID_K_EPS_ZETA_F
    case('LES')                   
      turbulence_model = LES
    case('DNS')                   
      turbulence_model = DNS
    case('DES_SPALART')           
      turbulence_model = DES_SPALART
    case('SPALART_ALLMAR AS')     
      turbulence_model = SPALART_ALLMARAS
    case('HANJALIC_JAKIRLIC')     
      turbulence_model = HANJALIC_JAKIRLIC
    case('REYNOLDS_STRESS') 
      turbulence_model = REYNOLDS_STRESS

    case default
      print *, '# Unknown turbulence model :', trim(val)
      print *, '# Exiting!'
      stop 

  end select

  end subroutine
