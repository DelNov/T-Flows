!==============================================================================!
  subroutine Control_Mod_Turbulence_Model(verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model from the control file.                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Turbulence_Mod, only: turbulence_model,          &
                            turbulence_statistics,     &
                            NONE,                      &
                            K_EPS,                     &
                            K_EPS_ZETA_F,              &
                            SMAGORINSKY,               &
                            DYNAMIC,                   &
                            WALE,                      &
                            DNS,                       &
                            DES_SPALART,               &
                            SPALART_ALLMARAS,          &
                            HANJALIC_JAKIRLIC,         &
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

  !-----------------------------!
  !   Select turbulence model   !
  !-----------------------------!
  select case(val)

    case('NONE')
      turbulence_model = NONE
    case('K_EPS')
      turbulence_model = K_EPS
    case('K_EPS_ZETA_F')
      turbulence_model = K_EPS_ZETA_F
    case('SMAGORINSKY')
      turbulence_model = SMAGORINSKY
    case('DYNAMIC')
      turbulence_model = DYNAMIC
    case('WALE')
      turbulence_model = WALE
    case('DNS')
      turbulence_model = DNS
    case('DES_SPALART')
      turbulence_model = DES_SPALART
    case('SPALART_ALLMARAS')
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

  !-------------------------------------------------------------------!
  !   For scale-resolving simulations, engage turbulence statistics   !
  !-------------------------------------------------------------------!
  if(turbulence_model .eq. SMAGORINSKY .or.  &
     turbulence_model .eq. DYNAMIC     .or.  &
     turbulence_model .eq. WALE        .or.  &
     turbulence_model .eq. DNS         .or.  &
     turbulence_model .eq. DES_SPALART) then

    print *, '# Scale resolving simulation used; ' // &
             'turbulence statistics engaged!'

    turbulence_statistics = .true.
  end if

  end subroutine
