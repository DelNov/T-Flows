!==============================================================================!
  subroutine Control_Mod_Turbulence_Model(verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model from the control file.                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod,      only: HUGE_INT
  use Comm_Mod,       only: this_proc, Comm_Mod_End
  use Turbulence_Mod, only: turbulence_model,          &
                            turbulence_statistics,     &
                            NONE,                      &
                            K_EPS,                     &
                            K_EPS_ZETA_F,              &
                            LES_SMAGORINSKY,           &
                            LES_DYNAMIC,               &
                            LES_WALE,                  &
                            DNS,                       &
                            DES_SPALART,               &
                            SPALART_ALLMARAS,          &
                            RSM_HANJALIC_JAKIRLIC,     &
                            RSM_MANCEAU_HANJALIC,      &
                            HYBRID_LES_RANS
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
  integer           :: n_times, n_stat
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
    case('LES_SMAGORINSKY')
      turbulence_model = LES_SMAGORINSKY
    case('LES_DYNAMIC')
      turbulence_model = LES_DYNAMIC
    case('LES_WALE')
      turbulence_model = LES_WALE
    case('DNS')
      turbulence_model = DNS
    case('DES_SPALART')
      turbulence_model = DES_SPALART
    case('SPALART_ALLMARAS')
      turbulence_model = SPALART_ALLMARAS
    case('RSM_HANJALIC_JAKIRLIC')
      turbulence_model = RSM_HANJALIC_JAKIRLIC
    case('RSM_MANCEAU_HANJALIC')
      turbulence_model = RSM_MANCEAU_HANJALIC
    case('HYBRID_LES_RANS')
      turbulence_model = HYBRID_LES_RANS

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown turbulence model :', trim(val)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End

  end select

  ! Does user want to gather statistics?
  call Control_Mod_Read_Int_Item('NUMBER_OF_TIME_STEPS', 0, n_times, .false.)
  call Control_Mod_Read_Int_Item('STARTING_TIME_STEP_FOR_STATISTICS',  &
                                 HUGE_INT, n_stat, verbose)

  !-------------------------------------------------------------------!
  !   For scale-resolving simulations, engage turbulence statistics   !
  !-------------------------------------------------------------------!
  if(turbulence_model .eq. LES_SMAGORINSKY .or.  &
     turbulence_model .eq. LES_DYNAMIC     .or.  &
     turbulence_model .eq. LES_WALE        .or.  &
     turbulence_model .eq. DNS             .or.  &
     turbulence_model .eq. DES_SPALART     .or.  &
     turbulence_model .eq. HYBRID_LES_RANS .or.  &
     n_times > n_stat) then  ! last line covers unsteady RANS models

    if(this_proc < 2) then
      print *, '# NOTE! Scale resolving simulation used; ' // &
               'turbulence statistics engaged!'
    end if

    turbulence_statistics = .true.
  end if

  end subroutine
