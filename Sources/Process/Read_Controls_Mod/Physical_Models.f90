!==============================================================================!
  subroutine Physical_Models(Rc, Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   Reads details about physical models from control file.                     !
!                                                                              !
!   Good practice: default values, outlined in Documents/all_control_keywords, !
!   should be defined only in Control_Mod, not be scattered around the code.   !
!   In other words, Control_Mod changes less frequently than other parts of    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc
  type(Field_Type), target              :: Flow
  type(Turb_Type),  target              :: Turb
  type(Vof_Type),   target              :: Vof
  type(Swarm_Type), target              :: Swarm
!----------------------------------[Locals]------------------------------------!
  type(Bulk_Type), pointer :: bulk
  character(SL)            :: name
  integer                  :: n_times, n_stat, n_stat_p
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
!==============================================================================!

  ! Give some sign
  if(First_Proc())  &
    print '(a)', ' # Reading about physical models'

  ! Take aliases
  bulk => Flow % bulk

  !--------------------------!
  !                          !
  !   Number of time steps   !
  !                          !
  !--------------------------!
  call Control % Read_Int_Item('NUMBER_OF_TIME_STEPS', 0, n_times, .false.)

  !-------------------------------------------!
  !                                           !
  !   Related to heat transfer and bouyancy   !
  !                                           !
  !-------------------------------------------!
  call Control % Heat_Transfer(Flow % heat_transfer, verbose = .true.)
  call Control % Gravitational_Vector(Flow % grav_x,  &
                                      Flow % grav_y,  &
                                      Flow % grav_z, .true.)
  call Control % Buoyancy(name, .true.)
  select case(name)
    case('NONE')
      Flow % buoyancy = NO_BUOYANCY
    case('DENSITY')
      Flow % buoyancy = DENSITY_DRIVEN
    case('THERMAL')
      Flow % buoyancy = THERMALLY_DRIVEN
    case default
      call Message % Error(60,                                       &
                           'Unknown buoyancy model: '//trim(name)//  &
                           '.  \n Exiting!')
  end select

  call Control % Reference_Density           (Flow % dens_ref, .true.)
  call Control % Reference_Temperature       (Flow % t_ref,    .true.)
  call Control % Volume_Expansion_Coefficient(Flow % beta,     .true.)
  call Control % Turbulent_Prandtl_Number    (pr_t)  ! default is (0.9)
  call Control % Turbulent_Schmidt_Number    (sc_t)  ! default is (0.9)
  call Control % Extrapolate_Temperature_Exp (Flow % exp_temp_wall, .true.)

  !---------------------------!
  !                           !
  !   Related to turbulence   !
  !                           !
  !---------------------------!
  call Control % Turbulence_Model(name, .true.)
  select case(name)

    case('NONE')
      Turb % model = NO_TURBULENCE_MODEL
    case('K_EPS')
      Turb % model = K_EPS
    case('K_EPS_ZETA_F')
      Turb % model = K_EPS_ZETA_F
    case('LES_SMAGORINSKY')
      Turb % model = LES_SMAGORINSKY
    case('HYBRID_LES_PRANDTL')
      Turb % model = HYBRID_LES_PRANDTL
    case('LES_DYNAMIC')
      Turb % model = LES_DYNAMIC
    case('LES_WALE')
      Turb % model = LES_WALE
    case('DNS')
      Turb % model = DNS
    case('DES_SPALART')
      Turb % model = DES_SPALART
    case('SPALART_ALLMARAS')
      Turb % model = SPALART_ALLMARAS
    case('RSM_HANJALIC_JAKIRLIC')
      Turb % model = RSM_HANJALIC_JAKIRLIC
    case('RSM_MANCEAU_HANJALIC')
      Turb % model = RSM_MANCEAU_HANJALIC
    case('HYBRID_LES_RANS')
      Turb % model = HYBRID_LES_RANS
    case('LES_TVM')
      Turb % model = LES_TVM

    case default
      call Message % Error(60,                                         &
                           'Unknown turbulence model: '//trim(name)//  &
                           '.  \n Exiting!')
  end select

  !---------------------------------------------------------!
  !   Turbulence model variant for Reynolds stress models   !
  !---------------------------------------------------------!
  if(Turb % model .eq. RSM_HANJALIC_JAKIRLIC .or.  &
     Turb % model .eq. RSM_MANCEAU_HANJALIC) then
    call Control % Turbulence_Model_Variant(name, .true.)
    if     (name .eq. 'NONE') then
      Turb % model_variant = NO_TURBULENCE_MODEL
    else if(name .eq. 'STABILIZED') then
      Turb % model_variant = STABILIZED
    else
      call Message % Error(72,                                         &
                   'Unknown turbulence model variant: '//trim(name)//  &
                   '.  \n Exiting!')
    end if
  end if

  !----------------------------!
  !   Rough or smooth walls?   !
  !----------------------------!
  call Control % Rough_Walls(Turb % rough_walls, .true.)

  !----------------------------!
  !   Monin-Obukov for ABL?    !
  !----------------------------!
  call Control % Monin_Obukov(Turb % monin_obukov, .true.)

  ! Does the user want to gather statistics?
  call Control % Read_Int_Item('STARTING_TIME_STEP_FOR_TURB_STATISTICS',  &
                               HUGE_INT, n_stat, .false.)

  !-------------------------------------------------------------------!
  !   For scale-resolving simulations, engage turbulence statistics   !
  !-------------------------------------------------------------------!
  if((Turb % model .eq. LES_SMAGORINSKY         .or.   &
      Turb % model .eq. LES_DYNAMIC             .or.   &
      Turb % model .eq. LES_WALE                .or.   &
      Turb % model .eq. LES_TVM                 .or.   &
      Turb % model .eq. DNS                     .or.   &
      Turb % model .eq. DES_SPALART             .or.   &
      Turb % model .eq. HYBRID_LES_PRANDTL      .or.   &
      Turb % model .eq. HYBRID_LES_RANS         .or.   &
      Turb % model .eq. K_EPS                   .or.   &
      Turb % model .eq. K_EPS_ZETA_F)           .and.  &
     n_times > n_stat) then  ! last line covers unsteady RANS models

    if(First_Proc()) then
      print *, '# NOTE! Scale resolving simulation used; ' // &
               'turbulence statistics engaged!'
    end if

    Turb % statistics = .true.
  end if

  !-------------------------------!
  !   Turbulent heat flux model   !
  !-------------------------------!
  call Control % Turbulent_Heat_Flux_Model(name, .true.)
  select case(name)
    case('SGDH')
      Turb % heat_flux_model = SGDH
    case('GGDH')
      Turb % heat_flux_model = GGDH
    case('AFM')
      Turb % heat_flux_model = AFM
    case default
      call Message % Error(64,                                          &
                   'Unknown turbulent heat flux model: '//trim(name)//  &
                   '.  \n Exiting!')
  end select

  !-------------------------------------------!
  !   Type of switching for hybrid LES/RANS   !
  !-------------------------------------------!
  if(Turb % model .eq. HYBRID_LES_RANS) then
    call Control % Hybrid_Les_Rans_Switch(name, .true.)
    select case(name)
      case('SWITCH_DISTANCE')
        Turb % hybrid_les_rans_switch = SWITCH_DISTANCE
      case('SWITCH_VELOCITY')
        Turb % hybrid_les_rans_switch = SWITCH_VELOCITY
      case default
        call Message % Error(72,                                            &
                  'Unknown type of hybrid LES/RANS switch: '//trim(name)//  &
                  '.  \n Exiting!')
    end select
  end if

  !-------------------------------------------------------------------------!
  !   Initialization of model constants depending on the turbulence model   !
  !-------------------------------------------------------------------------!
  if(Turb % model .eq. LES_SMAGORINSKY .or.  &
     Turb % model .eq. HYBRID_LES_PRANDTL) then
    call Control % Smagorinsky_Constant(c_smag, .true.)
  end if

  if(Turb % model .eq. K_EPS) then
    call Turb % Const_K_Eps()
  end if

  if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
    call Turb % Const_Manceau_Hanjalic()
  end if

  if(Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Turb % Const_Hanjalic_Jakirlic()
  end if

  if(Turb % model .eq. K_EPS_ZETA_F .or.  &
     Turb % model .eq. HYBRID_LES_RANS) then
    call Turb % Const_K_Eps_Zeta_F()
  end if

  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then
    call Turb % Const_Spalart_Allmaras()
  end if

  if(Turb % model .eq. LES_DYNAMIC .or.  &
     Turb % model .eq. LES_SMAGORINSKY) then
    call Turb % Const_Les()
  end if

  !------------------------------------!
  !                                    !
  !   Pressure drops and mass fluxes   !
  !                                    !
  !------------------------------------!
  call Control % Pressure_Drops(bulk % p_drop_x,  &
                                bulk % p_drop_y,  &
                                bulk % p_drop_z)
  call Control % Volume_Flow_Rates(bulk % flux_x_o,  &
                                   bulk % flux_y_o,  &
                                   bulk % flux_z_o)

  !-----------------------!
  !                       !
  !   Number of scalars   !
  !                       !
  !-----------------------!
  call Control % Number_Of_Scalars(Flow % n_scalars, verbose = .true.)

  !-----------------------------------!
  !                                   !
  !   Related to interface tracking   !
  !                                   !
  !-----------------------------------!
  call Control % Interface_Tracking(Flow % with_interface, .true.)

  if(Flow % with_interface) then
    call Control % Track_Front  (Vof % track_front,   .true.)
    call Control % Track_Surface(Vof % track_surface, .true.)
    call Control % Mass_Transfer(Flow % mass_transfer)
  end if

  !-----------------------!
  !                       !
  !   Particle tracking   !
  !                       !
  !-----------------------!
  call Control % Particle_Tracking(Flow % with_particles, .true.)

  if(Flow % with_particles) then

    call Control % Max_Particles (Swarm % max_particles, verbose = .true.)
    call Control % Swarm_Density (Swarm % density,       verbose = .true.)
    call Control % Swarm_Diameter(Swarm % diameter,      verbose = .true.)

    call Control % Swarm_Coefficient_Of_Restitution(Swarm % rst,         &
                                                    verbose = .true.)
    call Control % Number_Of_Swarm_Sub_Steps       (Swarm % n_sub_steps, &
                                                    verbose = .true.)

    call Control % Read_Int_Item('STARTING_TIME_STEP_FOR_SWARM_STATISTICS',  &
                                 HUGE_INT, n_stat_p, .false.)

    ! SGS models for particle
    call Control % Swarm_Subgrid_Scale_Model(name, verbose = .true.)
    select case(name)
      case('BROWNIAN_FUKAGATA')
           Swarm % subgrid_scale_model = BROWNIAN_FUKAGATA
      case('DISCRETE_RANDOM_WALK')
           Swarm % subgrid_scale_model = DISCRETE_RANDOM_WALK
    end select

    if(n_times > n_stat_p) then  ! last line covers unsteady RANS models
      if(First_Proc()) then
        print *, '# NOTE! Lagrangian particle tracking used; ' // &
                 'swarm statistics engaged!'                   // &
                 'and particle statistics begins at:', n_stat_p
      end if
      Swarm % statistics = .true.
    end if

  end if

  end subroutine
