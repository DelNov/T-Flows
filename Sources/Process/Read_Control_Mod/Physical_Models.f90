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
  class(Read_Control_Type) :: Rc
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
!----------------------------------[Locals]------------------------------------!
  type(Bulk_Type), pointer :: bulk
  character(SL)            :: name
  integer                  :: n_times, n_stat, n_stat_p
!==============================================================================!

  ! Take aliases
  bulk => Flow % bulk

  !-------------------------------------------!
  !                                           !
  !   Related to heat transfer and bouyancy   !
  !                                           !
  !-------------------------------------------!
  call Control_Mod_Heat_Transfer(Flow % heat_transfer, verbose = .true.)
  call Control_Mod_Gravitational_Vector(Flow % grav_x,  &
                                        Flow % grav_y,  &
                                        Flow % grav_z, .true.)
  call Control_Mod_Buoyancy(name, .true.)
  select case(name)
    case('NONE')
      Flow % buoyancy = NO_BUOYANCY
    case('DENSITY')
      Flow % buoyancy = DENSITY_DRIVEN
    case('THERMAL')
      Flow % buoyancy = THERMALLY_DRIVEN
    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown buoyancy model :', trim(name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End
      stop
  end select

  call Control_Mod_Reference_Density           (Flow % dens_ref, .true.)
  call Control_Mod_Reference_Temperature       (Flow % t_ref,    .true.)
  call Control_Mod_Volume_Expansion_Coefficient(Flow % beta,     .true.)
  call Control_Mod_Turbulent_Prandtl_Number    (pr_t)  ! default is (0.9)

  !---------------------------!
  !                           !
  !   Related to turbulence   !
  !                           !
  !---------------------------!
  call Control_Mod_Turbulence_Model(name, .true.)
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

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown turbulence model :', trim(name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End
      stop

  end select

  !---------------------------------------------------------!
  !   Turbulence model variant for Reynolds stress models   !
  !---------------------------------------------------------!
  if(Turb % model .eq. RSM_HANJALIC_JAKIRLIC .or.  &
     Turb % model .eq. RSM_MANCEAU_HANJALIC) then
    call Control_Mod_Turbulence_Model_Variant(name, .true.)
    if     (name .eq. 'NONE') then
      Turb % model_variant = NO_TURBULENCE_MODEL
    else if(name .eq. 'STABILIZED') then
      Turb % model_variant = STABILIZED
    else
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown turbulence model variant: ', trim(name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End
    end if
  end if

  !----------------------------!
  !   Rough or smooth walls?   !
  !----------------------------!
  call Control_Mod_Rough_Walls(Turb % rough_walls, .true.)

  ! Does the user want to gather statistics?
  call Control_Mod_Read_Int_Item('NUMBER_OF_TIME_STEPS',               &
                                 0, n_times, .false.)
  call Control_Mod_Read_Int_Item('STARTING_TIME_STEP_FOR_TURB_STATISTICS',  &
                                 HUGE_INT, n_stat, .false.)

  !-------------------------------------------------------------------!
  !   For scale-resolving simulations, engage turbulence statistics   !
  !-------------------------------------------------------------------!
  if((Turb % model .eq. LES_SMAGORINSKY         .or.   &
      Turb % model .eq. LES_DYNAMIC             .or.   &
      Turb % model .eq. LES_WALE                .or.   &
      Turb % model .eq. DNS                     .or.   &
      Turb % model .eq. DES_SPALART             .or.   &
      Turb % model .eq. HYBRID_LES_PRANDTL      .or.   &
      Turb % model .eq. HYBRID_LES_RANS         .or.   &
      Turb % model .eq. K_EPS_ZETA_F)           .and.  &
     n_times > n_stat) then  ! last line covers unsteady RANS models

    if(this_proc < 2) then
      print *, '# NOTE! Scale resolving simulation used; ' // &
               'turbulence statistics engaged!'
    end if

    Turb % statistics = .true.
  end if

  !-------------------------------!
  !   Turbulent heat flux model   !
  !-------------------------------!
  call Control_Mod_Turbulent_Heat_Flux_Model(name, .true.)
  select case(name)
    case('SGDH')
      Turb % heat_flux_model = SGDH
    case('GGDH')
      Turb % heat_flux_model = GGDH
    case('AFM')
      Turb % heat_flux_model = AFM
    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown turbulent heat flux model :', trim(name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End
  end select

  !-------------------------------------------!
  !   Type of switching for hybrid LES/RANS   !
  !-------------------------------------------!
  if(Turb % model .eq. HYBRID_LES_RANS) then
    call Control_Mod_Hybrid_Les_Rans_Switch(name, .true.)
    select case(name)
      case('SWITCH_DISTANCE')
        Turb % hybrid_les_rans_switch = SWITCH_DISTANCE
      case('SWITCH_VELOCITY')
        Turb % hybrid_les_rans_switch = SWITCH_VELOCITY
      case default
        if(this_proc < 2) then
          print *, '# ERROR!  Unknown type of hybrid LES/RANS switch:',  &
                   trim(name)
          print *, '# Exiting!'
        end if
        call Comm_Mod_End
    end select
  end if

  !-------------------------------------------------------------------------!
  !   Initialization of model constants depending on the turbulence model   !
  !-------------------------------------------------------------------------!
  if(Turb % model .eq. LES_SMAGORINSKY .or.  &
     Turb % model .eq. HYBRID_LES_PRANDTL) then
    call Control_Mod_Smagorinsky_Constant(c_smag, .true.)
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

  !------------------------------------!
  !                                    !
  !   Pressure drops and mass fluxes   !
  !                                    !
  !------------------------------------!
  call Control_Mod_Pressure_Drops(bulk % p_drop_x,  &
                                  bulk % p_drop_y,  &
                                  bulk % p_drop_z)
  call Control_Mod_Mass_Flow_Rates(bulk % flux_x_o,  &
                                   bulk % flux_y_o,  &
                                   bulk % flux_z_o)

  !-----------------------!
  !                       !
  !   Number of scalars   !
  !                       !
  !-----------------------!
  call Control_Mod_Number_Of_Scalars(Flow % n_scalars, verbose = .true.)

  !-----------------------------------!
  !                                   !
  !   Related to interface tracking   !
  !                                   !
  !-----------------------------------!
  call Control_Mod_Interface_Tracking(Flow % with_interface, .true.)

  if(Flow % with_interface) then
    call Control_Mod_Track_Front  (Vof % track_front,   .true.)
    call Control_Mod_Track_Surface(Vof % track_surface, .true.)
    call Control_Mod_Mass_Transfer(Flow % mass_transfer)
  end if

  !-----------------------!
  !                       !
  !   Particle tracking   !
  !                       !
  !-----------------------!
  call Control_Mod_Particle_Tracking(Flow % with_particles, .true.)

  if(Flow % with_particles) then

    call Control_Mod_Number_Of_Particles(Swarm % n_particles, verbose = .true.)
    call Control_Mod_Swarm_Density      (Swarm % density,     verbose = .true.)
    call Control_Mod_Swarm_Diameter     (Swarm % diameter,    verbose = .true.)

    call Control_Mod_Swarm_Coefficient_Of_Restitution(Swarm % rst,         &
                                                      verbose = .true.)
    call Control_Mod_Number_Of_Swarm_Sub_Steps       (Swarm % n_sub_steps, &
                                                      verbose = .true.)

    call Control_Mod_Read_Int_Item('STARTING_TIME_STEP_FOR_SWARM_STATISTICS',  &
                                   HUGE_INT, n_stat_p, .false.)

    ! SGS models for particle
    call Control_Mod_Swarm_Subgrid_Scale_Model(name, verbose = .true.)
    select case(name)
      case('BROWNIAN_FUKAGATA')
           Swarm % subgrid_scale_model = BROWNIAN_FUKAGATA
      case('DISCRETE_RANDOM_WALK')
           Swarm % subgrid_scale_model = DISCRETE_RANDOM_WALK
    end select

    if(n_times > n_stat_p) then  ! last line covers unsteady RANS models
      if(this_proc < 2) then
        print *, '# NOTE! Lagrangian particle tracking used; ' // &
                 'swarm statistics engaged!'                   // &
                 'and particle statistics begins at:', n_stat_p
      end if
      Swarm % statistics = .true.
    end if

  end if

  end subroutine
