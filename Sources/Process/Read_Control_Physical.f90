!==============================================================================!
  subroutine Read_Control_Physical(flow, turb, mult, swarm)
!------------------------------------------------------------------------------!
!   Reads details about physical models from control file.                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod,      only: HUGE_INT
  use Comm_Mod,       only: Comm_Mod_End, this_proc
  use Field_Mod,      only: Field_Type, buoyancy, heat_transfer, t_ref,  &
                            grav_x, grav_y, grav_z
  use Bulk_Mod,       only: Bulk_Type
  use Turb_Mod,       NO_TURBULENCE => NONE
  use Multiphase_Mod, NO_MULTIPHASE => NONE
  use Control_Mod
  use Swarm_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
!----------------------------------[Locals]------------------------------------!
  type(Bulk_Type), pointer :: bulk
  character(len=80)        :: name
  integer                  :: n_times, n_stat, n_stat_p
!==============================================================================!

  ! Take aliases
  bulk => flow % bulk

  !-------------------------------------------!
  !                                           !
  !   Related to heat transfer and bouyancy   !
  !                                           !
  !-------------------------------------------!
  call Control_Mod_Heat_Transfer(heat_transfer, verbose = .true.)
  call Control_Mod_Gravitational_Vector(grav_x,  &
                                        grav_y,  &
                                        grav_z, .true.)
  call Control_Mod_Buoyancy                    (buoyancy,     .true.)
  call Control_Mod_Reference_Temperature       (t_ref,        .true.)
  call Control_Mod_Volume_Expansion_Coefficient(flow % beta,  .true.)
  call Control_Mod_Turbulent_Prandtl_Number    (pr_t)  ! default is (0.9)

  !---------------------------!
  !                           !
  !   Related to turbulence   !
  !                           !
  !---------------------------!
  call Control_Mod_Turbulence_Model(name, .true.)
  select case(name)

    case('NONE')
      turb % model = NO_TURBULENCE
    case('K_EPS')
      turb % model = K_EPS
    case('K_EPS_ZETA_F')
      turb % model = K_EPS_ZETA_F
    case('LES_SMAGORINSKY')
      turb % model = LES_SMAGORINSKY
    case('HYBRID_LES_PRANDTL')
      turb % model = HYBRID_LES_PRANDTL
    case('LES_DYNAMIC')
      turb % model = LES_DYNAMIC
    case('LES_WALE')
      turb % model = LES_WALE
    case('DNS')
      turb % model = DNS
    case('DES_SPALART')
      turb % model = DES_SPALART
    case('SPALART_ALLMARAS')
      turb % model = SPALART_ALLMARAS
    case('RSM_HANJALIC_JAKIRLIC')
      turb % model = RSM_HANJALIC_JAKIRLIC
    case('RSM_MANCEAU_HANJALIC')
      turb % model = RSM_MANCEAU_HANJALIC
    case('HYBRID_LES_RANS')
      turb % model = HYBRID_LES_RANS

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown turbulence model :', trim(name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End

  end select

  ! Does the user want to gather statistics?
  call Control_Mod_Read_Int_Item('NUMBER_OF_TIME_STEPS',               &
                                 0, n_times, .false.)
  call Control_Mod_Read_Int_Item('STARTING_TIME_STEP_FOR_STATISTICS',  &
                                 HUGE_INT, n_stat, .false.)

  !-------------------------------------------------------------------!
  !   For scale-resolving simulations, engage turbulence statistics   !
  !-------------------------------------------------------------------!
  if(turb % model .eq. LES_SMAGORINSKY         .or.  &
     turb % model .eq. LES_DYNAMIC             .or.  &
     turb % model .eq. LES_WALE                .or.  &
     turb % model .eq. DNS                     .or.  &
     turb % model .eq. DES_SPALART             .or.  &
     turb % model .eq. HYBRID_LES_PRANDTL      .or.  &
     turb % model .eq. HYBRID_LES_RANS         .or.  &
     n_times > n_stat) then  ! last line covers unsteady RANS models

    if(this_proc < 2) then
      print *, '# NOTE! Scale resolving simulation used; ' // &
               'turbulence statistics engaged!'
    end if

    turb % statistics = .true.
  end if

  !------------------------------!
  !   Turbulence model variant   !
  !------------------------------!
  call Control_Mod_Turbulence_Model_Variant(name, .true.)
  if     (name .eq. 'NONE') then
    turb % model_variant = NO_TURBULENCE
  else if(name .eq. 'STABILIZED') then
    turb % model_variant = STABILIZED
  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown turbulence model variant: ', trim(name)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End
  end if

  !----------------------------!
  !   Rough or smooth walls?   !
  !----------------------------!
  call Control_Mod_Rough_Walls(turb % rough_walls, .true.)

  !-------------------------------!
  !   Turbulent heat flux model   !
  !-------------------------------!
  call Control_Mod_Turbulent_Heat_Flux_Model(name, .true.)
  select case(name)
    case('SGDH')
      turb % heat_flux_model = SGDH
    case('GGDH')
      turb % heat_flux_model = GGDH
    case('AFM')
      turb % heat_flux_model = AFM
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
  if(turb % model .eq. HYBRID_LES_RANS) then
    call Control_Mod_Hybrid_Les_Rans_Switch(name, .true.)
    select case(name)
      case('SWITCH_DISTANCE')
        turb % hybrid_les_rans_switch = SWITCH_DISTANCE
      case('SWITCH_VELOCITY')
        turb % hybrid_les_rans_switch = SWITCH_VELOCITY
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
  if(turb % model .eq. LES_SMAGORINSKY .or.  &
     turb % model .eq. HYBRID_LES_PRANDTL) then
    call Control_Mod_Smagorinsky_Constant(c_smag, .true.)
  end if

  if(turb % model .eq. K_EPS) then
    call Turb_Mod_Const_K_Eps(turb)
  end if

  if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
    call Turb_Mod_Const_Manceau_Hanjalic(turb)
  end if

  if(turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Turb_Mod_Const_Hanjalic_Jakirlic(turb)
  end if

  if(turb % model .eq. K_EPS_ZETA_F .or.  &
     turb % model .eq. HYBRID_LES_RANS) then
    call Turb_Mod_Const_K_Eps_Zeta_F(turb)
  end if

  if(turb % model .eq. SPALART_ALLMARAS .or.  &
     turb % model .eq. DES_SPALART) then
    call Turb_Mod_Const_Spalart_Allmaras(turb)
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

  !-------------------------------!
  !                               !
  !   Number of passive scalars   !
  !                               !
  !-------------------------------!
  call Control_Mod_Number_Of_Scalars(flow % n_scalars, verbose = .true.)

  !-------------------------------!
  !                               !
  !   Related to Multiphase Flow  !
  !                               !
  !-------------------------------!
  call Control_Mod_Multiphase_Model(name, .true.)

  if(name .eq. 'VOLUME_OF_FLUID' ) then
    mult % model = VOLUME_OF_FLUID
    call Control_Mod_Distance_Function(mult % d_func)
  else if(name .eq. 'LAGRANGIAN_PARTICLES' ) then
    mult % model = LAGRANGIAN_PARTICLES
  else if(name .eq. 'EULER_EULER' ) then
    mult % model = EULER_EULER
  else
    mult % model = NO_MULTIPHASE
  end if

  !-----------------------!
  !                       !
  !   Particle tracking   !
  !                       !
  !-----------------------!
  call Control_Mod_Number_Of_Particles(swarm % n_particles, verbose = .true.)
  call Control_Mod_Swarm_Density      (swarm % density,     verbose = .true.)
  call Control_Mod_Swarm_Diameter     (swarm % diameter,    verbose = .true.)

  call Control_Mod_Swarm_Coefficient_Of_Restitution(swarm % rst,         &
                                                    verbose = .true.)
  call Control_Mod_Number_Of_Swarm_Sub_Steps       (swarm % n_sub_steps, &
                                                    verbose = .true.)

  call Control_Mod_Read_Int_Item('STARTING_TIME_STEP_FOR_SWARM_STATISTICS',  &
                                 HUGE_INT, n_stat_p, .false.)

  ! SGS models for particle
  call Control_Mod_Swarm_Subgrid_Scale_Model(name, verbose = .true.)
  select case(name)
    case('BROWNIAN_FUKAGATA')
         swarm % subgrid_scale_model = BROWNIAN_FUKAGATA
    case('DISCRETE_RANDOM_WALK')
         swarm % subgrid_scale_model = DISCRETE_RANDOM_WALK
  end select

  if(n_times > n_stat_p) then  ! last line covers unsteady RANS models
    if(this_proc < 2) then
      print *, '# NOTE! Lagrangian particle tracking used; ' // &
               'swarm statistics engaged!'                   // &
               'and particle statistics begins at:', n_stat_p
    end if
    swarm % statistics = .true.
  end if
  end subroutine
