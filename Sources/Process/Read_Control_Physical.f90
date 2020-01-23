!==============================================================================!
  subroutine Read_Control_Physical(flow, mult, swarm, turb)
!------------------------------------------------------------------------------!
!   Reads details about physical models from control file.                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod,      only: HUGE_INT
  use Comm_Mod,       only: Comm_Mod_End, this_proc
  use Field_Mod,      only: Field_Type, buoyancy, heat_transfer, t_ref,  &
                            grav_x, grav_y, grav_z
  use Swarm_Mod,      only: Swarm_Type
  use Bulk_Mod,       only: Bulk_Type
  use Turb_Mod,       NO_TURBULENCE => NONE
  use Multiphase_Mod, NO_MULTIPHASE => NONE
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  type(Turb_Type),       target :: turb
!----------------------------------[Locals]------------------------------------!
  type(Bulk_Type), pointer :: bulk
  character(len=80)        :: name
  integer                  :: n_times, n_stat
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
      turbulence_model = NO_TURBULENCE
    case('K_EPS')
      turbulence_model = K_EPS
    case('K_EPS_ZETA_F')
      turbulence_model = K_EPS_ZETA_F
    case('LES_SMAGORINSKY')
      turbulence_model = LES_SMAGORINSKY
    case('HYBRID_LES_PRANDTL')
      turbulence_model = HYBRID_LES_PRANDTL
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
  if(turbulence_model .eq. LES_SMAGORINSKY         .or.  &
     turbulence_model .eq. LES_DYNAMIC             .or.  &
     turbulence_model .eq. LES_WALE                .or.  &
     turbulence_model .eq. DNS                     .or.  &
     turbulence_model .eq. DES_SPALART             .or.  &
     turbulence_model .eq. HYBRID_LES_PRANDTL      .or.  &
     turbulence_model .eq. HYBRID_LES_RANS         .or.  &
     n_times > n_stat) then  ! last line covers unsteady RANS models

    if(this_proc < 2) then
      print *, '# NOTE! Scale resolving simulation used; ' // &
               'turbulence statistics engaged!'
    end if

    turbulence_statistics = .true.
  end if

  !------------------------------!
  !   Turbulence model variant   !
  !------------------------------!
  call Control_Mod_Turbulence_Model_Variant(name, .true.)
  if     (name .eq. 'NONE') then
    turbulence_model_variant = NO_TURBULENCE
  else if(name .eq. 'STABILIZED') then
    turbulence_model_variant = STABILIZED
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
      turbulent_heat_flux_model = SGDH
    case('GGDH')
      turbulent_heat_flux_model = GGDH
    case('AFM')
      turbulent_heat_flux_model = AFM
    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown turbulent heat flux model :', trim(name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End
  end select

  !-------------------------------------------------------------------------!
  !   Initialization of model constants depending on the turbulence model   !
  !-------------------------------------------------------------------------!
  if(turbulence_model .eq. LES_SMAGORINSKY .or.  &
     turbulence_model .eq. HYBRID_LES_PRANDTL) then
    call Control_Mod_Smagorinsky_Constant(c_smag, .true.)
  end if

  if(turbulence_model .eq. K_EPS) then
    call Turb_Mod_Const_K_Eps(turb)
  end if

  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
    call Turb_Mod_Const_Manceau_Hanjalic(turb)
  end if

  if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Turb_Mod_Const_Hanjalic_Jakirlic(turb)
  end if

  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then
    call Turb_Mod_Const_K_Eps_Zeta_F(turb)
  end if

  if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
     turbulence_model .eq. DES_SPALART) then
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
    multiphase_model = VOLUME_OF_FLUID
    call Control_Mod_Distance_Function(mult % d_func)
  else if(name .eq. 'LAGRANGIAN_PARTICLES' ) then
    multiphase_model = LAGRANGIAN_PARTICLES
  else if(name .eq. 'EULER_EULER' ) then
    multiphase_model = EULER_EULER
  else
    multiphase_model = NO_MULTIPHASE
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

  end subroutine
