!==============================================================================!
  module Turb_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used for all turbulence modelling paradigms.       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Cpu_Timer_Mod
  use Info_Mod
  use Var_Mod
  use Face_Mod
  use Grid_Mod
  use Field_Mod
  use Solver_Mod
  use Matrix_Mod
  use Control_Mod
  use Numerics_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------------------!
  !                           !
  !      Turbulence type      !
  !   (variables needed for   !
  !   turbulence modelling)   !
  !                           !
  !---------------------------!
  type Turb_Type

    type(Grid_Type),  pointer :: pnt_grid  ! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  ! flow field for which it is defined

    !---------------!
    !   Variables   !
    !---------------!

    ! Single point closures
    type(Var_Type) :: kin
    type(Var_Type) :: eps
    type(Var_Type) :: zeta
    type(Var_Type) :: f22
    type(Var_Type) :: vis

    ! Reynolds stresses and turbulent heat fluxes
    type(Var_Type) :: uu, vv, ww
    type(Var_Type) :: uv, uw, vw
    type(Var_Type) :: ut, vt, wt, t2

    !----------------!
    !   Statistics   !
    !----------------!

    ! Time averaged momentum and energy equations
    real, allocatable :: u_mean(:), v_mean(:), w_mean(:), p_mean(:), t_mean(:)

    ! Time averaged modeled quantities
    ! (Time averages of modeled equations)
    real, allocatable :: kin_mean(:), eps_mean(:), zeta_mean(:), f22_mean(:)
    real, allocatable :: vis_mean(:)

    ! Resolved Reynolds stresses and heat fluxes
    ! (This is what you compute by gathering statistics
    !  with scale-resolving methods like LES and DES)
    real, allocatable :: uu_res(:), vv_res(:), ww_res(:)
    real, allocatable :: uv_res(:), vw_res(:), uw_res(:)
    real, allocatable :: ut_res(:), vt_res(:), wt_res(:), t2_res(:)

    ! Time averaged modelled Reynolds stresses and heat fluxes
    ! (Time averages of modeled equations for
    !  individual Reynolds stress components)
    real, allocatable :: uu_mean(:), vv_mean(:), ww_mean(:)
    real, allocatable :: uv_mean(:), vw_mean(:), uw_mean(:)
    real, allocatable :: ut_mean(:), vt_mean(:), wt_mean(:), t2_mean(:)

    ! Mean passive scalars
    real, allocatable :: scalar_mean(:,:)

    ! Non-dimensional distance
    real, allocatable :: y_plus(:)

    ! Production of turbulent kinetic energy and temperature fluctuations
    real, allocatable :: p_kin(:), p_t2(:)

    ! Buoyancy production for k-eps-zeta-f model
    real, allocatable :: g_buoy(:)

    ! Turbulent lenght and time Scales
    real, allocatable :: l_scale(:)
    real, allocatable :: t_scale(:)

    ! Turbulent viscosity
    real, allocatable :: vis_t(:)

    ! Effective turbulent viscosity
    real, allocatable :: vis_t_eff(:)
    real, allocatable :: vis_t_sgs(:)

    real, allocatable :: tau_wall(:)

    ! Wall viscosity and conductivity (wall function approach)
    real, allocatable :: vis_w(:)
    real, allocatable :: con_w(:)

    ! Scale-resolving simulations
    real, allocatable :: c_dyn(:)
    real, allocatable :: wale_v(:)

    ! Rough walls
    logical :: rough_walls

    ! Hydraulic roughness (constant and variable)
    real              :: z_o
    real, allocatable :: z_o_f(:)

    ! Various cell sizes for Spalart-Allmaras and DES models
    real, allocatable :: h_max(:)
    real, allocatable :: h_min(:)
    real, allocatable :: h_w(:)

    ! Variable for switch between RANS and LES
    real, allocatable :: alpha_l(:)  ! ratio of length scales
    real, allocatable :: alpha_u(:)  ! ratio of velocity scales

    ! For LES you need to know nearest wall cell
    integer, allocatable :: nearest_wall_cell(:)

  end type

  !--------------------------------------------------------!
  !   Parameters and variables defining turbulence model   !
  !--------------------------------------------------------!

  ! Variable holding the turbulence model; its variant and statistics
  integer :: turbulence_model
  integer :: turbulence_model_variant   ! STABILIZED or NONE
  integer :: turbulence_wall_treatment  ! HIGH_RE, LOW_RE, COMPOUND
  logical :: turbulence_statistics
  integer :: turbulent_heat_flux_model
  integer :: hybrid_les_rans_switch

  ! Parameters describing turbulence model choice
  ! (Prime numbers starting from 30000)
  integer, parameter :: NONE                  = 30011
  integer, parameter :: DNS                   = 30013
  integer, parameter :: LES_WALE              = 30029
  integer, parameter :: LES_DYNAMIC           = 30047
  integer, parameter :: LES_SMAGORINSKY       = 30059
  integer, parameter :: K_EPS                 = 30071
  integer, parameter :: K_EPS_ZETA_F          = 30089
  integer, parameter :: DES_SPALART           = 30091
  integer, parameter :: SPALART_ALLMARAS      = 30097
  integer, parameter :: RSM_HANJALIC_JAKIRLIC = 30103
  integer, parameter :: RSM_MANCEAU_HANJALIC  = 30109
  integer, parameter :: HYBRID_LES_RANS       = 30113
  integer, parameter :: HYBRID_LES_PRANDTL    = 30119

  ! Turbulence wall treatment
  integer, parameter :: STABILIZED = 30133

  ! Turbulent heat flux scheme
  integer, parameter :: SGDH = 30137
  integer, parameter :: GGDH = 30139
  integer, parameter :: AFM  = 30161

  ! Switching criteria for hybrid LES/RANS
  integer, parameter :: SWITCH_DISTANCE = 30169
  integer, parameter :: SWITCH_VELOCITY = 30181

  !--------------------------------!
  !   Turbulence model constants   !
  !--------------------------------!

  ! For the k-eps model:
  real :: c_1e, c_2e, c_3e, c_mu, c_mu25, c_mu75, kappa, e_log

  ! For the k-eps-v2f model:
  real :: c_mu_d, c_l, c_t, alpha, c_nu, c_f1, c_f2
  real :: g1, g1_star, g2, g3, g3_star, g4, g5, c_theta
  real :: c_mu_theta, c_mu_theta5, kappa_theta

  ! For the Spalart-Allmaras model:
  real :: c_b1, c_b2, c_w1, c_w2, c_w3, c_v1

  ! For scale-resolving models
  real :: c_smag

  !-----------------------------------!
  !   Auxiliary turbulent variables   !
  !-----------------------------------!

  ! Turbulent Prandtl and Schmidt numbers
  real :: pr_t, sc_t

  contains

  ! Logic of turbulence models
  include 'Turb_Mod/Init.f90'
  include 'Turb_Mod/Main.f90'

  ! The constructor-like
  include 'Turb_Mod/Allocate.f90'

  include 'Turb_Mod/Alias_K_Eps.f90'
  include 'Turb_Mod/Alias_K_Eps_Zeta_F.f90'
  include 'Turb_Mod/Alias_Heat_Fluxes.f90'
  include 'Turb_Mod/Alias_Stresses.f90'
  include 'Turb_Mod/Alias_T2.f90'

  include 'Turb_Mod/Calculate_Deltas.f90'
  include 'Turb_Mod/Calculate_Face_Stress.f90'
  include 'Turb_Mod/Calculate_Face_Vis.f90'
  include 'Turb_Mod/Calculate_Heat_Flux.f90'
  include 'Turb_Mod/Calculate_Mean.f90'
  include 'Turb_Mod/Calculate_Stress.f90'
  include 'Turb_Mod/Substract_Face_Stress.f90'

  ! Functions to set turbulence constants
  include 'Turb_Mod/Const_Hanjalic_Jakirlic.f90'
  include 'Turb_Mod/Const_K_Eps.f90'
  include 'Turb_Mod/Const_K_Eps_Zeta_F.f90'
  include 'Turb_Mod/Const_Manceau_Hanjalic.f90'
  include 'Turb_Mod/Const_Reynolds_Stress.f90'
  include 'Turb_Mod/Const_Spalart_Allmaras.f90'

  ! Computation of various turbulent quantities
  include 'Turb_Mod/Compute_F22.f90'
  include 'Turb_Mod/Compute_Stress.f90'
  include 'Turb_Mod/Compute_Variable.f90'

  ! Different sources
  include 'Turb_Mod/Src_Eps_K_Eps.f90'
  include 'Turb_Mod/Src_Eps_K_Eps_Zeta_F.f90'
  include 'Turb_Mod/Src_F22_K_Eps_Zeta_F.f90'
  include 'Turb_Mod/Src_F22_Rsm_Manceau_Hanjalic.f90'
  include 'Turb_Mod/Src_Kin_K_Eps.f90'
  include 'Turb_Mod/Src_Kin_K_Eps_Zeta_F.f90'
  include 'Turb_Mod/Src_Rsm_Hanjalic_Jakirlic.f90'
  include 'Turb_Mod/Src_Rsm_Manceau_Hanjalic.f90'
  include 'Turb_Mod/Src_T2.f90'
  include 'Turb_Mod/Src_Vis_Spalart_Almaras.f90'
  include 'Turb_Mod/Src_Zeta_K_Eps_Zeta_F.f90'

  ! Computation of turbulence viscosity
  include 'Turb_Mod/Vis_T_Dynamic.f90'
  include 'Turb_Mod/Vis_T_Hybrid_Les_Prandtl.f90'
  include 'Turb_Mod/Vis_T_Hybrid_Les_Rans.f90'
  include 'Turb_Mod/Vis_T_K_Eps.f90'
  include 'Turb_Mod/Vis_T_K_Eps_Zeta_F.f90'
  include 'Turb_Mod/Vis_T_Rsm.f90'
  include 'Turb_Mod/Vis_T_Smagorinsky.f90'
  include 'Turb_Mod/Vis_T_Spalart_Allmaras.f90'
  include 'Turb_Mod/Vis_T_Wale.f90'

  ! Other subroutines ellipitic blending, turbulent Prandtl number
  include 'Turb_Mod/Ebf_Momentum.f90'
  include 'Turb_Mod/Ebf_Scalar.f90'
  include 'Turb_Mod/Prandtl_Number.f90'

  end module
