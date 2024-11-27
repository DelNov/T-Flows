#include "../../Shared/Assert.h90"
#include "../../Shared/Browse.h90"
#include "../../Shared/Unused.h90"

!==============================================================================!
  module Turb_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used for all turbulence modelling paradigms.       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Info_Mod
  use Iter_Mod
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

    type(Grid_Type),   pointer :: pnt_grid    ! grid for which it is defined
    type(Field_Type),  pointer :: pnt_flow
    type(Matrix_Type), pointer :: pnt_matrix

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
    real, allocatable :: u_mean(:), v_mean(:), w_mean(:), p_mean(:)
    real, allocatable :: t_mean(:), q_mean(:)

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

    ! Mean passive scalars and scalar turbulent fluxes
    real, allocatable :: scalar_mean(:,:)
    real, allocatable :: uc(:), vc(:), wc(:)

    ! Non-dimensional distance
    real, allocatable :: y_plus(:)

    ! Production of turbulent kinetic energy and temperature fluctuations
    real, allocatable :: p_kin(:), p_t2(:) ! [m^2/s^3], [K^2/s]

    ! Buoyancy production for k-eps-zeta-f model
    real, allocatable :: g_buoy(:)

    ! Turbulent lenght and time Scales
    real, allocatable :: l_scale(:) ! [m]
    real, allocatable :: t_scale(:) ! [s]

    ! Turbulent viscosity
    real, allocatable :: vis_t(:) ! [kg/(m s)]

    ! Tensorial turbulent viscosity (for LES_TVM)
    real, allocatable :: ten_turb_11(:), ten_turb_12(:), ten_turb_13(:)
    real, allocatable :: ten_turb_21(:), ten_turb_22(:), ten_turb_23(:)
    real, allocatable :: ten_turb_31(:), ten_turb_32(:), ten_turb_33(:)

    ! Turbulent stress tensor (for LES_TVM)
    real, allocatable :: tau_11(:), tau_12(:), tau_13(:)
    real, allocatable :: tau_21(:), tau_22(:), tau_23(:)
    real, allocatable :: tau_31(:), tau_32(:), tau_33(:)

    ! Effective turbulent viscosity
    real, allocatable :: vis_t_eff(:) ! [kg/(m s)]
    real, allocatable :: vis_t_sgs(:) ! [kg/(m s)]

    real, allocatable :: tau_wall(:)  ! [kg/(m s^2)]

    ! Wall viscosity and conductivity (wall function approach)
    real, allocatable :: vis_w(:)   ! [kg/(m s)]
    real, allocatable :: con_w(:)   ! [W/(m K)]
    real, allocatable :: diff_w(:)  ! [m^2/s]

    ! Scale-resolving simulations
    real, allocatable :: c_dyn(:)
    real, allocatable :: wale_v(:)

    ! Are walls treated as rough, and array for hydraulic roughness
    logical           :: rough_walls
    real, allocatable :: z_o(:)

    ! Monin-Obukov similarity theory for atmospheric boundary layers
    logical           :: monin_obukov

    ! Various cell sizes for Spalart-Allmaras and DES models
    real, allocatable :: h_max(:)
    real, allocatable :: h_min(:)
    real, allocatable :: h_w(:)

    ! Variable for switch between RANS and LES
    real, allocatable :: alpha_l(:)  ! ratio of length scales
    real, allocatable :: alpha_u(:)  ! ratio of velocity scales

    ! Variable holding the turbulence model; its variant and statistics
    integer :: model
    integer :: model_variant   ! STABILIZED or NONE
    integer :: wall_treatment  ! HIGH_RE, LOW_RE, COMPOUND
    logical :: statistics
    integer :: heat_flux_model
    integer :: scalar_flux_model
    integer :: hybrid_les_rans_switch

    !--------------------------------!
    !   Turbulence model constants   !
    !--------------------------------!

    ! For wall function:
    real :: kappa, e_log, c_mu_theta, c_mu_theta5, kappa_theta

    ! For the k-eps model:
    real :: c_1e, c_2e, c_3e, c_mu, c_mu25, c_mu75

    ! For the k-eps-v2f model:
    real :: c_mu_d, c_l, c_t, alpha, c_nu, c_f1, c_f2

    ! For the Spalart-Allmaras model:
    real :: c_b1, c_b2, c_w1, c_w2, c_w3, c_v1

    ! For HJ and EBM Reynolds Stress Models:
    real :: g1, g1_star, g2, g3, g3_star, g4, g5

    ! For AFM turbulent flux model
    real :: afm_eta, afm_psi, c_theta

    ! For scale-resolving models
    real :: c_smag

    contains
      procedure :: Init_Turb
      procedure :: Main_Turb
      procedure :: Create_Turb

      procedure :: Alias_K_Eps
      procedure :: Alias_K_Eps_Zeta_F
      procedure :: Alias_Heat_Fluxes
      procedure :: Alias_Stresses
      procedure :: Alias_T2
      procedure :: Alias_Vis

      procedure :: Calculate_Deltas
      procedure :: Calculate_Heat_Flux
      procedure :: Calculate_Scalar_Flux
      procedure :: Calculate_Mean
      procedure :: Calculate_Stress
      procedure :: Face_Cond_And_Stress
      procedure :: Face_Diff_And_Stress
      procedure :: Face_Stress
      procedure :: Face_Vis
      procedure :: Substract_Face_Stress

      ! Functions to set turbulence constants
      ! They are called from Read_Command_Mod
      procedure :: Const_Hanjalic_Jakirlic
      procedure :: Const_K_Eps
      procedure :: Const_K_Eps_Zeta_F
      procedure :: Const_Les
      procedure :: Const_Manceau_Hanjalic
      procedure :: Const_Spalart_Allmaras

      ! Computation of various turbulent quantities
      procedure, private :: Compute_F22
      procedure, private :: Compute_Stress
      procedure, private :: Compute_Variable

      ! Different sources
      procedure, private :: Src_Eps_K_Eps
      procedure, private :: Src_Eps_K_Eps_Zeta_F
      procedure, private :: Src_F22_K_Eps_Zeta_F
      procedure, private :: Src_F22_Rsm_Manceau_Hanjalic
      procedure, private :: Src_Kin_K_Eps
      procedure, private :: Src_Kin_K_Eps_Zeta_F
      procedure, private :: Src_Rsm_Hanjalic_Jakirlic
      procedure, private :: Src_Rsm_Manceau_Hanjalic
      procedure, private :: Src_T2
      procedure, private :: Src_Vis_Spalart_Allmaras
      procedure, private :: Src_Zeta_K_Eps_Zeta_F

      ! Computation of turbulence viscosity
      procedure          :: Vis_T_Dynamic             ! also called from Swarm
      procedure, private :: Vis_T_Hybrid_Les_Prandtl
      procedure, private :: Vis_T_Hybrid_Les_Rans
      procedure, private :: Vis_T_K_Eps
      procedure, private :: Vis_T_K_Eps_Zeta_F
      procedure, private :: Vis_T_Rsm
      procedure, private :: Vis_T_Subgrid
      procedure, private :: Vis_T_Spalart_Allmaras
      procedure, private :: Vis_T_Wale
      procedure, private :: Vis_T_Tensorial

      procedure, private :: Beta_Scalar
      procedure, private :: Ebf_Momentum
      procedure, private :: Ebf_Scalar
      procedure          :: Prandtl_Turb

      procedure :: Y_Plus_Rough_Walls
      procedure :: Tau_Wall_Log_Law
      procedure :: U_Plus_Log_Law
      procedure :: Time_And_Length_Scale
      procedure :: Roughness_Coeff
      procedure :: Monin_Obukov_Momentum
      procedure :: Monin_Obukov_Thermal

      procedure :: Les
      procedure :: Rsm

  end type

  ! Parameters describing turbulence model choice
  ! (Prime numbers starting from 30000)
  integer, parameter :: NO_TURBULENCE_MODEL   = 30011
  integer, parameter :: DNS                   = 30013
  integer, parameter :: LES_SMAGORINSKY       = 30029
  integer, parameter :: LES_DYNAMIC           = 30047
  integer, parameter :: LES_WALE              = 30059
  integer, parameter :: LES_TVM               = 30071
  integer, parameter :: K_EPS                 = 30089
  integer, parameter :: K_EPS_ZETA_F          = 30091
  integer, parameter :: DES_SPALART           = 30097
  integer, parameter :: SPALART_ALLMARAS      = 30103
  integer, parameter :: RSM_HANJALIC_JAKIRLIC = 30109
  integer, parameter :: RSM_MANCEAU_HANJALIC  = 30113
  integer, parameter :: HYBRID_LES_RANS       = 30119
  integer, parameter :: HYBRID_LES_PRANDTL    = 30133

  ! Turbulence wall treatment
  integer, parameter :: STABILIZED = 30137

  ! Turbulent heat flux scheme
  integer, parameter :: SGDH = 30139
  integer, parameter :: GGDH = 30161
  integer, parameter :: AFM  = 30169

  ! Switching criteria for hybrid LES/RANS
  integer, parameter :: SWITCH_DISTANCE = 30181
  integer, parameter :: SWITCH_VELOCITY = 30187

  !-----------------------------------!
  !   Auxiliary turbulent variables   !
  !-----------------------------------!

  ! Turbulent Prandtl and Schmidt numbers
  real :: pr_t, sc_t

  contains

    ! Logic of turbulence models
#   include "Turb_Mod/Init_Turb.f90"
#   include "Turb_Mod/Main_Turb.f90"

    ! The constructor-like
#   include "Turb_Mod/Create_Turb.f90"

#   include "Turb_Mod/Alias_K_Eps.f90"
#   include "Turb_Mod/Alias_K_Eps_Zeta_F.f90"
#   include "Turb_Mod/Alias_Heat_Fluxes.f90"
#   include "Turb_Mod/Alias_Stresses.f90"
#   include "Turb_Mod/Alias_T2.f90"
#   include "Turb_Mod/Alias_Vis.f90"

#   include "Turb_Mod/Calculate_Deltas.f90"
#   include "Turb_Mod/Calculate_Heat_Flux.f90"
#   include "Turb_Mod/Calculate_Scalar_Flux.f90"
#   include "Turb_Mod/Calculate_Mean.f90"
#   include "Turb_Mod/Calculate_Stress.f90"
#   include "Turb_Mod/Face_Cond_And_Stress.f90"
#   include "Turb_Mod/Face_Diff_and_Stress.f90"
#   include "Turb_Mod/Face_Stress.f90"
#   include "Turb_Mod/Face_Vis.f90"
#   include "Turb_Mod/Substract_Face_Stress.f90"

    ! Functions to set turbulence constants
#   include "Turb_Mod/Const_Hanjalic_Jakirlic.f90"
#   include "Turb_Mod/Const_K_Eps.f90"
#   include "Turb_Mod/Const_K_Eps_Zeta_F.f90"
#   include "Turb_Mod/Const_Les.f90"
#   include "Turb_Mod/Const_Manceau_Hanjalic.f90"
#   include "Turb_Mod/Const_Spalart_Allmaras.f90"

    ! Computation of various turbulent quantities
#   include "Turb_Mod/Compute_F22.f90"
#   include "Turb_Mod/Compute_Stress.f90"
#   include "Turb_Mod/Compute_Variable.f90"

    ! Different sources
#   include "Turb_Mod/Src_Eps_K_Eps.f90"
#   include "Turb_Mod/Src_Eps_K_Eps_Zeta_F.f90"
#   include "Turb_Mod/Src_F22_K_Eps_Zeta_F.f90"
#   include "Turb_Mod/Src_F22_Rsm_Manceau_Hanjalic.f90"
#   include "Turb_Mod/Src_Kin_K_Eps.f90"
#   include "Turb_Mod/Src_Kin_K_Eps_Zeta_F.f90"
#   include "Turb_Mod/Src_Rsm_Hanjalic_Jakirlic.f90"
#   include "Turb_Mod/Src_Rsm_Manceau_Hanjalic.f90"
#   include "Turb_Mod/Src_T2.f90"
#   include "Turb_Mod/Src_Vis_Spalart_Allmaras.f90"
#   include "Turb_Mod/Src_Zeta_K_Eps_Zeta_F.f90"

    ! Computation of turbulence viscosity
#   include "Turb_Mod/Vis_T_Dynamic.f90"
#   include "Turb_Mod/Vis_T_Hybrid_Les_Prandtl.f90"
#   include "Turb_Mod/Vis_T_Hybrid_Les_Rans.f90"
#   include "Turb_Mod/Vis_T_K_Eps.f90"
#   include "Turb_Mod/Vis_T_K_Eps_Zeta_F.f90"
#   include "Turb_Mod/Vis_T_Rsm.f90"
#   include "Turb_Mod/Vis_T_Subgrid.f90"
#   include "Turb_Mod/Vis_T_Spalart_Allmaras.f90"
#   include "Turb_Mod/Vis_T_Wale.f90"
#   include "Turb_Mod/Vis_T_Tensorial.f90"

    ! Other subroutines ellipitic blending, turbulent Prandtl number
#   include "Turb_Mod/Beta_Scalar.f90"
#   include "Turb_Mod/Ebf_Momentum.f90"
#   include "Turb_Mod/Ebf_Scalar.f90"
#   include "Turb_Mod/Prandtl_Turb.f90"

#   include "Turb_Mod/Y_Plus_Rough_Walls.f90"
#   include "Turb_Mod/Tau_Wall_Log_Law.f90"
#   include "Turb_Mod/U_Plus_Log_Law.f90"
#   include "Turb_Mod/Time_And_Length_Scale.f90"
#   include "Turb_Mod/Roughness_Coeff.f90"
#   include "Turb_Mod/Monin_Obukov_Momentum.f90"
#   include "Turb_Mod/Monin_Obukov_Thermal.f90"

#   include "Turb_Mod/Les.f90"
#   include "Turb_Mod/Rsm.f90"

  end module
