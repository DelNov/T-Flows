!==============================================================================!
  module Turbulence_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used for all turbulence modelling paradigms.       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Var_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !--------------------------------------------------------!
  !   Parameters and variables defining turbulence model   !
  !--------------------------------------------------------!

  ! Variable holding the turbulence model; its variant and statistics
  integer :: turbulence_model
  integer :: turbulence_model_variant   ! STABILIZED or NONE
  integer :: turbulence_wall_treatment  ! HIGH_RE, LOW_RE, COMPOUND
  logical :: turbulence_statistics
  integer :: turbulent_heat_flux_model

  ! Parameters describing turbulence model choice
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

  ! Turbulence wall treatment
  integer, parameter :: STABILIZED = 30119

  ! Turbulent heat flux scheme
  integer, parameter :: SGDH = 30121
  integer, parameter :: GGDH = 30123
  integer, parameter :: AFM  = 30125

  logical :: rough_walls

  !---------------------------------!
  !   Turbulence models variables   !
  !---------------------------------!

  ! Variables for single-point closures 
  type(Var_Type), target :: kin
  type(Var_Type), target :: eps
  type(Var_Type), target :: zeta
  type(Var_Type), target :: f22
  type(Var_Type), target :: vis
  type(Var_Type), target :: t2

  ! Reynolds stresses
  type(Var_Type), target :: uu
  type(Var_Type), target :: vv
  type(Var_Type), target :: ww
  type(Var_Type), target :: uv
  type(Var_Type), target :: uw
  type(Var_Type), target :: vw

  ! Turbulent heat fluxes
  type(Var_Type), target :: tt
  type(Var_Type), target :: ut
  type(Var_Type), target :: vt
  type(Var_Type), target :: wt

  !--------------------------------!
  !   Turbulence model constants   !
  !--------------------------------!

  ! For the k-eps model:
  real :: c_1e, c_2e, c_3e, c_mu, c_mu25, c_mu75, kappa, e_log

  ! For the k-eps-v2f model:
  real :: c_mu_d, c_l, c_t, alpha, c_nu, c_f1, c_f2
  real :: g1, g1_star, g2, g3, g3_star, g4, g5, c_theta

  ! For the Spalart-Allmaras model:
  real :: c_b1, c_b2, c_w1, c_w2, c_w3, c_v1

  ! For scale-resolving models
  real              :: c_smag
  real, allocatable :: c_dyn(:)
  real, allocatable :: h_max(:)                                                 
  real, allocatable :: h_min(:)                                                 
  real, allocatable :: h_w(:)                                                   

  !-----------------------------------!
  !   Auxiliary turbulent variables   !
  !-----------------------------------!

  ! Shear and wall stress are used in a number of turbulence models
  ! (These two would make more sense in the Field_Mod)
  real, allocatable :: shear(:)
  real, allocatable :: vort(:)

  ! Turbulent viscosity
  real, allocatable :: vis_t(:)

  real, allocatable :: wale_v(:)

  ! For LES you need to know nearest wall cell
  integer, allocatable :: nearest_wall_cell(:)

  ! Wall viscosity (wall function approuch)
  real, allocatable :: vis_wall(:)
  real, allocatable :: con_wall(:)

  ! Non-dimensional distance
  real, allocatable :: y_plus(:)

  ! Friction at the wall and velocity
  real, allocatable :: u_tau(:)
  real, allocatable :: tau_wall(:)

  ! Effective turbulent viscosity
  real, allocatable :: vis_t_eff(:)
  real, allocatable :: vis_t_sgs(:)

  ! Lenght and Time Scales
  real, allocatable :: l_scale(:)
  real, allocatable :: t_scale(:)

  ! Production of turbulent kinetic energy
  real,allocatable :: p_kin(:), p_t2(:)

  ! Hydraulic roughness (constant and variable)
  real              :: z_o 
  real, allocatable :: z_o_f(:)

  ! Buoyancy production for k-eps-zeta-f model
  ! (bouy_beta is only set to 1 and used as such.  Is it needed?)
  real, allocatable :: g_buoy(:)
  real, allocatable :: buoy_beta(:)
  real, allocatable :: g_kin(:)

  ! Turbulent Prandtl and Schmidt numbers
  real :: pr_t, sc_t

  end module
