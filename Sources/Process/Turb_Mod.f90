!==============================================================================!
  module Turb_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used for all turbulence modelling paradigms.       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Info_Mod
  use Var_Mod,    only: Var_Type,                    &
                        Var_Mod_Allocate_New_Only,   &
                        Var_Mod_Allocate_Solution,   &
                        Var_Mod_Allocate_Statistics
  use Grid_Mod
  use Grad_Mod
  use Field_Mod,  only: Field_Type,                                            &
                        viscosity, density, buoyancy, conductivity, capacity,  &
                        grav_x,  grav_y,  grav_z,                              &
                        omega_x, omega_y, omega_z,                             &
                        heat_transfer, t_ref
  use Solver_Mod, only: Solver_Type, Bicg, Cg, Cgs
  use Matrix_Mod, only: Matrix_Type
  use Numerics_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------------------!
  !      Turbulence type      !
  !   (variables needed for   !
  !   turbulence modelling)   !
  !---------------------------!
  type Turb_Type

    type(Grid_Type),  pointer :: pnt_grid  ! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  ! flow field for which it is defined

    ! Variables for single-point closures
    ! (These are modelled values)
    type(Var_Type) :: kin
    type(Var_Type) :: eps
    type(Var_Type) :: zeta
    type(Var_Type) :: f22

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

  contains

  ! The constructor-like
  include 'Turb_Mod/Allocate.f90'

  ! Functions to set turbulence constants
  include 'Turb_Mod/Const_Hanjalic_Jakirlic.f90'
  include 'Turb_Mod/Const_K_Eps.f90'
  include 'Turb_Mod/Const_K_Eps_Zeta_F.f90'
  include 'Turb_Mod/Const_Manceau_Hanjalic.f90'
  include 'Turb_Mod/Const_Reynolds_Stress.f90'
  include 'Turb_Mod/Const_Spalart_Allmaras.f90'

  ! Computation of various turbulent quantities
  include 'Turb_Mod/Compute_F22.f90'
  include 'Turb_Mod/Compute_Stresses.f90'
  include 'Turb_Mod/Compute_Turbulent.f90'

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
  include 'Turb_Mod/Vis_T_Hybrid.f90'
  include 'Turb_Mod/Vis_T_K_Eps.f90'
  include 'Turb_Mod/Vis_T_K_Eps_Zeta_F.f90'
  include 'Turb_Mod/Vis_T_Rsm.f90'
  include 'Turb_Mod/Vis_T_Smagorinsky.f90'
  include 'Turb_Mod/Vis_T_Spalart_Allmaras.f90'
  include 'Turb_Mod/Vis_T_Wale.f90'


  end module
