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

  ! Reynolds stresses
  type(Var_Type) :: uu
  type(Var_Type) :: vv
  type(Var_Type) :: ww
  type(Var_Type) :: uv
  type(Var_Type) :: uw
  type(Var_Type) :: vw
 
  ! Turbulent heat fluxes
  type(Var_Type) :: tt
  type(Var_Type) :: ut
  type(Var_Type) :: vt
  type(Var_Type) :: wt
 
  ! Shear and wall stress are used in a number of turbulence models
  real, allocatable :: shear(:)
  real, allocatable :: vort(:)

  ! Turbulent viscosity
  real, allocatable :: vis_t(:)

  ! Wall viscosity (wall function approuch)
  real, allocatable :: vis_wall(:)
  real, allocatable :: con_wall(:)

  ! Non-dimensional distance
  real, allocatable :: y_plus(:)
 
  ! Friction at the wall and velocity
  real, allocatable :: u_tau(:)
  real, allocatable :: tau_wall(:)

  ! Turbulent Prandtl and Schmidt numbers
  real :: pr_t, sc_t
 
  end module 
