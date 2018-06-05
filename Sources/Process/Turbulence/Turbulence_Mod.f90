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
  integer :: turbulence_wall_treatment
  integer :: turbulence_statistics     ! can be YES or NO

  ! Parameters describing turbulence model choice
  integer, parameter :: NONE              = 30011
  integer, parameter :: DNS               = 30013
  integer, parameter :: WALE              = 30029
  integer, parameter :: DYNAMIC           = 30047
  integer, parameter :: SMAGORINSKY       = 30059
  integer, parameter :: K_EPS             = 30071
  integer, parameter :: K_EPS_ZETA_F      = 30089
  integer, parameter :: DES_SPALART       = 30091
  integer, parameter :: SPALART_ALLMARAS  = 30097  
  integer, parameter :: HANJALIC_JAKIRLIC = 30103
  integer, parameter :: REYNOLDS_STRESS   = 30109

  integer :: rough_walls

  ! Turbulence wall treatment 
  integer, parameter :: LOW_RE      = 30509
  integer, parameter :: HIGH_RE     = 30517
  integer, parameter :: COMPOUND    = 30529

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
  real,allocatable :: vis_t(:)

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
