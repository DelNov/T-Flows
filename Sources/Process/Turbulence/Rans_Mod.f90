!==============================================================================!
  module Rans_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used by RANS turbulence models.                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Var_Mod
  use Turbulence_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Turbulence models variables
  type(Var_Type) :: kin
  type(Var_Type) :: eps
  type(Var_Type) :: zeta
  type(Var_Type) :: f22
  type(Var_Type) :: vis

  ! Constants for the k-eps model:
  real :: c_1e, c_2e, c_3e, c_mu, c_mu25, c_mu75, kappa, e_log, Zo
 
  ! Constants for the k-eps-v2f model:
  real :: c_mu_d, c_l, c_t, alpha, Cni, c_f1, c_f2
  real :: g1, g1_star, g2, g3, g3_star, g4, g5 

  ! Constants for the Spalart-Allmaras model:
  real :: c_b1, c_b2, c_w1, c_w2, c_w3, c_v1

  ! Total dissipation in 'HJ' model
  real,allocatable :: eps_tot(:)

  ! Vorticity
  real,allocatable :: vort(:)
  real,allocatable :: vort_mean(:)

  ! Effective turbulent viscosity
  real,allocatable :: vis_t_eff(:)
 
  ! Lenght and Time Scales
  real,allocatable :: l_scale(:)
  real,allocatable :: t_scale(:)   

  ! Production of turbulent kinetic energy
  real,allocatable :: p_kin(:)

  ! Buoyancy production for k-eps-zeta-f model
  ! (bouy_beta is only set to 1 and used as such.  Is it needed?)
  real,allocatable :: g_buoy(:)
  real,allocatable :: buoy_beta(:)
  real,allocatable :: p_buoy(:)
 
  ! Gravity
  real :: grav_x, grav_y, grav_z

  end module 
