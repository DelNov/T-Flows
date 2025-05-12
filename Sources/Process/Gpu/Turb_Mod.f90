#include "../../Shared/Assert.h90"
#include "../../Shared/Browse.h90"
#include "../../Shared/Macros.h90"
#include "../../Shared/Unused.h90"

!==============================================================================!
  module Turb_Mod
!------------------------------------------------------------------------------!
!>  Definition of variables used for all turbulence modelling paradigms.       !
!>
!>  This is s a simplified version from the same module in Process_Cpu.
!>  Hopefully, as more modules are ported to Process_Gpu, this module will
!>  get closer and closer to its origin from Process_Cpu.  In any case, it
!>  shoiuld be our utmost priority to keep this file and its counterpart
!>  in Process_Cpu as close to one another as possible.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Info_Mod
  use Iter_Mod
  use Gpu_Mod
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

    !---------------!
    !   Variables   !
    !---------------!

    ! Single point closures
    type(Var_Type) :: vis

    ! Non-dimensional distance
    real, allocatable :: y_plus(:)

    ! Turbulent viscosity
    real, allocatable :: vis_t(:) ! [kg/(m s)]

    ! Wall viscosity and conductivity (wall function approach)
    real, allocatable :: vis_w(:)   ! [kg/(m s)]
    real, allocatable :: con_w(:)   ! [W/(m K)]
    real, allocatable :: diff_w(:)  ! [m^2/s]

    ! Scale-resolving simulations
    real, allocatable :: wale_v(:)

    ! Are walls treated as rough, and array for hydraulic roughness
    logical           :: rough_walls
    real, allocatable :: z_o(:)

    ! Various cell sizes for Spalart-Allmaras and DES models
    real, allocatable :: h_max(:)
    real, allocatable :: h_min(:)
    real, allocatable :: h_w(:)

    ! Variable holding the turbulence model; its variant and statistics
    integer :: model

    !--------------------------------!
    !   Turbulence model constants   !
    !--------------------------------!

    ! For wall function:
    real :: kappa, e_log, c_mu_theta, c_mu_theta5, kappa_theta

    ! For the Spalart-Allmaras model:
    real :: c_b1, c_b2, c_w1, c_w2, c_w3, c_v1

    ! For AFM turbulent flux model
    real :: afm_eta, afm_psi, c_theta

    ! For scale-resolving models
    real :: c_smag

    contains
      procedure :: Init_Turb
      procedure :: Main_Turb
      procedure :: Create_Turb

      ! Functions to set turbulence constants
      ! They are called from Read_Command_Mod
      procedure :: Const_Spalart_Allmaras
      ! Computation of various turbulent quantities
      procedure, private :: Compute_Variable
      procedure, private :: Form_Variable_Matrix
      procedure, private :: Insert_Variable_Bc

      ! Different sources
      procedure, private :: Src_Vis_Spalart_Allmaras

      ! Computation of turbulence viscosity
      procedure, private :: Vis_T_Subgrid
      procedure, private :: Vis_T_Spalart_Allmaras
      procedure, private :: Vis_T_Wale

      procedure, private :: Beta_Scalar
      procedure, private :: Ebf_Momentum
      procedure, private :: Ebf_Scalar
      procedure          :: Prandtl_Turb

      procedure :: Y_Plus_Rough_Walls
      procedure :: U_Plus_Log_Law
      procedure :: Roughness_Coeff

      procedure :: Les

      ! Procedures to copy turbulence variables to and from the device (GPU)
      procedure :: Copy_Turb_To_Device
      procedure :: Destroy_Turb_On_Device
      procedure :: Update_Turb_On_Host

  end type

  ! Parameters describing turbulence model choice
  ! (Prime numbers starting from 30000)
  integer, parameter :: NO_TURBULENCE_MODEL   = 30011
  integer, parameter :: LES_SMAGORINSKY       = 30029
  integer, parameter :: LES_DYNAMIC           = 30047
  integer, parameter :: LES_WALE              = 30059
  integer, parameter :: LES_TVM               = 30071
  integer, parameter :: DES_SPALART           = 30097
  integer, parameter :: SPALART_ALLMARAS      = 30103

  !-----------------------------------!
  !   Auxiliary turbulent variables   !
  !-----------------------------------!

  ! Turbulent Prandtl and Schmidt numbers
  real :: sc_t

  contains

    ! Logic of turbulence models
#   include "Turb_Mod/Init_Turb.f90"
#   include "Turb_Mod/Main_Turb.f90"

    ! The constructor-like
#   include "Turb_Mod/Create_Turb.f90"


    ! Functions to set turbulence constants
#   include "Turb_Mod/Const_Les.f90"
#   include "Turb_Mod/Const_Spalart_Allmaras.f90"

    ! Computation of various turbulent quantities
#   include "Turb_Mod/Compute_Variable.f90"
#   include "Turb_Mod/Form_Variable_Matrix.f90"
#   include "Turb_Mod/Insert_Variable_Bc.f90"

    ! Different sources
#   include "Turb_Mod/Src_Vis_Spalart_Allmaras.f90"

    ! Computation of turbulence viscosity
#   include "Turb_Mod/Vis_T_Subgrid.f90"
#   include "Turb_Mod/Vis_T_Spalart_Allmaras.f90"
#   include "Turb_Mod/Vis_T_Wale.f90"

    ! Other subroutines ellipitic blending, turbulent Prandtl number
#   include "Turb_Mod/Beta_Scalar.f90"
#   include "Turb_Mod/Ebf_Momentum.f90"
#   include "Turb_Mod/Ebf_Scalar.f90"
#   include "Turb_Mod/Prandtl_Turb.f90"

#   include "Turb_Mod/Y_Plus_Rough_Walls.f90"
#   include "Turb_Mod/U_Plus_Log_Law.f90"
#   include "Turb_Mod/Roughness_Coeff.f90"

#   include "Turb_Mod/Les.f90"

    ! Procedures to copy turbulence variables to and from the device (GPU)
#   include "Turb_Mod/Gpu/Copy_To_Device.f90"
#   include "Turb_Mod/Gpu/Destroy_On_Device.f90"
#   include "Turb_Mod/Gpu/Update_Host.f90"

  end module
