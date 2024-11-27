#include "../../Shared/Assert.h90"
#include "../../Shared/Browse.h90"
#include "../../Shared/Unused.h90"

!==============================================================================!
  module Read_Controls_Mod
!------------------------------------------------------------------------------!
!>  This module consolidates functions for reading and interpreting the control
!>  file used in the Process program of T-Flows.  This module streamlines the
!>  configuration of simulations by organizing parameters related to boundary
!>  conditions, physical models, numerical schemes, linear solvers' details
!>  and physical properties.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Swarm_Mod
  use Eddies_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Interfaces]---------------------------------!
  interface
    include '../../Shared/Key_Ind.h90'
  end interface
!==============================================================================!

  !-----------------------!
  !   Read control type   !
  !-----------------------!
  !> This type encapsulates various procedures, each responsible for parsing
  !> specific sections of the control file.  Specific sections here mean
  !> sections for boundary conditions, numerical schemes, linear solvers'
  !> setting, physical models and physical properties.
  type Read_Controls_Type

    contains
      procedure          :: Boundary_Conditions
      procedure          :: Iterations
      procedure, private :: Native_Solvers
      procedure          :: Numerical_Schemes
      procedure, private :: Petsc_Solvers
      procedure          :: Physical_Models
      procedure          :: Physical_Properties
      procedure          :: Solvers             ! on top of native & PETSc

  end type

  type(Read_Controls_Type) :: Read_Control  !! singleton type object of the
    !! Read_Controls_Type, defined for easier access to its procedures

  contains

    ! Member function
#   include "Read_Controls_Mod/Boundary_Conditions.f90"
#   include "Read_Controls_Mod/Iterations.f90"
#   include "Read_Controls_Mod/Native_Solvers.f90"
#   include "Read_Controls_Mod/Numerical_Schemes.f90"
#   include "Read_Controls_Mod/Petsc_Solvers.f90"
#   include "Read_Controls_Mod/Physical_Models.f90"
#   include "Read_Controls_Mod/Physical_Properties.f90"
#   include "Read_Controls_Mod/Solvers.f90"

  end module
