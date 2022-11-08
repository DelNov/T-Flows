!==============================================================================!
  module Read_Controls_Mod
!------------------------------------------------------------------------------!
!   Module containig functions for reading control file.  They are, to some    !
!   extent, independent from each other, but they all do similar tasks, are    !
!   named simularly, so why to have them scattered in the working directory?   !
!   Besides, being inside a module will enforce argument checking.             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Swarm_Mod
  use Eddies_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------------!
  !   Read control type   !
  !-----------------------!
  type Read_Controls_Type

    contains
      procedure          :: Boundary_Conditions
      procedure          :: Linear_Solvers       ! on top of native & PETSc
      procedure, private :: Native_Solvers
      procedure          :: Numerical_Schemes
      procedure, private :: Petsc_Solvers
      procedure          :: Physical_Models
      procedure          :: Physical_Properties

  end type

  type(Read_Controls_Type) :: Read_Control

  contains

    ! Member function
#   include "Read_Controls_Mod/Boundary_Conditions.f90"
#   include "Read_Controls_Mod/Linear_Solvers.f90"
#   include "Read_Controls_Mod/Native_Solvers.f90"
#   include "Read_Controls_Mod/Numerical_Schemes.f90"
#   include "Read_Controls_Mod/Petsc_Solvers.f90"
#   include "Read_Controls_Mod/Physical_Models.f90"
#   include "Read_Controls_Mod/Physical_Properties.f90"
#   include "Key_Ind.f90"

  end module
