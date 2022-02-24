!==============================================================================!
  module Read_Control_Mod
!------------------------------------------------------------------------------!
!   Module containig functions for reading control file.  They are, to some    !
!   extent, independent from each other, but they all do similar tasks, are    !
!   named simularly, so why to have them scattered in the working directory?   !
!   Besides, being inside a module will enforce argument checking.             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod,    only: Field_Type
  use Var_Mod,      only: Var_Type
  use Vof_Mod,      only: Vof_Type
  use Turb_Mod
  use Swarm_Mod
  use Control_Mod
  use Numerics_Mod
  use Eddies_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------------!
  !   Read control type   !
  !-----------------------!
  type Read_Control_Type

    contains
      procedure          :: Boundary_Conditions
      procedure          :: Linear_Solvers       ! on top of native & PETSc
      procedure, private :: Native_Solvers
      procedure          :: Numerical_Schemes
      procedure, private :: Petsc_Solvers
      procedure          :: Physical_Models
      procedure          :: Physical_Properties

  end type

  type(Read_Control_Type) :: Read_Control

  contains

  ! Member function
  include 'Read_Control_Mod/Boundary_Conditions.f90'
  include 'Read_Control_Mod/Linear_Solvers.f90'
  include 'Read_Control_Mod/Native_Solvers.f90'
  include 'Read_Control_Mod/Numerical_Schemes.f90'
  include 'Read_Control_Mod/Petsc_Solvers.f90'
  include 'Read_Control_Mod/Physical_Models.f90'
  include 'Read_Control_Mod/Physical_Properties.f90'

  end module
