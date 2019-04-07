!==============================================================================!
  module Var_Mod
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  !--------------!
  !   Var type   !
  !--------------!
  type Var_Type

    type(Grid_Type), pointer :: pnt_grid  ! grid for which it is defined

    character(len=4)     :: name                  ! variable name, always
                                                  ! uppercase and very short
    character(len=4)     :: flux_name             ! variable flux name, always
                                                  ! uppercase and very short
    real, allocatable    :: n(:)                  ! new value
    real, allocatable    :: o(:), oo(:)           ! old and older then old
    real, allocatable    :: a(:)                  ! advection fluxes
    real, allocatable    :: c(:)                  ! cross-difusion fluxes
    real, allocatable    :: x(:), y(:), z(:)      ! gradient components
    real, allocatable    :: q(:)                  ! flux of a variable
    real                 :: sigma                 ! sigma
    real                 :: res                   ! residual after lin. solver
    real                 :: units(5)              ! mass, length, time,
                                                  ! temperature, angle
    ! Boundary cell type (important for scalars, since they
    ! can have different boundary conditions at the walls)
    integer, allocatable :: bnd_cell_type(:)

    ! Parameters for numerical solution of the variable
    character(len=80)    :: precond
    integer              :: adv_scheme  ! advection scheme
    real                 :: blend       ! blending (1.0 central; 0.0 upwind)
    integer              :: td_scheme   ! time-disretization
    real                 :: tol         ! linear solver tolerance
    real                 :: urf         ! under-relaxation factor
    integer              :: niter       ! number of iterations
    real, allocatable    :: max(:)      ! max and min around a face ...
    real, allocatable    :: min(:)      ! important for advection schemes
  end type

  contains

  include 'Var_Mod/Allocate_New_Only.f90'
  include 'Var_Mod/Allocate_Solution.f90'
  include 'Var_Mod/Bnd_Cell_Type.f90'

  end module
