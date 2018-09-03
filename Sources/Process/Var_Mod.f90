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
    real, allocatable    :: a(:), a_o(:), a_oo(:) ! advection fluxes
    real, allocatable    :: d_o(:), d_oo(:)       ! difussion fluxes
    real, allocatable    :: c(:), c_o(:), c_oo(:) ! cross-difusion
    real, allocatable    :: mean(:)               ! time average
    real, allocatable    :: x(:), y(:), z(:)      ! gradient components
    real, allocatable    :: q(:)                  ! flux of a variable
    real                 :: sigma                 ! sigma
    real                 :: res                   ! residual after lin. solver
    real                 :: units(5)              ! mass, length, time,
                                                  ! temperature, angle
    ! Boundary cell type (important for scalars, since they
    ! can have different boundary conditions at the walls)
    integer, allocatable :: bnd_cell_type(:)
    
  end type

  contains 

  include 'Var_Mod/Allocate_New_Only.f90'
  include 'Var_Mod/Allocate_Solution.f90'
  include 'Var_Mod/Allocate_Statistics.f90'
  include 'Var_Mod/Allocate_Gradients.f90'
  include 'Var_Mod/Bnd_Cell_Type.f90'

  end module
