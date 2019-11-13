!==============================================================================!
  module Eddies_Mod
!------------------------------------------------------------------------------!
!   Module to define eddy type used in synthetic eddy generation               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Math_Mod
  use Grid_Mod
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !               !
  !   Eddy type   !
  !               !
  !---------------!
  type Eddy_Type
    real :: x, y ,z  ! eddy's position
    real :: length   ! eddy's length
    real :: radius   ! eddy's radius
    real :: sgn      ! eddy's sign (sense of rotation)
  end type

  !-----------------!
  !                 !
  !   Eddies type   !
  !                 !
  !-----------------!
  type Eddies_Type
    type(Grid_Type),  pointer :: pnt_grid  ! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  ! flow for which it is defined
    character(len=80)         :: bc_name   ! boundary condition name

    integer                      :: n_eddies
    type(Eddy_Type), allocatable :: eddy(:)

    real :: max_radius

    integer              :: n_bnd_cells
    integer, allocatable :: bnd_cell(:)
  end type

  contains

  include 'Eddies_Mod/Allocate.f90'
  include 'Eddies_Mod/Place_Eddy.f90'
  include 'Eddies_Mod/Superimpose.f90'

  end module
