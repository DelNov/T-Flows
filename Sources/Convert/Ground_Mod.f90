!==============================================================================!
  module Ground_Mod
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use File_Mod
  use Math_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Facet type   !
  !----------------!
  type Facet_Type
    real :: x(3), y(3), z(3)  ! vertex's coordinates
    real :: area_z
  end type

  !-----------------!
  !   Ground Type   !
  !-----------------!
  type Ground_Type
    integer                       :: n_facets
    type(Facet_Type), allocatable :: facet(:)
  end type

  contains

  include 'Ground_Mod/Read_Stl.f90'

  end module
