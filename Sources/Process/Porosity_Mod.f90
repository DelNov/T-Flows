!==============================================================================!
  module Porosity_Mod
!------------------------------------------------------------------------------!
!   Porosity                                                                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Control_Mod
  use Work_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!


  !---------------------------------!
  !   Facet (in an STL file) type   !
  !---------------------------------!
  type Facet_Type
    real :: x(3), y(3), z(3)  ! vertex's coordinates
    real :: xc,   yc,   zc    ! facet's centroid
    real :: nx,   ny,   nz    ! normal vector
  end type

  !--------------!
  !   Stl type   !
  !--------------!
  type Stl_Type
    integer                       :: n_facets
    type(Facet_Type), allocatable :: facet(:)
  end type

  !------------------------!
  !   Porous region type   !
  !------------------------!
  type Porous_Region_Type

    ! Convex STL file which defines the porosity shape (must be convex)
    character(SL) :: stl_file

    ! Porosity coefficients
    real :: c1_x = 0.0, c2_x = 0.0;
    real :: c1_y = 0.0, c2_y = 0.0;
    real :: c1_z = 0.0, c2_z = 0.0;

    ! Store porosity for all cells
    logical, allocatable :: cell_porous(:)

  end type

  !-------------------!
  !   Porosity type   !
  !-------------------!
  type Porosity_Type

    ! Pointers to grids for which the porosity is defined
    type(Grid_Type), pointer :: pnt_grid

    ! Number of regions, and regions themlselves
    integer                               :: n_regions = 0
    type(Porous_Region_Type), allocatable :: region(:)

    contains
       procedure          :: Create_Porosity
       procedure, private :: Set_Porosity_In_Cells

  end type

  contains

#   include "Porosity_Mod/Create_Porosity.f90"
#   include "Porosity_Mod/Set_Porosity_In_Cells.f90"

  end module
