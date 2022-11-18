!==============================================================================!
!   Introduce new types to be used with User_Mod                               !
!==============================================================================!

  ! Face in an STL file
  type Facet_Type
    real :: x(3), y(3), z(3)  ! vertex's coordinates
    real :: xc,   yc,   zc    ! facet's centroid
    real :: nx, ny, nz        ! normal vector
  end type

  ! Forrest type
  type Forrest_Type
    integer                       :: n_facets
    type(Facet_Type), allocatable :: facet(:)
  end type

  ! Porosity
  logical, allocatable :: porosity(:)

