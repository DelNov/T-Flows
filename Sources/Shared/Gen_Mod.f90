!==============================================================================!
  module Gen_Mod
!------------------------------------------------------------------------------!
!   Global variable definitions for the mesh generator.                        !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Twin nodes
  integer, allocatable :: twin_n(:,:)

  ! New numbers for nodes, cells and faces
  integer, allocatable :: new_n(:)
  integer, allocatable :: new_c(:)
  integer, allocatable :: new_f(:)

  ! For each face, which neighbours are to each other cells which meet there
  ! The name derives from "face cell to cell"
  integer, allocatable :: face_c_to_c(:,:)
                                                
  ! Variables for refinement
  integer, parameter   :: ELIPSOID  = 3
  integer, parameter   :: RECTANGLE = 4
  integer, parameter   :: PLANE     = 5
  integer              :: n_refine_levels
  integer, allocatable :: level(:)                ! refinement level
  integer, allocatable :: n_refined_regions(:)    ! number of refin. regions
  real,    allocatable :: refined_regions(:,:,:)  ! levels, regions
  logical, allocatable :: cell_marked(:)          ! true if cell markered 

  ! Variables for smoothing
  integer              :: n_smoothing_regions   
  logical, allocatable :: smooth_in_x(:),  &
                          smooth_in_y(:),  &
                          smooth_in_z(:)   
  integer, allocatable :: smooth_iters(:)   
  real,    allocatable :: smooth_regions(:,:), smooth_relax(:)

  ! Periodic and copy boundaries
  integer              :: n_periodic_cond,      n_copy_cond
  integer, allocatable ::   periodic_cond(:,:),   copy_cond(:,:)

  ! First and last wall face - should speed up the cal-
  ! culation of distance to walls to a certain degree
  integer :: first_wall_face, last_wall_face

end module
