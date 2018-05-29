!==============================================================================!
  module Grid_Mod
!------------------------------------------------------------------------------!
!   Grids module is used throughout all programs                               !
!   (that means in "Generate", "Divide", "Convert", "Process".                 !
!------------------------------------------------------------------------------!
  use Material_Mod
  use Bnd_Cond_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !               !
  !   Grid type   !
  !               !
  !---------------!
  type Grid_Type

    ! Number of ...
    integer :: n_nodes      ! ... nodes
    integer :: n_cells      ! ... cells
    integer :: n_faces      ! ... faces
    integer :: n_bnd_cells  ! ... boundary cells
    integer :: n_per_faces  ! ... periodic faces (shadows)
    integer :: n_materials  ! ... materials
    integer :: n_bnd_cond   ! ... boundary conditions
    integer :: n_copy       ! ... copy cells and faces
    integer :: n_sh         ! ... shadow faces           

    !-------------------------!
    !  Cell-based variables   !
    !-------------------------!

    ! Cell center coordinates
    real, allocatable :: xc(:), yc(:), zc(:)  

    ! Cell volumes
    real, allocatable :: vol(:)
  
    ! General cell size (max(dx,dy,dz) or maybe volume/area?)
    real, allocatable :: delta(:)
    
    ! Wall distance - distance from the nearest wall
    real, allocatable :: wall_dist(:)

    ! True if cell is near wall.  Used in Process for some turblence models.
    logical, allocatable :: cell_near_wall(:)
    
    ! Cells' nodes and neigboring cells
    integer, allocatable :: cells_n(:,:)      
    integer, allocatable :: cells_c(:,:)

    ! Number of nodes at each cell (determines cell's shape really)
    integer, allocatable :: cells_n_nodes(:)

    ! For each cell; type of the boundary condition in a given direction
    integer, allocatable :: cells_bnd_color(:,:)

    ! Material for each cell 
    integer, allocatable :: material(:)

    !-------------------------!
    !  Face-based variables   !
    !-------------------------!

    ! Number of nodes at each face (determines face's shape really)
    integer, allocatable :: faces_n_nodes(:)

    ! Faces' nodes and neigboring cells
    integer, allocatable :: faces_n(:,:)
    integer, allocatable :: faces_c(:,:)

    ! Face surface areas (si), total surface (s) 
    ! and distances between cells (di)
    real, allocatable :: sx(:), sy(:), sz(:), s(:)
    real, allocatable :: dx(:), dy(:), dz(:)

    ! Face coordinates 
    real, allocatable :: xf(:), yf(:), zf(:)
    
    ! Face weight-factor
    real, allocatable :: f(:)
    
    !-------------------------!
    !  Node-based variables   !
    !-------------------------!

    ! Node coordinates
    real, allocatable :: xn(:), yn(:), zn(:)
    
    type(Material_Type), allocatable :: materials(:)
    type(Bnd_Cond_Type)              :: bnd_cond

    !  Maximum number of cells, boundary cells and faces
    ! (Used for tentative memory allocation in Generator)
    integer :: max_n_nodes
    integer :: max_n_bnd_cells
    integer :: max_n_faces

  end type

  contains
 
  include 'Grid_Mod/Allocate_Cells.f90'
  include 'Grid_Mod/Allocate_Faces.f90'
  include 'Grid_Mod/Allocate_Nodes.f90'
  include 'Grid_Mod/Bnd_Cond_Type.f90'
  include 'Grid_Mod/Print_Bnd_Cond_Info.f90'
  include 'Grid_Mod/Sort_Cells_By_Index.f90'
  include 'Grid_Mod/Sort_Faces_By_Index.f90'

  end module
