!==============================================================================!
  module Grid_Mod
!------------------------------------------------------------------------------!
!   Grids module is used throughout all programs                               !
!   (that means in "Generate", "Divide", "Convert", "Process".                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Math_Mod
  use Comm_Mod
  use File_Mod
  use Grid_Level_Mod
  use Material_Mod
  use Bnd_Cond_Mod
  use Metis_Options_Mod
  use Sort_Mod
  use Div_Mod
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
    integer :: n_bnd_cond   ! ... boundary conditions
    integer :: n_levels     ! ... multigrid levels

    ! Periodic span
    real :: per_x, per_y, per_z

    ! Minimum, Maximum and total volumes
    real :: min_vol, max_vol, tot_vol

    !-------------------------!
    !  Cell-based variables   !
    !-------------------------!

    ! Cell center coordinates
    real, allocatable :: xc(:), yc(:), zc(:)

    ! Cell volumes
    real, allocatable :: vol(:)

    ! Wall distance - distance from the nearest wall
    real, allocatable :: wall_dist(:)

    ! True if cell is near wall.  Used in Process for some turblence models.
    logical, allocatable :: cell_near_wall(:)

    ! Number of nodes at each cell (determines cell's shape really)
    integer, allocatable :: cells_n_nodes(:)

    ! Number of faces surrounding each cell
    integer, allocatable :: cells_n_faces(:)

    ! Number of cells surrounding each cell
    integer, allocatable :: cells_n_cells(:)

    ! Cells' nodes, faces and neigboring cells
    integer, allocatable :: cells_n(:,:)
    integer, allocatable :: cells_f(:,:)
    integer, allocatable :: cells_c(:,:)

    ! For boundary cells, store corresponding face
    integer, allocatable :: cells_bnd_face(:)

    ! For each cell; type of the boundary condition in a given direction
    integer, allocatable :: cells_bnd_color(:,:)

    ! Coarser levels for the grid
    type(Grid_Level_Type) :: level(MAX_MG_LEVELS)

    !-------------------------!
    !  Face-based variables   !
    !-------------------------!

    ! Number of nodes at each face (determines face's shape really)
    integer, allocatable :: faces_n_nodes(:)

    ! Faces' nodes and neigboring cells
    integer, allocatable :: faces_n(:,:)
    integer, allocatable :: faces_c(:,:)

    ! Periodic faces
    integer, allocatable :: per_faces(:)

    ! Face surface areas (si), total surface (s) 
    ! and distances between cells (di)
    real, allocatable :: sx(:), sy(:), sz(:), s(:)
    real, allocatable :: dx(:), dy(:), dz(:)

    ! Face coordinates 
    real, allocatable :: xf(:), yf(:), zf(:)

    ! Face weight-factors: purely geometrical (f) and
    ! adapted to near wall cells in the fluid phase (fw)
    real, allocatable :: f(:)
    real, allocatable :: fw(:)

    !-------------------------!
    !  Node-based variables   !
    !-------------------------!

    ! Node coordinates
    real, allocatable :: xn(:), yn(:), zn(:)

    type(Material_Type) :: material
    type(Bnd_Cond_Type) :: bnd_cond

    !  Maximum number of cells, boundary cells and faces
    ! (Used for tentative memory allocation in Generator)
    integer :: max_n_nodes
    integer :: max_n_bnd_cells
    integer :: max_n_faces

    ! New numbers for nodes, cells and faces
    integer, allocatable :: new_n(:)
    integer, allocatable :: new_c(:)
    integer, allocatable :: new_f(:)

    ! Number of cells surrounding each node
    integer, allocatable :: nodes_n_cells(:)

    ! List of cells surrounding each node
    integer, allocatable :: nodes_c(:,:)

    !------------------------------------------!
    !   Variables important for parallel run   ! 
    !------------------------------------------!
    type(Comm_Type) :: comm

  end type

  contains

  include 'Grid_Mod/Allocate_Cells.f90'
  include 'Grid_Mod/Allocate_Comm.f90'
  include 'Grid_Mod/Allocate_Faces.f90'
  include 'Grid_Mod/Allocate_Levels.f90'
  include 'Grid_Mod/Allocate_New_Numbers.f90'
  include 'Grid_Mod/Allocate_Nodes.f90'
  include 'Grid_Mod/Bnd_Cond_Name.f90'
  include 'Grid_Mod/Bnd_Cond_Type.f90'
  include 'Grid_Mod/Bnd_Cond_Ranges.f90'
  include 'Grid_Mod/Calculate_Face_Geometry.f90'
  include 'Grid_Mod/Calculate_Global_Volumes.f90'
  include 'Grid_Mod/Calculate_Wall_Distance.f90'
  include 'Grid_Mod/Check_Levels.f90'
  include 'Grid_Mod/Coarsen.f90'
  include 'Grid_Mod/Correction_Periodicity.f90'
  include 'Grid_Mod/Create_Buffers.f90'
  include 'Grid_Mod/Create_Levels.f90'
  include 'Grid_Mod/Decompose.f90'
  include 'Grid_Mod/Estimate_Big_And_Small.f90'
  include 'Grid_Mod/Exchange_Int.f90'
  include 'Grid_Mod/Exchange_Real.f90'
  include 'Grid_Mod/Find_Cells_Faces.f90'
  include 'Grid_Mod/Find_Nodes_Cells.f90'
  include 'Grid_Mod/Find_Periodic_Faces.f90'
  include 'Grid_Mod/Get_C1_And_C2_At_Level.f90'
  include 'Grid_Mod/Load_Cns.f90'
  include 'Grid_Mod/Load_Geo.f90'
  include 'Grid_Mod/Load_Maps.f90'
  include 'Grid_Mod/Print_Bnd_Cond_Info.f90'
  include 'Grid_Mod/Save_Cns.f90'
  include 'Grid_Mod/Save_Geo.f90'
  include 'Grid_Mod/Sort_Cells_By_Index.f90'
  include 'Grid_Mod/Sort_Cells_Smart.f90'
  include 'Grid_Mod/Sort_Faces_By_Index.f90'
  include 'Grid_Mod/Sort_Faces_Smart.f90'

  end module
