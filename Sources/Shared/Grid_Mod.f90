#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"

!==============================================================================!
  module Grid_Mod
!------------------------------------------------------------------------------!
!>  The Grid_Mod module is a cornerstone in the T-Flows program, serving as the
!>  central repository for data structures and operations related to the
!>  computational grid. It is utilized across all sub-programs of T-Flows
!>  including Convert, Generate, Divide, and Process. This module's primary
!>  construct, Grid_Type, houses the necessary data for representing a
!>  computational grid and provides a comprehensive suite of methods for
!>  grid manipulation. These methods encompass a wide range of functionalities
!>  such as geometrical calculations, connectivity analysis, mesh division,
!>  and grid data export in various formats. The module's design emphasizes
!>  versatility and robustness, enabling efficient handling of diverse grid
!>  configurations and facilitating complex fluid dynamics simulations.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Memory_Mod
  use Omp_Mod
  use Profiler_Mod
  use File_Mod
  use Region_Mod
  use Vtk_Mod
  use Metis_Mod
  use Sort_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !               !
  !   Grid type   !
  !               !
  !---------------!
  !> Grid_Type is used throughout all sub-programs ("Generate", "Divide",
  !> "Convert", "Process") and holds the data which defines a grid as well
  !> as the procedures needed to manipulate the grids.
  type Grid_Type

    ! Stores the name of this domain
    character(SL) :: name                  !! name of the grid
    logical       :: polyhedral = .false.  !! true if grid is polyhedral

    ! Number of ...
    integer :: n_nodes       = 0  !! number of nodes in a grid
    integer :: n_cells       = 0  !! number of cells in a grid
    integer :: n_faces       = 0  !! number of faces in a grid
    integer :: n_bnd_cells   = 0  !! number of boundary cells
    integer :: n_bnd_regions = 0  !! number of boundary conditions
    integer :: n_regions     = 0  !! number of all conditions (bnd, inside, per)
    integer :: n_shadows     = 0  !! number of shadow (periodic) faces
    integer :: n_edges       = 0  !! number of edges (when creating dual grid)
    integer :: per_x_reg     = 0  !! periodic x region
    integer :: per_y_reg     = 0  !! periodic y region
    integer :: per_z_reg     = 0  !! periodic z region

    ! Rank (in case of simulations with multiple domains)
    integer :: rank = 0  !! grid rank (it's number)

    ! Periodic span
    real :: per_x, per_y, per_z  !! periodic distance in coordinate direction

    ! Minimum, Maximum and total volumes
    real :: min_vol  !! volume of the smallest cell in the grid
    real :: max_vol  !! volume of the biggest cell in the grid
    real :: tot_vol  !! total volume of the grid

    !-------------------------!
    !  Cell-based variables   !
    !-------------------------!

    ! Cell center coordinates
    real, allocatable :: xc(:), yc(:), zc(:)  !! cell center coordinates

    ! Cell volumes
    real, allocatable :: vol(:)  !! cell volumes

    ! Cells' tensor of inertia
    real, allocatable :: ixx(:), iyy(:), izz(:), ixy(:), ixz(:), iyz(:)

    ! Fractional cell volumes around faces
    real, allocatable :: dv1(:)
    real, allocatable :: dv2(:)

    ! Wall distance - distance from the nearest wall
    real, allocatable :: wall_dist(:)  !! distance from the nearest wall

    ! True if cell is near wall.  Used in Process for some turblence models.
    logical, allocatable :: cell_near_wall(:)  !! true if cell is near wall

    ! Number of nodes at each cell (determines cell's shape really)
    integer, allocatable :: cells_n_nodes(:)  !! number of nodes in each cell

    ! Number of faces surrounding each cell
    integer, allocatable :: cells_n_faces(:)  !! number of faces in each cell

    ! Number of cells surrounding each cell
    integer, allocatable :: cells_n_cells(:)  !! number of cells around cell

    ! Cells' nodes, faces, and neigboring cells
    integer, allocatable :: cells_n(:,:)  !! list of cells' nodes
    integer, allocatable :: cells_f(:,:)  !! list of cells' faces
    integer, allocatable :: cells_c(:,:)  !! list of cells' neighbouring cells
                                          !! (needed in Generate and Convert)

    ! Weights for interpolation from nodes
    real, allocatable :: weight_n2c(:,:)

    ! For boundary cells, store corresponding face
    integer, allocatable :: cells_bnd_face(:)

    ! For each cell; type of the boundary condition in a given direction
    integer, allocatable :: cells_bnd_region(:,:)

    ! For each cell: to each porous region it belongs
    ! (Introduced with Insert_Buildings in Convert)
    integer, allocatable :: por(:)  !! porous region to which each cell belongs

    !-------------------------!
    !  Face-based variables   !
    !-------------------------!

    ! Number of nodes at each face (determines face's shape really)
    integer, allocatable :: faces_n_nodes(:)  !! number of nodes for each face

    ! Faces' nodes, neigboring cells and shadows
    integer, allocatable :: faces_n(:,:)  !! faces' nodes
    integer, allocatable :: faces_c(:,:)  !! faces' cells (the two around it)
    integer, allocatable :: faces_s(:)    !! faces' shadows (periodic pairs)

    ! Face surface areas (si), total surface (s) 
    ! and distances between cells (di)
    real, allocatable :: sx(:), sy(:), sz(:), s(:)
    real, allocatable :: dx(:), dy(:), dz(:), d(:)

    ! Face coordinates
    real, allocatable :: xf(:), yf(:), zf(:)  !! face center coordinate

    ! Vectors connecting face center with face cell centers connection
    real, allocatable :: rx(:), ry(:), rz(:)

    ! Face weight-factors: purely geometrical (f) and
    ! adapted to near wall cells in the fluid phase (fw)
    real, allocatable :: f (:)
    real, allocatable :: fw(:)

    ! For each face, which neighbours are to each other cells which meet there.
    ! The name derives from "face cell to cell".  Used in Generate only
    integer, allocatable :: face_c_to_c(:,:)

    !-------------------------!
    !  Node-based variables   !
    !-------------------------!

    ! Node coordinates
    real, allocatable :: xn(:), yn(:), zn(:)  !! node coordinate

    type(Region_Type) :: region  !! boundary condition regions

    ! Twin (periodic) nodes, and node_positioned are used only in Generate
    integer, allocatable :: twin_n(:,:)         !! for each node, list of
                                                !! twins, periodic nodes
    logical, allocatable :: node_positioned(:)  !! true if node is positioned

    ! New numbers for nodes, cells, faces and edges
    integer, allocatable :: new_n(:)  !! new node numbers (when renumbering)
    integer, allocatable :: new_c(:)  !! new cell numbers (when renumbering)
    integer, allocatable :: new_f(:)  !! new face numbers (when renumbering)
    integer, allocatable :: new_e(:)  !! new edge numbers
                                      !! (used for dual grid creation)
    ! Old numbers for cells and faces
    integer, allocatable :: old_n(:)  !! old node numbers (when renumbering)
    integer, allocatable :: old_c(:)  !! old cell numbers (when renumbering)
    integer, allocatable :: old_f(:)  !! old face numbers (when renumbering)

    ! Number of cells surrounding each node
    integer, allocatable :: nodes_n_cells(:)  !! number of cells surrounding
                                              !! each node
    ! List of cells surrounding each node ...
    ! ... and weights for cell to node interpolation
    integer, allocatable :: nodes_c(:,:)     !! cells surrounding each node
    real,    allocatable :: weight_c2n(:,:)  !! weights for cell to node
                                             !! interpolation

    ! Edge-base variables
    integer, allocatable :: edges_n (:,:)  !! edges' nodes
    integer, allocatable :: edges_bc(:,:)  !! edges' boundary conditions
    integer, allocatable :: edges_fb(:,:)  !! edges' faces on boundaries

    !-------------------------------------------!
    !   Communication class for parallel runs   !
    !-------------------------------------------!
    type(Comm_Type) :: Comm  !! module for MPI communication, local to grid

    !---------------------------------!
    !   OMP class for manycore runs   !
    !---------------------------------!
    type(Omp_Type) :: Omp  !! used in OMP vectorization

    contains
      procedure :: Are_Nodes_Twins
      procedure :: Allocate_Cells
      procedure :: Allocate_Faces
      procedure :: Allocate_Nodes
      procedure :: Allocate_Regions
      procedure :: Bnd_Cond_Name_At_Cell
      procedure :: Bnd_Cond_Name_At_Face
      procedure :: Bnd_Cond_Type
      procedure :: Bounding_Box
      procedure :: Calculate_Cell_Centers
      procedure :: Calculate_Cell_Inertia
      procedure :: Calculate_Cell_Volumes
      procedure :: Calculate_Face_Centers
      procedure :: Calculate_Face_Geometry
      procedure :: Calculate_Face_Interpolation
      procedure :: Calculate_Face_Surfaces
      procedure :: Calculate_Global_Volumes
      procedure :: Calculate_Wall_Distance
      procedure :: Calculate_Weights_Cells_To_Nodes
      procedure :: Calculate_Weights_Nodes_To_Cells
      procedure :: Check_Cells_Closure
      procedure :: Correct_Face_Surfaces
      procedure :: Decompose
      procedure :: Determine_Regions_Ranges
      procedure :: Determine_Threads
      procedure :: Exchange_Cells_Int
      procedure :: Exchange_Cells_Log
      procedure :: Exchange_Cells_Real
      procedure :: Exchange_Inside_Cells_Real
      procedure :: Face_Normal
      procedure :: Faces_Surface
      procedure :: Find_Cells_Faces
      procedure :: Find_Nodes_Cells
      procedure :: Form_Cells_Comm
      procedure :: Form_Maps_For_Backup
      procedure :: Initialize_New_Numbers
      procedure :: Is_Face_In_Cell
      procedure :: Is_Point_In_Cell
      procedure :: Load_And_Prepare_For_Processing
      procedure :: Load_Cfn
      procedure :: Load_Dim
      procedure :: Merge_Duplicate_Nodes
      procedure :: Print_Regions_List
      procedure :: Print_Grid_Statistics
      procedure :: Save_Cfn
      procedure :: Save_Dim
      procedure :: Save_Debug_Vtu
      procedure :: Save_Vtk_Cell
      procedure :: Save_Vtk_Face
      procedure :: Save_Vtu_Cells
      procedure :: Save_Vtu_Edges
      procedure :: Save_Vtu_Faces
      procedure :: Search_Coordinate_Clusters
      procedure :: Sort_Cells_By_Coordinates
      procedure :: Sort_Cells_By_Thread
      procedure :: Sort_Faces_By_Index
      procedure :: Sort_Faces_By_Region
      procedure :: Sort_Nodes_By_Coordinates
      procedure :: Write_Template_Control_File

  end type

  contains

#   include "Grid_Mod/Boundary/Bnd_Cond_Name_At_Cell.f90"
#   include "Grid_Mod/Boundary/Bnd_Cond_Name_At_Face.f90"
#   include "Grid_Mod/Boundary/Bnd_Cond_Type.f90"
#   include "Grid_Mod/Calculate/Bounding_Box.f90"
#   include "Grid_Mod/Calculate/Calculate_Cell_Centers.f90"
#   include "Grid_Mod/Calculate/Calculate_Cell_Inertia.f90"
#   include "Grid_Mod/Calculate/Calculate_Cell_Volumes.f90"
#   include "Grid_Mod/Calculate/Calculate_Face_Centers.f90"
#   include "Grid_Mod/Calculate/Calculate_Face_Geometry.f90"
#   include "Grid_Mod/Calculate/Calculate_Face_Interpolation.f90"
#   include "Grid_Mod/Calculate/Calculate_Faces_Surface.f90"
#   include "Grid_Mod/Calculate/Calculate_Face_Surfaces.f90"
#   include "Grid_Mod/Calculate/Calculate_Global_Volumes.f90"
#   include "Grid_Mod/Calculate/Calculate_Wall_Distance.f90"
#   include "Grid_Mod/Calculate/Calculate_Weights_Cells_To_Nodes.f90"
#   include "Grid_Mod/Calculate/Calculate_Weights_Nodes_To_Cells.f90"
#   include "Grid_Mod/Connectivity/Are_Nodes_Twins.f90"
#   include "Grid_Mod/Connectivity/Check_Cells_Closure.f90"
#   include "Grid_Mod/Connectivity/Correct_Face_Surfaces.f90"
#   include "Grid_Mod/Connectivity/Determine_Regions_Ranges.f90"
#   include "Grid_Mod/Connectivity/Face_Normal.f90"
#   include "Grid_Mod/Connectivity/Find_Cells_Faces.f90"
#   include "Grid_Mod/Connectivity/Find_Nodes_Cells.f90"
#   include "Grid_Mod/Connectivity/Initialize_New_Numbers.f90"
#   include "Grid_Mod/Connectivity/Is_Face_In_Cell.f90"
#   include "Grid_Mod/Connectivity/Is_Point_In_Cell.f90"
#   include "Grid_Mod/Connectivity/Merge_Duplicate_Nodes.f90"
#   include "Grid_Mod/Connectivity/Search_Coordinate_Clusters.f90"
#   include "Grid_Mod/Input_Output/Load_And_Prepare_For_Processing.f90"
#   include "Grid_Mod/Input_Output/Load_Cfn.f90"
#   include "Grid_Mod/Input_Output/Load_Dim.f90"
#   include "Grid_Mod/Input_Output/Print_Grid_Statistics.f90"
#   include "Grid_Mod/Input_Output/Print_Regions_List.f90"
#   include "Grid_Mod/Input_Output/Save_Cfn.f90"
#   include "Grid_Mod/Input_Output/Save_Debug_Vtu.f90"
#   include "Grid_Mod/Input_Output/Save_Dim.f90"
#   include "Grid_Mod/Input_Output/Save_Vtk_Cell.f90"
#   include "Grid_Mod/Input_Output/Save_Vtk_Face.f90"
#   include "Grid_Mod/Input_Output/Save_Vtu_Cells.f90"
#   include "Grid_Mod/Input_Output/Save_Vtu_Edges.f90"
#   include "Grid_Mod/Input_Output/Save_Vtu_Faces.f90"
#   include "Grid_Mod/Input_Output/Write_Template_Control_File.f90"
#   include "Grid_Mod/Memory/Allocate_Cells.f90"
#   include "Grid_Mod/Memory/Allocate_Faces.f90"
#   include "Grid_Mod/Memory/Allocate_Nodes.f90"
#   include "Grid_Mod/Memory/Allocate_Regions.f90"
#   include "Grid_Mod/Parallel/Decompose.f90"
#   include "Grid_Mod/Parallel/Determine_Threads.f90"
#   include "Grid_Mod/Parallel/Exchange_Cells_Int.f90"
#   include "Grid_Mod/Parallel/Exchange_Cells_Log.f90"
#   include "Grid_Mod/Parallel/Exchange_Cells_Real.f90"
#   include "Grid_Mod/Parallel/Exchange_Inside_Cells_Real.f90"
#   include "Grid_Mod/Parallel/Form_Cells_Comm.f90"
#   include "Grid_Mod/Parallel/Form_Maps_For_Backup.f90"
#   include "Grid_Mod/Sorting/Sort_Cells_By_Coordinates.f90"
#   include "Grid_Mod/Sorting/Sort_Cells_By_Thread.f90"
#   include "Grid_Mod/Sorting/Sort_Faces_By_Index.f90"
#   include "Grid_Mod/Sorting/Sort_Faces_By_Region.f90"
#   include "Grid_Mod/Sorting/Sort_Nodes_By_Coordinates.f90"

  end module
