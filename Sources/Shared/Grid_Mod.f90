#include "../Shared/Assert.h90"

!==============================================================================!
  module Grid_Mod
!------------------------------------------------------------------------------!
!   Grids module is used throughout all programs                               !
!   (that means in "Generate", "Divide", "Convert", "Process".                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Assert_Mod
  use Profiler_Mod
  use File_Mod
  use Boundary_Mod
  use Vtk_Mod
  use Metis_Mod
  use Sort_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Interfaces]---------------------------------!
  interface
    include '../Shared/Adjust_First_Dim.h90'
  end interface
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: INWARDS  = -1  ! for faces wrapping cells, can ...
  integer, parameter :: OUTWARDS = +1  ! ... be pointing inwards or outwards
!==============================================================================!

  !---------------!
  !               !
  !   Grid type   !
  !               !
  !---------------!
  type Grid_Type

    ! Stores the name of this domain
    character(SL) :: name
    logical       :: polyhedral = .false.

    ! Number of ...
    integer :: n_nodes     = 0  ! ... nodes
    integer :: n_cells     = 0  ! ... cells
    integer :: n_faces     = 0  ! ... faces
    integer :: n_bnd_cells = 0  ! ... boundary cells
    integer :: n_bnd_cond  = 0  ! ... boundary conditions
    integer :: n_shadows   = 0  ! ... shadow faces
    integer :: n_edges     = 0  ! ... edges (needed to create dual grid)

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

    ! Cells' tensor of inertia
    real, allocatable :: ixx(:), iyy(:), izz(:), ixy(:), ixz(:), iyz(:)

    ! Fractional cell volumes around faces
    real, allocatable :: dv1(:)
    real, allocatable :: dv2(:)

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

    ! Cells' nodes, faces, and neigboring cells
    integer, allocatable :: cells_n(:,:)
    integer, allocatable :: cells_f(:,:)
    integer, allocatable :: cells_c(:,:)

    ! Orientation (INWARDS, or OUTWARDS) of faces sourrounding the cells
    integer, allocatable :: cells_f_orient(:,:)

    ! Weights for interpolation from nodes
    real, allocatable :: weight_n2c(:,:)

    ! For boundary cells, store corresponding face
    integer, allocatable :: cells_bnd_face(:)

    ! For each cell; type of the boundary condition in a given direction
    integer, allocatable :: cells_bnd_color(:,:)

    !-------------------------!
    !  Face-based variables   !
    !-------------------------!

    ! Number of nodes at each face (determines face's shape really)
    integer, allocatable :: faces_n_nodes(:)

    ! Faces' nodes, neigboring cells and shadows
    integer, allocatable :: faces_n(:,:)
    integer, allocatable :: faces_c(:,:)
    integer, allocatable :: faces_s(:)

    ! Face surface areas (si), total surface (s) 
    ! and distances between cells (di)
    real, allocatable :: sx(:), sy(:), sz(:), s(:)
    real, allocatable :: dx(:), dy(:), dz(:), d(:)

    ! Face coordinates
    real, allocatable :: xf(:), yf(:), zf(:)

    ! Face-based interserction with surface
    real, allocatable :: xs(:), ys(:), zs(:)

    ! Vectors connecting face center with face cell centers connection
    real, allocatable :: rx(:), ry(:), rz(:)

    ! Face weight-factors: purely geometrical (f) and
    ! adapted to near wall cells in the fluid phase (fw)
    real, allocatable :: f (:)
    real, allocatable :: fw(:)

    !-------------------------!
    !  Node-based variables   !
    !-------------------------!

    ! Node coordinates
    real, allocatable :: xn(:), yn(:), zn(:)

    type(Boundary_Type) :: bnd_cond

    !  Maximum number of cells, boundary cells and faces
    ! (Used for tentative memory allocation in Generator)
    integer :: max_n_nodes
    integer :: max_n_bnd_cells
    integer :: max_n_faces

    ! New numbers for nodes, cells, faces and edges
    integer, allocatable :: new_n(:)
    integer, allocatable :: new_c(:)
    integer, allocatable :: new_f(:)
    integer, allocatable :: new_e(:)  ! used for dual grid creation

    ! Old numbers for cells and faces
    integer, allocatable :: old_n(:)
    integer, allocatable :: old_c(:)
    integer, allocatable :: old_f(:)

    ! Number of cells surrounding each node
    integer, allocatable :: nodes_n_cells(:)

    ! List of cells surrounding each node ...
    ! ... and weights for cell to node interpolation
    integer, allocatable :: nodes_c(:,:)
    real,    allocatable :: weight_c2n(:,:)

    ! Edge-base variables
    integer, allocatable :: edges_n (:,:)  ! edges' nodes
    integer, allocatable :: edges_bc(:,:)  ! edges' boundary conditions
    integer, allocatable :: edges_fb(:,:)  ! edges' faces on boundaries

    !------------------------------------------!
    !   Communication class for parallel run   !
    !------------------------------------------!
    type(Comm_Type) :: Comm

    ! User arrays.  I am neither sure if this is the ...
    ! ... best place for them nor do I need them at all?
    integer           :: n_user_arrays
    real, allocatable :: user_array(:,:)

    contains
      procedure :: Allocate_Cells
      procedure :: Allocate_Faces
      procedure :: Allocate_Nodes
      procedure :: Bnd_Cond_Name
      procedure :: Bnd_Cond_Type
      procedure :: Bnd_Cond_Ranges
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
      procedure :: Decompose
      procedure :: Determine_Face_Orientation
      procedure :: Exchange_Cells_Int
      procedure :: Exchange_Cells_Log
      procedure :: Exchange_Cells_Real
      procedure :: Face_In_Cell
      procedure :: Face_Normal
      procedure :: Faces_Surface
      procedure :: Find_Cells_Faces
      procedure :: Find_Nodes_Cells
      procedure :: Form_Cells_Comm
      procedure :: Form_Maps
      procedure :: Initialize_New_Numbers
      procedure :: Is_Point_In_Cell
      procedure :: Load_Cfn
      procedure :: Load_Dim
      procedure :: Merge_Duplicate_Nodes
      procedure :: Print_Bnd_Cond_List
      procedure :: Print_Grid_Statistics
      procedure :: Save_Cfn
      procedure :: Save_Dim
      procedure :: Save_Debug_Vtu
      procedure :: Save_Vtu_Cells
      procedure :: Save_Vtu_Edges
      procedure :: Save_Vtu_Faces
      procedure :: Sort_Cells_By_Index
      procedure :: Sort_Cells_Smart
      procedure :: Sort_Faces_By_Index
      procedure :: Sort_Faces_Smart
      procedure :: Write_Template_Control_File

  end type

  contains

#   include "Grid_Mod/Allocate_Cells.f90"
#   include "Grid_Mod/Allocate_Faces.f90"
#   include "Grid_Mod/Allocate_Nodes.f90"
#   include "Grid_Mod/Bnd_Cond_Name.f90"
#   include "Grid_Mod/Bnd_Cond_Type.f90"
#   include "Grid_Mod/Bnd_Cond_Ranges.f90"
#   include "Grid_Mod/Bounding_Box.f90"
#   include "Grid_Mod/Calculate_Cell_Centers.f90"
#   include "Grid_Mod/Calculate_Cell_Inertia.f90"
#   include "Grid_Mod/Calculate_Cell_Volumes.f90"
#   include "Grid_Mod/Calculate_Face_Centers.f90"
#   include "Grid_Mod/Calculate_Face_Geometry.f90"
#   include "Grid_Mod/Calculate_Face_Interpolation.f90"
#   include "Grid_Mod/Calculate_Face_Surfaces.f90"
#   include "Grid_Mod/Calculate_Faces_Surface.f90"
#   include "Grid_Mod/Calculate_Global_Volumes.f90"
#   include "Grid_Mod/Calculate_Wall_Distance.f90"
#   include "Grid_Mod/Calculate_Weights_Cells_To_Nodes.f90"
#   include "Grid_Mod/Calculate_Weights_Nodes_To_Cells.f90"
#   include "Grid_Mod/Decompose.f90"
#   include "Grid_Mod/Determine_Face_Orientation.f90"
#   include "Grid_Mod/Exchange_Cells_Int.f90"
#   include "Grid_Mod/Exchange_Cells_Log.f90"
#   include "Grid_Mod/Exchange_Cells_Real.f90"
#   include "Grid_Mod/Face_In_Cell.f90"
#   include "Grid_Mod/Face_Normal.f90"
#   include "Grid_Mod/Find_Cells_Faces.f90"
#   include "Grid_Mod/Find_Nodes_Cells.f90"
#   include "Grid_Mod/Form_Cells_Comm.f90"
#   include "Grid_Mod/Form_Maps.f90"
#   include "Grid_Mod/Initialize_New_Numbers.f90"
#   include "Grid_Mod/Is_Point_In_Cell.f90"
#   include "Grid_Mod/Load_Cfn.f90"
#   include "Grid_Mod/Load_Dim.f90"
#   include "Grid_Mod/Merge_Duplicate_Nodes.f90"
#   include "Grid_Mod/Print_Bnd_Cond_List.f90"
#   include "Grid_Mod/Print_Grid_Statistics.f90"
#   include "Grid_Mod/Save_Cfn.f90"
#   include "Grid_Mod/Save_Dim.f90"
#   include "Grid_Mod/Save_Debug_Vtu.f90"
#   include "Grid_Mod/Save_Vtu_Cells.f90"
#   include "Grid_Mod/Save_Vtu_Edges.f90"
#   include "Grid_Mod/Save_Vtu_Faces.f90"
#   include "Grid_Mod/Sort_Cells_By_Index.f90"
#   include "Grid_Mod/Sort_Cells_Smart.f90"
#   include "Grid_Mod/Sort_Faces_By_Index.f90"
#   include "Grid_Mod/Sort_Faces_Smart.f90"
#   include "Grid_Mod/Write_Template_Control_File.f90"

  end module
