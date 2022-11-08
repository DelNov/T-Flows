!==============================================================================!
  module Convert_Mod
!------------------------------------------------------------------------------!
!   Collection of functions used in the Convert program.  In honesty, it was   !
!   introduced to get rid of the Fortran header files with interfaces which,   !
!   in effect, was needed for Intel compiler to work in the debug mode.        !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Convert type   !
  !------------------!
  type Convert_Type
    contains
      procedure :: Allocate_Memory
      procedure :: Calculate_Geometry
      procedure :: Create_Dual
      procedure :: Find_Faces
      procedure :: Find_Parents
      procedure :: Grid_Topology
      procedure :: Guess_Format
      procedure :: Insert_Buildings
      procedure :: Load_Fluent
      procedure :: Load_Gambit
      procedure :: Load_Gmsh
      procedure :: Logo_Con
      procedure :: N_Bnd_Cells_In_Color
      procedure :: N_Edges_In_Bnd_Color
      procedure :: N_Nodes_In_Bnd_Color
      procedure :: N_Sharp_Corners
      procedure :: N_Sharp_Edges
      procedure :: Sort_Face_Nodes
      procedure :: Triangle_Area_Z
  end type

  type(Convert_Type) :: Convert

  contains

#   include "Convert_Mod/Allocate_Memory.f90"
#   include "Convert_Mod/Calculate_Geometry.f90"
#   include "Convert_Mod/Create_Dual.f90"
#   include "Convert_Mod/Find_Faces.f90"
#   include "Convert_Mod/Find_Parents.f90"
#   include "Convert_Mod/Grid_Topology.f90"
#   include "Convert_Mod/Guess_Format.f90"
#   include "Convert_Mod/Insert_Buildings.f90"
#   include "Convert_Mod/Load_Fluent.f90"
#   include "Convert_Mod/Load_Gambit.f90"
#   include "Convert_Mod/Load_Gmsh.f90"
#   include "Convert_Mod/Logo_Con.f90"
#   include "Convert_Mod/N_Bnd_Cells_In_Color.f90"
#   include "Convert_Mod/N_Edges_In_Bnd_Color.f90"
#   include "Convert_Mod/N_Nodes_In_Bnd_Color.f90"
#   include "Convert_Mod/N_Sharp_Corners.f90"
#   include "Convert_Mod/N_Sharp_Edges.f90"
#   include "Convert_Mod/Sort_Face_Nodes.f90"
#   include "Convert_Mod/Triangle_Area_Z.f90"

  end module
