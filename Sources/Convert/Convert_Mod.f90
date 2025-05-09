#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Convert_Mod
!------------------------------------------------------------------------------!
!>  Convert_Mod is a module designed for the Convert program in the T-Flows
!>  project. It provides a collection of procedures for various conversion and
!>  geometry processing tasks. The module integrates functions for handling
!>  different mesh formats and geometrical operations, enhancing the
!>  versatility and capability of the Convert program.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Pattern_Mod
  use Stl_Mod
# ifdef __INTEL_COMPILER
  use Ifport              ! Intel's module for fseek and ftell
# endif
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Convert type   !
  !------------------!
  !> Convert_Type encapsulates a number of procedures for various
  !> tasks such as geometry calculation, grid topology analysis,
  !> and handling different mesh formats.
  type Convert_Type

    contains
      procedure, private :: Allocate_Memory
      procedure          :: Calculate_Geometry
      procedure          :: Create_Dual
      procedure          :: Find_Boundary_Faces
      procedure          :: Find_Inside_Faces
      procedure          :: Find_Parents
      procedure          :: Guess_Format
      procedure          :: Insert_Buildings
      procedure          :: Load_Fluent
      procedure          :: Load_Forrest
      procedure          :: Load_Gambit
      procedure          :: Load_Gmsh
      procedure          :: Load_Obj
      procedure          :: Logo_Con
      procedure, private :: N_Bnd_Cells_In_Region
      procedure, private :: N_Edges_In_Region
      procedure, private :: N_Nodes_In_Region
      procedure, private :: N_Sharp_Corners
      procedure, private :: N_Sharp_Edges
      procedure          :: Save_Fluent
      procedure, private :: Sort_Face_Nodes
  end type

  ! Singleton Convert object
  type(Convert_Type) :: Convert  !! singleton Convert object for easy
                                 !! access to Conver_Type's functions

  !----------------------------------!
  !   Declarations for cell shapes   !
  !   (Was "Cells_Face_Nodes.f90")   !
  !----------------------------------!

  !====================================================================!
  !                                 hex                                !
  !                                                                    !
  !                 8-----------7             8-----------7            !
  !                /|          /|            /|          /|            !
  !               /           / |           /    (6)    / |            !
  !              /  |    (3) /  |          /  |        /  |            !
  !             5-----------6   |         5-----------6   |            !
  !             |(4)|       |(2)|         |   |       |   |            !
  !             |   4- - - -|- -3         |   4- - - -|- -3            !
  !             |  / (1)    |  /          |  /        |  /             !
  !             |           | /           |      (5)  | /              !
  !             |/          |/            |/          |/               !
  !             1-----------2             1-----------2                !
  !                                                                    !
  !             (3) and (4) are behind                                 !
  !--------------------------------------------------------------------!
  integer, dimension(6,4), parameter :: HEX   &
      = transpose(reshape( (/ 1, 5, 6, 2,  &
                              2, 6, 7, 3,  &
                              3, 7, 8, 4,  &
                              4, 8, 5, 1,  &
                              1, 2, 3, 4,  &
                              5, 8, 7, 6  /), (/4, 6/) ))

  !====================================================================!
  !                                 tet                                !
  !                                                                    !
  !                     4-----3                   4-----3              !
  !                    / \  .'|                  / \  .'|              !
  !                   /(4)\(3)|                 /   \'  |              !
  !                  /  .' \  |                /  .' \  |              !
  !                 / .' (1)\ |               / .(2)  \ |              !
  !                /.'       \|              /.'       \|              !
  !               1-----------2             1-----------2              !
  !                                                                    !
  !               (1) and (4) are behind                               !
  !--------------------------------------------------------------------!
  integer, dimension(6,4), parameter :: TET   &
      = transpose(reshape( (/ 1, 2, 3,-1,  &
                              1, 4, 2,-1,  &
                              2, 4, 3,-1,  &
                              3, 4, 1,-1,  &
                             -1,-1,-1,-1,  &
                             -1,-1,-1,-1  /), (/4, 6/) ))

  !====================================================================!
  !                                 wed                                !
  !                                                                    !
  !                   6.                        6.                     !
  !                  /| `.                     /| `.                   !
  !                 /     `.                  / (5) `.                 !
  !                /  |     `.               /  |     `.               !
  !               4-----------5             4-----------5              !
  !               |(3)|  (2)  |             |   |       |              !
  !               |   3.      |             |   3.      |              !
  !               |  / (1)    |             |  /  `.    |              !
  !               |       `.  |             |   (4) `.  |              !
  !               |/        `.|             |/        `.|              !
  !               1-----------2             1-----------2              !
  !                                                                    !
  !               (2), (3) and (4) are behind                          !
  !                                                                    !
  !--------------------------------------------------------------------!
  integer, dimension(6,4), parameter :: WED   &
      = transpose(reshape( (/ 1, 4, 5, 2,  &
                              2, 5, 6, 3,  &
                              3, 6, 4, 1,  &
                              1, 2, 3,-1,  &
                              4, 6, 5,-1,  &
                             -1,-1,-1,-1  /), (/4, 6/) ))

  !====================================================================!
  !                                pyr                                 !
  !                                                                    !
  !                                .5.                                 !
  !                              .'/ \`.                               !
  !                            .' /   \ `.                             !
  !                          .'  / (4) \  `.                           !
  !                        .'(5)/       \(3)`.                         !
  !                      4' - -/- - - - -\- - `3                       !
  !                      |    /    (2)    \    |                       !
  !                      |   /             \   |                       !
  !                      |  /      (1)      \  |                       !
  !                      | /                 \ |                       !
  !                      |/                   \|                       !
  !                      1---------------------2                       !
  !                                                                    !
  !                      (4) is behind                                 !
  !                                                                    !
  !--------------------------------------------------------------------!
  integer, dimension(6,4), parameter :: PYR   &
      = transpose(reshape( (/ 1, 2, 3, 4,  &
                              1, 5, 2,-1,  &
                              2, 5, 3,-1,  &
                              3, 5, 4,-1,  &
                              4, 5, 1,-1,  &
                             -1,-1,-1,-1  /), (/4, 6/) ))

  contains

#   include "Convert_Mod/Allocate_Memory.f90"
#   include "Convert_Mod/Calculate_Geometry.f90"
#   include "Convert_Mod/Create_Dual.f90"
#   include "Convert_Mod/Find_Boundary_Faces.f90"
#   include "Convert_Mod/Find_Inside_Faces.f90"
#   include "Convert_Mod/Find_Parents.f90"
#   include "Convert_Mod/Guess_Format.f90"
#   include "Convert_Mod/Insert_Buildings.f90"
#   include "Convert_Mod/Load_Fluent.f90"
#   include "Convert_Mod/Load_Forrest.f90"
#   include "Convert_Mod/Load_Gambit.f90"
#   include "Convert_Mod/Load_Gmsh.f90"
#   include "Convert_Mod/Load_Obj.f90"
#   include "Convert_Mod/Logo_Con.f90"
#   include "Convert_Mod/N_Bnd_Cells_In_Region.f90"
#   include "Convert_Mod/N_Edges_In_Region.f90"
#   include "Convert_Mod/N_Nodes_In_Region.f90"
#   include "Convert_Mod/N_Sharp_Corners.f90"
#   include "Convert_Mod/N_Sharp_Edges.f90"
#   include "Convert_Mod/Save_Fluent.f90"
#   include "Convert_Mod/Sort_Face_Nodes.f90"

  end module
