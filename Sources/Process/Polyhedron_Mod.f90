!==============================================================================!
  module Polyhedron_Mod
!------------------------------------------------------------------------------!
!   This module is to deal with single polyhedron as used in Isoap "library"   !
!   More details here: https://data.mendeley.com/datasets/4rcf98s74c           !
!                                                                              !
!   Conditional compilation allows to test Isoap outside of T-Flows.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
# if T_FLOWS_COMPILATION == 1
    use Grid_Mod
    use Work_Mod
# endif
  use Iso_Polygons_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  integer, parameter :: MAX_ISOAP_FACES = 200
  integer, parameter :: MAX_ISOAP_VERTS = 240

  !---------------------!
  !   Polyhedron type   !
  !---------------------!
  type Polyhedron_Type

    integer                              :: n_nodes        ! ntp
    integer                              :: n_faces        ! nts
    integer, allocatable, dimension(:)   :: faces_n_nodes  ! nipv (ns)
    integer, allocatable, dimension(:,:) :: faces_n        ! ipv  (ns,nv)
    real,    allocatable, dimension(:,:) :: nodes_xyz      ! vertp(nv,3)
    real,    allocatable, dimension(:)   :: phi            ! (nv)
    real                                 :: phiiso

    contains
#     if T_FLOWS_COMPILATION == 1
        procedure        :: Create_From_Grid_Cell  ! for coupling with VOF
#     endif
      procedure          :: Create_Complexcell
      procedure          :: Create_Cube
      procedure          :: Create_Cutcube
      procedure          :: Create_Distortedcube
      procedure          :: Create_Dodecahedron
      procedure          :: Create_Drilledcube
      procedure          :: Create_Hollowedcube
      procedure          :: Create_Icosahedron
      procedure          :: Create_Nchexahedron
      procedure          :: Create_Pentapyramid
      procedure          :: Create_Scube
      procedure          :: Create_Sdodecahedron
      procedure          :: Create_Sicosahedron
      procedure          :: Create_Tetrahedron
      procedure          :: Create_Zigzagcell
      procedure, private :: Func_1
      procedure, private :: Func_2
      procedure, private :: Func_3
      procedure          :: Pick_A_Test_Case
      procedure          :: Plot_Polyhedron_Vtk

  end type

  ! Singleton type Polyhedron and Iso_Polygons objects
  type(Polyhedron_Type)   :: Polyhedron
  type(Iso_Polygons_Type) :: Iso_Polygons

  contains

#   if T_FLOWS_COMPILATION == 1
#   include "Polyhedron_Mod/Create_From_Grid_Cell.f90"
#   endif
#   include "Polyhedron_Mod/Create_Complexcell.f90"
#   include "Polyhedron_Mod/Create_Cube.f90"
#   include "Polyhedron_Mod/Create_Cutcube.f90"
#   include "Polyhedron_Mod/Create_Distortedcube.f90"
#   include "Polyhedron_Mod/Create_Dodecahedron.f90"
#   include "Polyhedron_Mod/Create_Drilledcube.f90"
#   include "Polyhedron_Mod/Create_Hollowedcube.f90"
#   include "Polyhedron_Mod/Create_Icosahedron.f90"
#   include "Polyhedron_Mod/Create_Nchexahedron.f90"
#   include "Polyhedron_Mod/Create_Pentapyramid.f90"
#   include "Polyhedron_Mod/Create_Scube.f90"
#   include "Polyhedron_Mod/Create_Sdodecahedron.f90"
#   include "Polyhedron_Mod/Create_Sicosahedron.f90"
#   include "Polyhedron_Mod/Create_Tetrahedron.f90"
#   include "Polyhedron_Mod/Create_Zigzagcell.f90"
#   include "Polyhedron_Mod/Func_1.f90"
#   include "Polyhedron_Mod/Func_2.f90"
#   include "Polyhedron_Mod/Func_3.f90"
#   include "Polyhedron_Mod/Pick_A_Test_Case.f90"
#   include "Polyhedron_Mod/Plot_Polyhedron_Vtk.f90"

  end module

