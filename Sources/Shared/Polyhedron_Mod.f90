#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Polyhedron_Mod
!------------------------------------------------------------------------------!
!>  This module integrates with the Isoap library, focusing on single polyhedron
!>  handling. It defines Polyhedron_Type with data members and procedures for
!>  easier interaction with Isoap's polyhedron structures. This module also
!>  supports conditional compilation for testing outside of T-Flows, enhancing
!>  its versatility and ease of testing and verification.
!------------------------------------------------------------------------------!
!   Note 1: Isoap library is here https://data.mendeley.com/datasets/4rcf98s74c!
!   Note 2: Conditional compilation allows to test Isoap outside of T-Flows.   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
# if T_FLOWS_COMPILATION == 1
  use Work_Mod
# endif
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  integer, parameter :: MAX_ISOAP_FACES = 200
  integer, parameter :: MAX_ISOAP_VERTS = 240

  !---------------------!
  !   Polyhedron type   !
  !---------------------!
  type Polyhedron_Type

    integer              :: n_nodes           !! ntp          in Isoap
    integer              :: n_faces           !! nts          in Isoap
    integer, allocatable :: faces_n_nodes(:)  !! nipv (ns)    in Isoap
    integer, allocatable :: faces_n    (:,:)  !! ipv  (ns,nv) in Isoap
    real,    allocatable :: nodes_xyz  (:,:)  !! vertp(nv,3)  in Isoap
    real,    allocatable :: phi          (:)  !! (nv)
    real                 :: phi_iso
    integer, allocatable :: global_node  (:)
    logical              :: allocated = .false.

    contains
      procedure          :: Allocate_Polyhedron
#   if T_FLOWS_COMPILATION == 1
      procedure, private :: Calculate_Cell_Centroid
      procedure          :: Calculate_Cell_Volume
      procedure, private :: Calculate_Face_Centroid
      procedure          :: Create_From_Polyhedron
#   else
      procedure          :: Create_Complexcell    ! exclude from FORD
      procedure          :: Create_Cube           ! exclude from FORD
      procedure          :: Create_Cutcube        ! exclude from FORD
      procedure          :: Create_Distortedcube  ! exclude from FORD
      procedure          :: Create_Dodecahedron   ! exclude from FORD
      procedure          :: Create_Drilledcube    ! exclude from FORD
      procedure          :: Create_Hollowedcube   ! exclude from FORD
      procedure          :: Create_Icosahedron    ! exclude from FORD
      procedure          :: Create_Nchexahedron   ! exclude from FORD
      procedure          :: Create_Pentapyramid   ! exclude from FORD
      procedure          :: Create_Scube          ! exclude from FORD
      procedure          :: Create_Sdodecahedron  ! exclude from FORD
      procedure          :: Create_Sicosahedron   ! exclude from FORD
      procedure          :: Create_Tetrahedron    ! exclude from FORD
      procedure          :: Create_Zigzagcell     ! exclude from FORD
#   endif
#   if T_FLOWS_COMPILATION == 1
      procedure          :: Extract_From_Grid
#   else
      procedure, private :: Func_1                ! exclude from FORD
      procedure, private :: Func_2                ! exclude from FORD
      procedure, private :: Func_3                ! exclude from FORD
      procedure          :: Pick_A_Test_Case      ! exclude from FORD
#   endif
      procedure          :: Plot_Polyhedron_Vtk

  end type

  ! Singleton type Polyhedron object
  type(Polyhedron_Type), target :: Polyhedron !! singleton object for
              !! easier access to inerfaces between T-Flows and Isoap

  contains

#   include "Polyhedron_Mod/Allocate_Polyhedron.f90"
# if T_FLOWS_COMPILATION == 1
#   include "Polyhedron_Mod/Calculate_Cell_Centroid.f90"
#   include "Polyhedron_Mod/Calculate_Cell_Volume.f90"
#   include "Polyhedron_Mod/Calculate_Face_Centroid.f90"
#   include "Polyhedron_Mod/Create_From_Polyhedron.f90"
# else
#   include "Polyhedron_Mod/Create_Complexcell.f90"    ! exclude from FORD
#   include "Polyhedron_Mod/Create_Cube.f90"           ! exclude from FORD
#   include "Polyhedron_Mod/Create_Cutcube.f90"        ! exclude from FORD
#   include "Polyhedron_Mod/Create_Distortedcube.f90"  ! exclude from FORD
#   include "Polyhedron_Mod/Create_Dodecahedron.f90"   ! exclude from FORD
#   include "Polyhedron_Mod/Create_Drilledcube.f90"    ! exclude from FORD
#   include "Polyhedron_Mod/Create_Hollowedcube.f90"   ! exclude from FORD
#   include "Polyhedron_Mod/Create_Icosahedron.f90"    ! exclude from FORD
#   include "Polyhedron_Mod/Create_Nchexahedron.f90"   ! exclude from FORD
#   include "Polyhedron_Mod/Create_Pentapyramid.f90"   ! exclude from FORD
#   include "Polyhedron_Mod/Create_Scube.f90"          ! exclude from FORD
#   include "Polyhedron_Mod/Create_Sdodecahedron.f90"  ! exclude from FORD
#   include "Polyhedron_Mod/Create_Sicosahedron.f90"   ! exclude from FORD
#   include "Polyhedron_Mod/Create_Tetrahedron.f90"    ! exclude from FORD
#   include "Polyhedron_Mod/Create_Zigzagcell.f90"     ! exclude from FORD
# endif
# if T_FLOWS_COMPILATION == 1
#   include "Polyhedron_Mod/Extract_From_Grid.f90"
# else
#   include "Polyhedron_Mod/Func_1.f90"                ! exclude from FORD
#   include "Polyhedron_Mod/Func_2.f90"                ! exclude from FORD
#   include "Polyhedron_Mod/Func_3.f90"                ! exclude from FORD
#   include "Polyhedron_Mod/Pick_A_Test_Case.f90"      ! exclude from FORD
# endif
#   include "Polyhedron_Mod/Plot_Polyhedron_Vtk.f90"

  end module

