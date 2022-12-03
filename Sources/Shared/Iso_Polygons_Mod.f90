!==============================================================================!
  module Iso_Polygons_Mod
!------------------------------------------------------------------------------!
!   This module is to deal with iso-surface as they are used in Isoap library  !
!   More details here: https://data.mendeley.com/datasets/4rcf98s74c           !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------------!
  !   Iso-polygon type   !
  !----------------------!
  type Iso_Polygons_Type

    integer                              :: n_polys        ! niso
    integer, allocatable, dimension(:)   :: polys_n_verts  ! nipviso (ns)
    integer, allocatable, dimension(:,:) :: polys_v        ! ipviso  (ns,nv)
    real,    allocatable, dimension(:,:) :: verts_xyz      ! vertiso (nv,3)
    integer, allocatable, dimension(:)   :: face_index     ! isoeface(nv)
    integer, allocatable, dimension(:)   :: b_node_1       ! ipia0
    integer, allocatable, dimension(:)   :: b_node_2       ! ipia1
    logical                              :: allocated = .false.

    contains
      procedure :: Allocate_Iso_Polygons
      procedure :: Plot_Iso_Polygons_Vtk

  end type

  ! Singleton type Iso_Polygons object
  type(Iso_Polygons_Type), target :: Iso_Polygons

  contains

#   include "Iso_Polygons_Mod/Allocate_Iso_Polygons.f90"
#   include "Iso_Polygons_Mod/Plot_Iso_Polygons_Vtk.f90"

  end module
