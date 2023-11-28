!==============================================================================!
  module Iso_Polygons_Mod
!------------------------------------------------------------------------------!
!>  This module is designed to facilitate interaction between Isoap library
!>  and T-Flows.  It defines the Iso_Polygons_Type which contains all the data
!>  members an Isoap's iso-polygon would have, but with less cryptic names.
!>  Furthermore, the Iso_Polygons_Type defines two member functions which
!>  allocate memory for iso-polygons, and plot it in .vtk format for visual
!>  inspection of created polygons.
!------------------------------------------------------------------------------!
!   Isoap library is here: https://data.mendeley.com/datasets/4rcf98s74c       !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------------!
  !   Iso-polygon type   !
  !----------------------!
  !> Encapsulates data members to facilitate interaction with Isoap library.
  !> All data members correspond to variables in the Isoap library, but with
  !> more descriptive and user-friendly names.  The only new data member is
  !> the flag 'allocated', set to true when object is allocated.
  type Iso_Polygons_Type

    integer              :: n_polys              !! niso            in Isaop
    integer, allocatable :: polys_n_verts(:)     !! nipviso (ns)    in Isoap
    integer, allocatable :: polys_v    (:,:)     !! ipviso  (ns,nv) in Isoap
    real,    allocatable :: verts_xyz  (:,:)     !! vertiso (nv,3)  in Isoap
    integer, allocatable :: face_index   (:)     !! isoeface(nv)    in Isoap
    integer, allocatable :: b_node_1     (:)     !! ipia0           in Isoap
    integer, allocatable :: b_node_2     (:)     !! ipia1           in Isoap
    logical              :: allocated = .false.  !! true if allocated

    contains
      procedure :: Allocate_Iso_Polygons
      procedure :: Plot_Iso_Polygons_Vtk

  end type

  ! Singleton type Iso_Polygons object
  type(Iso_Polygons_Type), target :: Iso_Polygons  !! singleton object for
                   !! easier access to inerfaces between T-Flows and Isoap

  contains

#   include "Iso_Polygons_Mod/Allocate_Iso_Polygons.f90"
#   include "Iso_Polygons_Mod/Plot_Iso_Polygons_Vtk.f90"

  end module
