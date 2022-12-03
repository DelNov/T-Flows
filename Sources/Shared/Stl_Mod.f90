#include "../Shared/Assert.h90"

!==============================================================================!
  module Stl_Mod
!------------------------------------------------------------------------------!
!   Module to describe geometry in STL format                                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!==============================================================================!

  !--------------!
  !   Stl type   !
  !--------------!
  type, extends(Grid_Type) :: Stl_Type
    real,    allocatable :: nx(:), ny(:), nz(:)  ! facet normal
    integer, allocatable :: body_c(:)            ! body at facets
    integer, allocatable :: body_n(:)            ! body at verticesl

    contains
      procedure, private :: Allocate_Stl
      procedure          :: Create_From_File
      procedure          :: Facet_Coords
      procedure          :: Facet_Normal
      procedure          :: Facets_Vert_Coords
      procedure          :: N_Facets
      procedure, private :: Read_Stl_Ascii
      procedure, private :: Read_Stl_Binary

  end type

  contains
#   include "Stl_Mod/Allocate_Stl.f90"
#   include "Stl_Mod/Create_From_File.f90"
#   include "Stl_Mod/Facet_Coords.f90"
#   include "Stl_Mod/Facet_Normal.f90"
#   include "Stl_Mod/Facets_Vert_Coords.f90"
#   include "Stl_Mod/N_Facets.f90"
#   include "Stl_Mod/Read_Stl_Ascii.f90"
#   include "Stl_Mod/Read_Stl_Binary.f90"

  end module
