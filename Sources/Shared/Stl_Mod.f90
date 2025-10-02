#include "../Shared/Assert.h90"

!==============================================================================!
  module Stl_Mod
!------------------------------------------------------------------------------!
!>  Stl_Mod is a specialized module in T-Flows for handling Stereolithography
!>  (STL) file formats. It extends Grid_Mod to include as much functionalities
!>  This module is crucial for defining porous regions and interfaces in VOF
!>  simulations, or any other processing which starts from an STL file.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!==============================================================================!

  !--------------!
  !   Stl type   !
  !--------------!
  !> Encapsulates data and functionality to manipulate and process STL files.
  type, extends(Grid_Type) :: Stl_Type
    real,    allocatable, private :: nx(:)      !! facet normal's x component
    real,    allocatable, private :: ny(:)      !! facet normal's y component
    real,    allocatable, private :: nz(:)      !! facet normal's z component
    integer                       :: n_boddies  !! number of independent boddies
    integer, allocatable          :: body_c(:)  !! body at facets
    integer, allocatable          :: body_n(:)  !! body at vertices

    contains
      procedure, private :: Allocate_Stl
      procedure          :: Create_From_File
      procedure          :: Facet_Body
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
#   include "Stl_Mod/Facet_Body.f90"
#   include "Stl_Mod/Facet_Coords.f90"
#   include "Stl_Mod/Facet_Normal.f90"
#   include "Stl_Mod/Facets_Vert_Coords.f90"
#   include "Stl_Mod/N_Facets.f90"
#   include "Stl_Mod/Read_Stl_Ascii.f90"
#   include "Stl_Mod/Read_Stl_Binary.f90"

  end module
