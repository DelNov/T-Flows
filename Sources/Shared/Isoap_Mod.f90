!==============================================================================!
  module Isoap_Mod
!------------------------------------------------------------------------------!
!>  The Isoap_Mod module in T-Flows serves as an interface to the Isoap
!>  library, facilitating the extraction of iso-surfaces from CFD simulations.
!>  (In most cases these surfaces are the interfaces between two phases in VOF
!>  simulations.)  The module defines the Isoap_Type and the singleton object
!>  Isoap for easy access to its functions.
!------------------------------------------------------------------------------!
!   Note: Isoap library is here: https://data.mendeley.com/datasets/4rcf98s74c !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Polyhedron_Mod    ! Polyhedron is sent to Isoap as input
  use Iso_Polygons_Mod  ! Iso_Polygons is the result of a call to Isoap
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Parameters used inside the module
  integer, parameter :: NS=200, NV=240

  !----------------!
  !   Isoap type   !
  !----------------!
  !> A placeholder for two original Isoap functions (Main_Isoap and Isopol)
  !> and Extract_Iso_Polygons, the heart of interfacing between Isoap and
  !> T-Flows.
  type Isoap_Type

    contains
      procedure :: Main_Isoap  !! main Isoap function
      procedure :: Isopol      !! utility from Isoap to arrange iso-vertices
#     if T_FLOWS_COMPILATION == 1
      procedure :: Extract_Iso_Polygons  !! the heart of T-Flows
#     endif                              !! and Isoap interaction

  end type

  ! Singleton object
  type(Isoap_Type) :: Isoap !! singleton object for
              !! easier access to Isoap's functions

  contains

#   include "Isoap_Mod/Main_Isoap.f90"
#   include "Isoap_Mod/Isopol.f90"
#   if T_FLOWS_COMPILATION == 1
#   include "Isoap_Mod/Extract_Iso_Polygons.f90"
#   endif

  end module
