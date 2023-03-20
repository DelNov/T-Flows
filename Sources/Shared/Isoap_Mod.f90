!==============================================================================!
  module Isoap_Mod
!------------------------------------------------------------------------------!
!   Procedures which embed Isoap algorithm                                     !
!                                                                              !
!   I made no attempts to change the code authors have written to comply with  !
!   T-Flows coding standards.  This is in case they release new versions or    !
!   fixes, or what have you.  I don't fully understand how these algorithms    !
!   work, so it is probably better to leave them as they are.                  !
!                                                                              !
!   More details here: https://data.mendeley.com/datasets/4rcf98s74c           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Polyhedron_Mod
  use Iso_Polygons_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Parameters used inside the module
  integer, parameter :: NS=200, NV=240

  !----------------!
  !   Isoap type   !
  !----------------!
  type Isoap_Type

    contains
      procedure :: Main_Isoap
      procedure :: Isopol
#     if T_FLOWS_COMPILATION == 1
      procedure :: Extract_Iso_Polygons
#     endif

  end type

  ! Singleton object
  type(Isoap_Type) :: Isoap

  contains

#   include "Isoap_Mod/Main_Isoap.f90"
#   include "Isoap_Mod/Isopol.f90"
#   if T_FLOWS_COMPILATION == 1
#   include "Isoap_Mod/Extract_Iso_Polygons.f90"
#   endif

  end module
