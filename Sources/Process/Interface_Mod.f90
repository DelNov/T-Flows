!==============================================================================!
  module Interface_Mod
!------------------------------------------------------------------------------!
!   Interface between two grids.                                               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
  use Control_Mod
  use Cpu_Timer_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !--------------------!
  !   Interface type   !
  !--------------------!
  type Interface_Type

    integer              :: n_faces        ! number of faces at that interface
    integer, allocatable :: faces_1   (:)  ! list of faces in domain 1
    integer, allocatable :: close_in_2(:)  ! list of closest faces in domain 2
    integer, allocatable :: faces_2   (:)  ! list of faces in domain 2
    integer, allocatable :: close_in_1(:)  ! list of closest faces in domain 1

  end type

  contains

  include 'Interface_Mod/Create.f90'
  include 'Interface_Mod/Exchange.f90'

  end module
