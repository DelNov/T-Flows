!==============================================================================!
  module Interface_Mod
!------------------------------------------------------------------------------!
!   Interface between two grids.                                               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
  use Control_Mod
  use Cpu_Timer_Mod
  use Work_Mod,     only: xf_1 => r_face_01,  &  ! face coordinates ...
                          yf_1 => r_face_02,  &  ! ... on the each side ...
                          zf_1 => r_face_03,  &  ! ... of the interface
                          xf_2 => r_face_04,  &
                          yf_2 => r_face_05,  &
                          zf_2 => r_face_06
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !--------------------!
  !   Interface type   !
  !--------------------!
  type Interface_Type

    ! Global number of faces at that interface
    integer :: n_faces

    ! Buffers for storing interface values
    real, allocatable :: phi_1(:)
    real, allocatable :: phi_2(:)

    ! Global cell numbers on each side of the interface
    integer, allocatable :: glo_1(:), glo_2(:)
    integer, allocatable :: bnd_1(:), bnd_2(:)

  end type

  contains

  include 'Interface_Mod/Create.f90'
  include 'Interface_Mod/Exchange.f90'

  end module
