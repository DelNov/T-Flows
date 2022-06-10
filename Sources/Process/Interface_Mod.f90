!==============================================================================!
  module Interface_Mod
!------------------------------------------------------------------------------!
!   Interface between two grids.                                               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Turb_Mod
  use Control_Mod
  use Cpu_Timer_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !--------------------!
  !   Interface type   !
  !--------------------!
  type Interface_Type

    ! Pointers to grids surrounding this interface
    type(Grid_Type), pointer :: pnt_grid1
    type(Grid_Type), pointer :: pnt_grid2

    ! Global number of faces at that interface
    integer :: n_tot
    integer :: n1_sub
    integer :: n2_sub

    ! Buffers for storing interface values
    real, allocatable :: phi_1(:,:)
    real, allocatable :: phi_2(:,:)

    integer, allocatable :: cell_1(:)
    integer, allocatable :: face_1(:)
    integer, allocatable :: bcel_1(:)
    integer, allocatable :: cell_2(:)
    integer, allocatable :: face_2(:)
    integer, allocatable :: bcel_2(:)
  end type

  contains

  include 'Interface_Mod/Create.f90'
  include 'Interface_Mod/To_Buffer.f90'

  end module
