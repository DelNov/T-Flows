!==============================================================================!
  module Interface_Mod
!------------------------------------------------------------------------------!
!>  Interface_Mod provides a structured way to manage data exchange between
!>  distinct computational grids within a single simulation environment. This
!>  module is particularly useful in multi-domain simulations, such as
!>  conjugate heat transfer scenarios, or cases where one domain serves as
!>  the inflow generator for another.  Since the number of combinations of
!>  variables or values use might want to exchange between two domains is
!>  infinite, this module gives the mechanism to exchange the values between
!>  domains through the so-called buffer, but a specific exchanges is performed
!>  through user function User_Mod_Interface_Exchange.  Check, for an example:
!>  with conjugate heat transfer this file:
!>  [root]/Sources/Process/User_Mod/Interface_Exchange.f90
!>  and for an example of pre-cursor domain:
!>  [root]/Tests/Laminar/Copy_Inlet/User_Mod/Interface_Exchange.f90
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Turb_Mod
  use Control_Mod
  use Profiler_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !--------------------!
  !   Interface type   !
  !--------------------!
  !> The Interface_Type encapsulates data structures and functionalities
  !> necessary for interfacing between two different grids. It includes
  !> pointers to the grids, global face counts, and buffers for
  !> inter-grid value transfer.
  type Interface_Type

    ! Pointers to grids surrounding this interface
    type(Grid_Type), pointer :: pnt_grid1  !! pointer to first grid
    type(Grid_Type), pointer :: pnt_grid2  !! pointer to second grid

    ! Global number of faces at that interface
    integer :: n_tot   !! global number of faces at the interface
    integer :: n1_sub  !! number of faces in first grid in this processor
    integer :: n2_sub  !! number of faces in second grid in this processor

    ! Buffers for storing interface values
    real, allocatable :: phi_1(:,:)  !! buffers in the first grid
    real, allocatable :: phi_2(:,:)  !! buffers in the second grid

    integer, allocatable :: cell_1(:)  !! mapping for cells in first grid
    integer, allocatable :: face_1(:)  !! mapping for faces in first grid
    integer, allocatable :: bcel_1(:)  !! mapping for bnd. cells in first grid
    integer, allocatable :: cell_2(:)  !! mapping for cells in second grid
    integer, allocatable :: face_2(:)  !! mapping for faces in second grid
    integer, allocatable :: bcel_2(:)  !! mapping for bnd. cells in second grid
  end type

  contains

#   include "Interface_Mod/Create.f90"
#   include "Interface_Mod/Exchange.f90"

  end module
