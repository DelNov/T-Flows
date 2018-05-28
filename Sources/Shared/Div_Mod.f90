!==============================================================================!
  module Div_Mod
!------------------------------------------------------------------------------!
!   Variable definitions to be used in Divisor.                                !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Division algorithm
  integer            :: division_algorithm
  integer, parameter :: COORDINATE = 20011
  integer, parameter :: INERTIAL   = 20021
  integer, parameter :: METIS      = 20023

  ! Number of sub-divisions
  integer :: n_sub

  ! Buffer send index and buffer receive index.  
  ! Used for plotting dcomposed grids with links.
  integer, allocatable :: sub_n_cells(:)
  integer, allocatable :: buf_send_ind(:), buf_recv_ind(:), buf_pos(:)

  ! Processor i.d.
  integer, allocatable :: proces(:)  

  ! Axes of innertial
  integer, allocatable :: ix(:), iy(:), iz(:), iin(:)

  ! Division criterion (I believe)
  real, allocatable :: criter(:) 

  end module
