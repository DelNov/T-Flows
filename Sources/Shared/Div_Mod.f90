!==============================================================================!
  module Div_Mod
!------------------------------------------------------------------------------!
!   Variable definitions to be used in Divisor.  There are so few in fact      !
!   that there is very little reason for this module to exist at all.          !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Buffer send index and buffer receive index.  
  ! Used for plotting dcomposed grids with links.
  integer, allocatable :: buf_send_ind(:), buf_recv_ind(:), buf_pos(:)

  end module
