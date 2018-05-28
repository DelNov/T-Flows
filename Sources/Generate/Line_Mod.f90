!==============================================================================!
  module Line_Mod
!------------------------------------------------------------------------------!
!   Lines defining blocks used in "Generator"                                  !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Line type   !
  !---------------!
  type Line_Type

    integer           :: points(2)
    real, allocatable :: x(:)
    real, allocatable :: y(:)
    real, allocatable :: z(:)
    real              :: weight      ! line weights for node clustering
    integer           :: resolution

  end type

  end module
