!==============================================================================!
  module Side_Mod
!------------------------------------------------------------------------------!
!   Storage for Side_Type used by Front_Mod and Surf_mod                       !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Side type   !
  !---------------!
  type Side_Type

    integer :: ei, ea, eb  ! element undefined, elements left and right
    integer :: a, b, c, d  ! a and b make sense only for triangular surface
    real    :: length
    logical :: boundary

  end type

  end module
