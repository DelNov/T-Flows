!==============================================================================!
  subroutine Refines_Mod_Allocate_Cells(ref, nb, nc)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod, only: HUGE
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Refines_Type) :: ref     ! smoothing regions
  integer            :: nb, nc  ! number of boundary cells and cells
!==============================================================================!

  allocate (ref % cell_marked(-nb:nc));  ref % cell_marked(:) = .false.
  allocate (ref % cell_level (    nc));  ref % cell_level (:) = 0

  end subroutine
