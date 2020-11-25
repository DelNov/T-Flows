!==============================================================================!
  subroutine Refines_Mod_Allocate_Levels(ref, nl)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Refines_Type) :: ref  ! smoothing regions
  integer            :: nl   ! number of levels
!==============================================================================!

  ref % n_levels  = nl

  allocate(ref % n_regions(nl))
  allocate(ref % region(nl, MAX_SHAPES))

  end subroutine
