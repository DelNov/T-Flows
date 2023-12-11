!==============================================================================!
  subroutine Refines_Mod_Allocate_Levels(ref, nl)
!------------------------------------------------------------------------------!
!>  Initializes the number of refinement levels in the Refines_Type
!>  structure (ref), and allocates arrays for storing range information
!>  associated with each level, based on the provided number of levels (nl).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Refines_Type)  :: ref  !! refinement region
  integer, intent(in) :: nl   !! number of levels
!==============================================================================!

  ref % n_levels  = nl

  allocate(ref % n_ranges(nl))
  allocate(ref % range(nl, MAX_SHAPES))

  end subroutine
