!==============================================================================!
  subroutine Refines_Mod_Allocate_Cells(ref, nb, nc)
!------------------------------------------------------------------------------!
!>  Allocates arrays for cell marking and cell level in the Refines_Type
!>  structure (ref), based on the given number of boundary cells (nb) and
!>  cells (nc).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Refines_Type) :: ref     !! refinement regions
  integer            :: nb, nc  ! !number of boundary cells and cells
!==============================================================================!

  allocate (ref % cell_marked(-nb:nc));  ref % cell_marked(:) = .false.
  allocate (ref % cell_level (    nc));  ref % cell_level (:) = 0

  end subroutine
