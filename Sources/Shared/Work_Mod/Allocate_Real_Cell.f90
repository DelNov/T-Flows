!==============================================================================!
  subroutine Allocate_Real_Cell(Work, r)
!------------------------------------------------------------------------------!
!>  Allocates memory for real-typed working arrays associated with cells
!>  in the Work object.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work     !! parent; the singleton Work object
  integer          :: r        !! number of real cell arrays
!==============================================================================!

  if(r .eq. 0) return

  allocate(Work % r_cell(r) % array(-Work % max_nb : Work % max_nc))
  Work % r_cell(r) % array(:) = 0.0

  call Gpu % Vector_Real_Create_On_Device(Work % r_cell(r) % array)

  end subroutine
