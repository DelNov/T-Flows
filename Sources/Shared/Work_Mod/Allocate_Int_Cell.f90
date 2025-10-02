!==============================================================================!
  subroutine Allocate_Int_Cell(Work, r)
!------------------------------------------------------------------------------!
!>  Allocates memory for integer-typed working arrays associated with cells
!>  in the Work object.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work     !! parent; the singleton Work object
  integer          :: r        !! number of integer cell arrays
!==============================================================================!

  if(r .eq. 0) return

  allocate(Work % i_cell(r) % array(-Work % max_nb : Work % max_nc))
  Work % i_cell(r) % array(:) = 0

  call Gpu % Vector_Int_Create_On_Device(Work % i_cell(r) % array)

  end subroutine
