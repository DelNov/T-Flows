!==============================================================================!
  subroutine Allocate_Real_Face(Work, r)
!------------------------------------------------------------------------------!
!>  Allocates memory for real-typed working arrays associated with faces
!>  in the Work object.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work     !! parent; the singleton Work object
  integer          :: r        !! number of real face arrays
!==============================================================================!

  if(r .eq. 0) return

  allocate(Work % r_face(r) % array(1 : Work % max_nf))
  Work % r_face(r) % array(:) = 0.0

  call Gpu % Vector_Real_Create_On_Device(Work % r_face(r) % array)

  end subroutine
