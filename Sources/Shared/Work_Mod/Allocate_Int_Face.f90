!==============================================================================!
  subroutine Allocate_Int_Face(Work, r)
!------------------------------------------------------------------------------!
!>  Allocates memory for integer-typed working arrays associated with faces
!>  in the Work object.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work     !! parent; the singleton Work object
  integer          :: r        !! number of integer face arrays
!==============================================================================!

  if(r .eq. 0) return

  allocate(Work % i_face(r) % array(1 : Work % max_nf))
  Work % i_face(r) % array(:) = 0

  call Gpu % Vector_Int_Create_On_Device(Work % i_face(r) % array)

  end subroutine
