!==============================================================================!
  subroutine Allocate_Int_Face(Work, Grid, n)
!------------------------------------------------------------------------------!
!>  Allocates memory for integer-typed working arrays associated with faces
!>  in the Work object.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work     !! parent; the singleton Work object
  type(Grid_Type)  :: Grid(:)  !! grids on which the Work will be used
  integer          :: n        !! number of integer face arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nf, i
!==============================================================================!

  if(n .eq. 0) return

  ! Get number of faces
  nf = maxval(Grid(1:size(Grid)) % n_faces)

  allocate(Work % i_face(n))

  do i = 1, n
    allocate(Work % i_face(i) % array(1:nf))
    Work % i_face(i) % array(:) = 0
  end do

  Work % last_i_face = 0

  end subroutine
