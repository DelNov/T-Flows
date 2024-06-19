!==============================================================================!
  subroutine Allocate_Real_Face(Work, Grid, n)
!------------------------------------------------------------------------------!
!>  Allocates memory for real-typed working arrays associated with faces
!>  in the Work object.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work     !! parent; the singleton Work object
  type(Grid_Type)  :: Grid(:)  !! grids on which the Work will be used
  integer          :: n        !! number of real face arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nf, i
!==============================================================================!

  if(n .eq. 0) return

  ! Get number of faces
  nf = maxval(Grid(1:size(Grid)) % n_faces)

  allocate(Work % r_face(n))

  do i = 1, n
    allocate(Work % r_face(i) % array(1:nf))
    Work % r_face(i) % array(:) = 0.0
  end do

  Work % last_r_face = 0

  end subroutine
