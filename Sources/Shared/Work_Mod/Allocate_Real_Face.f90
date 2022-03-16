!==============================================================================!
  subroutine Allocate_Real_Face(Work, Grid, n)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work
  type(Grid_Type)  :: Grid(:)
  integer          :: n     ! number of real cell arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nf, i
!==============================================================================!

  if(n .eq. 0) return

  ! Get number of cells and boundary cells
  nf = maxval(Grid(1:size(Grid)) % n_faces)

  allocate(Work % r_face(n))

  do i = 1, n
    allocate(Work % r_face(i) % ptr(1:nf))
    Work % r_face(i) % ptr(:) = 0.0
  end do

  Work % last_r_face = 0

  end subroutine
