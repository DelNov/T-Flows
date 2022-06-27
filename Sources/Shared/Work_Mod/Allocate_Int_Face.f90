!==============================================================================!
  subroutine Allocate_Int_Face(Work, Grid, n)
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

  allocate(Work % i_face(n))

  do i = 1, n
    allocate(Work % i_face(i) % ptr(1:nf))
    Work % i_face(i) % ptr(:) = 0
  end do

  Work % last_i_face = 0

  end subroutine
