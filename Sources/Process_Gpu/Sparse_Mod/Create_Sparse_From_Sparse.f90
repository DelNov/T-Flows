!==============================================================================!
  subroutine Create_Sparse_From_Sparse(A, B)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sparse_Type)        :: A
  type(Sparse_Type), target :: B
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: i
!==============================================================================!

  ! Fetch the alias
  Grid => B % pnt_grid

  ! Store pointer to the grid
  A % pnt_grid => Grid

  print '(a)', ' # Creating a sparse matrix from another matrix'

  ! Copy the number of nonzeros
  A % nonzeros = B % nonzeros

  ! Allocatte the memory for all significant fields from the matrix ...

  allocate (A % val(A % nonzeros))
  allocate (A % fc (Grid % n_faces))
  allocate (A % col(A % nonzeros))
  allocate (A % row(Grid % n_cells+1))
  allocate (A % dia(Grid % n_cells))
  allocate (A % pos(2, Grid % n_faces))
  allocate (A % d_inv(Grid % n_cells))
  allocate (A % v_m  (Grid % n_cells))

  ! ... and then copy them from one matrix to another

# ifdef __INTEL_COMPILER
  do i = 1, A % nonzeros
    A % val(i) = B % val(i)
    A % col(i) = B % col(i)
  end do
  do i = 1, Grid % n_faces
    A % fc(i)    = B % fc(i)
    A % pos(1,i) = B % pos(1,i)
    A % pos(2,i) = B % pos(2,i)
  end do
  do i = 1, Grid % n_cells
    A % row  (i) = B % row  (i)
    A % dia  (i) = B % dia  (i)
    A % d_inv(i) = B % d_inv(i)
    A % v_m  (i) = B % v_m  (i)
  end do
  Assert(i .eq. Grid % n_cells + 1)
  A % row(i) = B % row(i)
# else
  A % val(:)   = B % val(:)
  A % fc(:)    = B % fc(:)
  A % col(:)   = B % col(:)
  A % row(:)   = B % row(:)
  A % dia(:)   = B % dia(:)
  A % pos(:,:) = B % pos(:,:)
  A % d_inv(:) = B % d_inv(:)
  A % v_m  (:) = B % v_m  (:)
# endif

  end subroutine
