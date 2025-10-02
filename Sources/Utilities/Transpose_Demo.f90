!==============================================================================!
!
! Imagine this matrix:
!
!   a   = |  40 -15   .  -5 |
!         |  -2  60  -1   . |
!         |   .  -3  50  -7 |
!         |   .  -8  -6  70 |
!
!   n  =  4
!   nz = 12
!
!   a_val    (1:nz)  = | 40 -15  -5  60  -2  -1  50  -3  -7  70  -8  -6 |
!
!   col_idx  (1:nz)  = |  1   2   4   2   1   3   3   2   4   4   2   3 |
!
!   row_ptr  (1:n+1) = |  1           4           7          10          13  |
!
! The transpose of that matrix reads
!
!   a^t = |  40  -2   .   . |
!         | -15  60  -3  -8 |
!         |   .  -1  50  -6 |
!         |  -5   .  -7  70 |
!
!   a^t_val  (1:nz)  = | 40  -2  60 -15  -3  -8  50  -1  -6  70  -5  -7 |
!
!   col^t_idx(1:nz)  = |  1   2   2   1   3   4   3   2   4   4   1   3 |
!
!   row^t_ptr(1:n+1) = |  1       3               7          10          13  |
!
!------------------------------------------------------------------------------!
  program Transpose_Demo

  implicit none

  integer :: i, j, k, n, nnz

  !---------------------------------------------------------------------!
  !                                                                     !
  !   Define a local CRS matrix and allocate memory for its transpose   !
  !                                                                     !
  !---------------------------------------------------------------------!

  ! Original system (unused)
  ! integer :: a(4, 4) = reshape([ &
  !   40,-15,  0, -5,  &  ! row 1
  !   -2, 60, -1,  0,  &  ! row 2
  !    0, -3, 50, -7,  &  ! row 3
  !    0, -8, -6, 70   &  ! row 4
  ! ], [4, 4])
  integer :: a_val  (12) = [40,-15, -5, 60, -2, -1, 50, -3, -7, 70, -8, -6]
  integer :: col_idx(12) = [ 1,  2,  4,  2,  1,  3,  3,  2,  4,  4,  2,  3]
  integer :: row_ptr (5) = [ 1,          4,          7,         10,         13]

  ! Transposed system (unused)
  ! integer :: a_t(4, 4) = reshape([ &
  !    0,  0,  0,  0,  &  ! row 1
  !    0,  0,  0,  0,  &  ! row 2
  !    0,  0,  0,  0,  &  ! row 3
  !    0,  0,  0,  0   &  ! row 4
  ! ], [4, 4])
  integer :: a_t_val  (12) = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  integer :: col_t_idx(12) = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  integer :: row_t_ptr (5) = [ 0, 0, 0, 0, 0]

  ! Helping arrays
  integer, allocatable          :: counter(:)
  double precision, allocatable :: x(:), y(:), z(:)

  !----------------------------------------------!
  !                                              !
  !   Create transpose of the local CRS matrix   !
  !                                              !
  !----------------------------------------------!

  !-------------------------------!
  !   Phase 0 - Initializations   !
  !-------------------------------!

  ! Fetch the number of unknowns and non-zeroes
  n   = size(row_ptr) - 1
  nnz = size(a_val)

  ! Allocate memory for helping arrays
  allocate(counter(n))
  print *, "Size of original matrix: ", n
  print *, "Number of non-zeroes   : ", nnz

  !----------------------------------------------------!
  !   Phase 1 - count the columns in original matrix   !
  !----------------------------------------------------!
  counter(:) = 0  ! acts as column counter

  do i = 1, n                      ! rows in the original matrix
    do k = row_ptr(i), row_ptr(i+1) - 1
      j = col_idx(k)               ! columns in the original matrix
      counter(j) = counter(j) + 1  ! increase the column count
    end do
  end do

  print *, "Size of each column:"
  print '(64i4)', counter(1:n)

  !---------------------------------------------------!
  !   Phase 2 - form the row_t_ptr based on counter   !
  !---------------------------------------------------!
  row_t_ptr(1) = 1  ! must start at 1
  do i = 2, n + 1   ! ends at n+1 by convention
    row_t_ptr(i) = row_t_ptr(i-1) + counter(i-1)
  end do

  print *, "Row pointer for transposed matrix:"
  print '(64i4)', row_t_ptr(1:n+1)

  !-------------------------------------------!
  !   Phase 3 - fill the rest of the matrix   !
  !-------------------------------------------!
  counter(:) = 0      ! re-initialize the column/row counter

  ! Phase 3.1 - fill the diagonals up
  do i = 1, n       ! rows in the original matrix
    k = row_ptr(i)
    j = col_idx(k)  ! column in the original matrix
                    ! but row in the transposed
    a_t_val  (row_t_ptr(j)) = a_val(k)
    col_t_idx(row_t_ptr(j)) = i
    counter(j) = 1  ! first entry is the diagonal
  end do

  ! Phase 3.2 - populate the off-diagonal terms
  do i = 1, n         ! rows in the original matrix
    do k = row_ptr(i) + 1, row_ptr(i+1) - 1
      j = col_idx(k)  ! columns in the original matrix
                      ! but rows in the transposed
      a_t_val  (row_t_ptr(j) + counter(j)) = a_val(k)
      col_t_idx(row_t_ptr(j) + counter(j)) = i
      counter(j) = counter(j) + 1  ! update the counter
    end do
  end do

  print '(64i4)', a_t_val  (1:nnz)
  print '(64i4)', col_t_idx(1:nnz)

  !----------------------------!
  !                            !
  !   Validate the transpose   !
  !                            !
  !----------------------------!

  allocate(x(n), y(n), z(n))

  ! Fill x with something non-trivial
  do i = 1, n
    x(i) = dble(i) * 0.1d0
  end do

  ! Compute y = A*x
  do i = 1, n
    do k = row_ptr(i), row_ptr(i+1)-1
      j = col_idx(k)
      y(i) = y(i) + a_val(k) * x(j)
    end do
  end do

  ! Compute z = A^T*x
  do i = 1, n
    do k = row_t_ptr(i), row_t_ptr(i+1)-1
      j = col_t_idx(k)
      z(i) = z(i) + a_t_val(k) * x(j)
    end do
  end do

  ! Now compare inner products
  print *, "dot((A^T x), x) =", dot_product(y, x)
  print *, "dot(x, A x)     =", dot_product(x, z)

  ! Optionally print difference
  print *, "Difference      =", abs(dot_product(y, x) - dot_product(x, z))

  end program
