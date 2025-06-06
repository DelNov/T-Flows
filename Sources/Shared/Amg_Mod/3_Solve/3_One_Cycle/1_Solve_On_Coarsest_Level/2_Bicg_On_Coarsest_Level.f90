! #define DEBUG  1

!==============================================================================!
  subroutine bicg_on_coarsest_level(amg, level, irel,  &
                                    a, u, f, ia, ja,   &  ! defining system
                                    iw, icg)
!------------------------------------------------------------------------------!
!   Bi-conjugate gradient on the coarsest level level
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: level, irel
  double precision :: a(:), u(:), f(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:), icg(:)
!-----------------------------------[locals]-----------------------------------!
  double precision :: s
  integer          :: ig, iaux, jg
!---------------------------------[new locals]---------------------------------!
  integer                           :: i, iter, j, k, n, nnz
  double precision,     allocatable :: m(:), b(:), x(:)
  double precision,     allocatable :: p(:),   q(:),   r(:),   z(:)
  double precision,     allocatable :: p_t(:), q_t(:), r_t(:), z_t(:)
  double precision,     allocatable :: a_val(:)
  integer, allocatable              :: row_ptr(:), col_idx(:)
  double precision,     allocatable :: a_t_val(:)
  integer, allocatable              :: row_t_ptr(:), col_t_idx(:)
  integer, allocatable              :: counter(:)
  double precision                  :: alpha, beta
  double precision                  :: pq, rho_0, rho_old, rho_new
!==============================================================================!

  call amg % timer_start()

  ! See comment in source "coarsening.f90" at line 180
  iaux = ia(amg % imax(level)+1)
  ia(amg % imax(level)+1) = iw(amg % iminw(level))


  !---------------------------------------------!
  !                                             !
  !                                             !
  !   Create local linear system of equations   !
  !                                             !
  !                                             !
  !---------------------------------------------!

  ! Number of unknowns and non-zeros on this level
  n = amg % imax(level) - amg % imin(level) + 1
  nnz = ia(amg % imax(level)+1) - ia(amg % imin(level))

  !-------------------------------------------------------!
  !                                                       !
  !   Create local CRS matrix from the global ia, ja, a   !
  !                                                       !
  !-------------------------------------------------------!
  allocate(a_val(nnz));   a_val(:)   = 0.0
  allocate(col_idx(nnz)); col_idx(:) = 0
  allocate(row_ptr(n+1)); row_ptr(:) = 0
  do ig = amg % imin(level), amg % imax(level)
    i = ig - amg % imin(level) + 1                   ! local unknown number
    row_ptr(i) = ia(ig) - ia(amg % imin(level)) + 1  ! local row pointer
    do jg = ia(ig), ia(ig+1) - 1                     ! browse through row ig
      j          = jg - ia(amg % imin(level)) + 1    ! local nonzero index
      col_idx(j) = ja(jg) - amg % imin(level) + 1    ! local column number
      a_val(j)   = a(jg)                             ! local matrix value
    end do
  end do
  ! Final row pointer
  row_ptr(n+1) = ia(amg % imax(level)+1) - ia(amg % imin(level)) + 1

  !------------------------------------------------------!
  !                                                      !
  !   Create local b and x vectors (Stueben's f and u)   !
  !                                                      !
  !------------------------------------------------------!
  allocate(b(n))
  allocate(x(n))
  do ig = amg % imin(level), amg % imax(level)
    i = ig - amg % imin(level) + 1  ! local unkown number
    b(i) = f(ig)
    x(i) = u(ig)
  end do

  !---------------------------------------------------------!
  !                                                         !
  !   Allocate memory for local vectors for the algorithm   !
  !                                                         !
  !---------------------------------------------------------!
  allocate(m(n))
  allocate(p(n),   q(n),   r(n),   z(n))
  allocate(p_t(n), q_t(n), r_t(n), z_t(n))
  m(:) = 0.0
  p(:) = 0.0;  p_t(:) = 0.0
  q(:) = 0.0;  q_t(:) = 0.0
  r(:) = 0.0;  r_t(:) = 0.0
  z(:) = 0.0;  z_t(:) = 0.0

  !----------------------------------------------!
  !                                              !
  !   Create transpose of the local CRS matrix   !
  !                                              !
  !----------------------------------------------!
  allocate(a_t_val(nnz));   a_t_val(:)   = 0.0
  allocate(col_t_idx(nnz)); col_t_idx(:) = 0
  allocate(row_t_ptr(n+1)); row_t_ptr(:) = 0

  !----------------------------------------------------!
  !   Phase 1 - count the columns in original matrix   !
  !----------------------------------------------------!
  allocate(counter(nnz))
  counter(:) = 0  ! acts as column counter

  do i = 1, n                      ! rows in the original matrix
    do k = row_ptr(i), row_ptr(i+1) - 1
      j = col_idx(k)               ! columns in the original matrix
      counter(j) = counter(j) + 1  ! increase the column count
    end do
  end do

# ifdef DEBUG
    print *, "Size of each column:"
    print '(64i4)', counter(1:n)
# endif

  !---------------------------------------------------!
  !   Phase 2 - form the row_t_ptr based on counter   !
  !---------------------------------------------------!
  row_t_ptr(1) = 1  ! must start at 1
  do i = 2, n + 1   ! ends at n+1 by convention
    row_t_ptr(i) = row_t_ptr(i-1) + counter(i-1)
  end do

# ifdef DEBUG
    print *, "Row pointer for transposed matrix:"
    print '(64i4)', row_t_ptr(1:n+1)
# endif

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

!#  !--------------------------------------!
!#  !   Phase 4 - validate the transpose   !
!#  !--------------------------------------!
!#
!#  ! Fill p with something non-trivial
!#  do i = 1, n
!#    p(i) = dble(i) * 0.1d0
!#  end do
!#
!#  ! Compute q = A*p
!#  do i = 1, n
!#    do k = row_ptr(i), row_ptr(i+1)-1
!#      j = col_idx(k)
!#      q(i) = q(i) + a_val(k) * p(j)
!#    end do
!#  end do
!#
!#  ! Compute r = A^T*p
!#  do i = 1, n
!#    do k = row_t_ptr(i), row_t_ptr(i+1)-1
!#      j = col_t_idx(k)
!#      r(i) = r(i) + a_t_val(k) * p(j)
!#    end do
!#  end do
!#
!#  ! Now compare inner products
!#  print *, "dot((A^T p), p) =", dot_product(q, p)
!#  print *, "dot(p, A p)     =", dot_product(p, r)
!#
!#  ! Optionally print difference
!#  print *, "Difference      =", abs(dot_product(q, p) - dot_product(p, r))


  !---------------------------------!
  !   Form preconditioning matrix   !   (it is the same for a and a_t)
  !---------------------------------!
  do i = 1, n
    m(i) = 1.0d0 / a_val(row_ptr(i))  ! diagonal is at a_val(row_ptr(i))
  end do

  !--------------------------------------------------------------!
  !   Compute r^(0) = b - A x^(0) for some initial guess x^(0)   !
  !           _                  _                               !
  !   Choose  r^(0) (for example r^(0) = r^(0))                  !
  !--------------------------------------------------------------!
  rho_0 = 0.0
  do i = 1, n
    s = b(i)
    do k = row_ptr(i), row_ptr(i+1) - 1
      j = col_idx(k)
      s = s - a_val(k) * x(j)
    end do
    r  (i) = s
    r_t(i) = r(i)
    rho_0 = rho_0 + r(i) * r(i)
  end do

# ifdef DEBUG
    print *, "r_ini = ", sqrt(rho_0)
# endif

  ! Set the value which will be used in the loop
  ! (Probably not needed)
  rho_old = rho_0

  !--------------------!
  !                    !
  !   Iteration loop   !
  !                    !
  !--------------------!
  do iter = 1, 60

    !---------------------------------!
    !   Solve M   z^(i-1) = r^(i-1)   !
    !             _         _         !
    !   Solve M^T z^(i-1) = r^(i-1)   !
    !---------------------------------!
    do i = 1, n
      z  (i) = m(i) * r  (i)
      z_t(i) = m(i) * r_t(i)
    end do

    !---------------_-------------------!
    !   rho_(i-1) = r^(i-1)^T z^(i-1)   !
    !-----------------------------------!
    rho_new = 0.0d0
    do i = 1, n
      rho_new = rho_new + r_t(i) * z(i)
    end do

    if(iter .eq. 1) then
      !-------------------!
      !   p^(1) = z^(0)   !
      !   _       _       !
      !   p^(1) = z^(0)   !
      !-------------------!
      do i = 1, n
        p  (i) = z  (i)
        p_t(i) = z_t(i)
      end do
    else
      !--------------------------------------------!
      !   beta_(i-1) = rho_(i-1) / rho_(i-2)       !
      !                                            !
      !   p^(i) = z^(i-1) + beta_(i-1) * p^(i-1)   !
      !   _       _                      _         !
      !   p^(i) = z^(i-1) + beta_(i-1) * p^(i-1)   !
      !--------------------------------------------!
      beta = rho_new / rho_old
      do i = 1, n
        p  (i) = z  (i) + beta * p  (i)
        p_t(i) = z_t(i) + beta * p_t(i)
      end do
    end if

    !-----------------------!
    !   q^(i) = A * p^(i)   !
    !-----------------------!
    do i = 1, n
      s = 0.0d0
      do k = row_ptr(i), row_ptr(i+1) - 1
        j = col_idx(k)
        s = s + a_val(k) * p(j)
      end do
      q(i) = s
    end do

    !---_-------------_-------!
    !   q^(i) = A^T * p^(i)   !
    !-------------------------!
    do i = 1, n
      s = 0.0d0
      do k = row_t_ptr(i), row_t_ptr(i+1) - 1
        j = col_t_idx(k)
        s = s + a_t_val(k) * p_t(j)
      end do
      q_t(i) = s
    end do

    !---------------------_---------------!
    !   alpha = rho_(i-1)/p^(i)^T q^(i)   !
    !-------------------------------------!
    pq = 0.0d0
    do i = 1, n
      pq = pq + p_t(i) * q(i)
    end do
    alpha = rho_new / pq

    !-------------------------------------!
    !   x^(i) = x^(i-1) + alpha * p^(i)   !
    !                                     !
    !   r^(i) = r^(i-1) - alpha * q^(i)   !
    !   _       _                 _       !
    !   r^(i) = r^(i-1) - alpha * q^(i)   !
    !-------------------------------------!
    do i = 1, n
      x  (i) = x  (i) + alpha * p  (i)
      r  (i) = r  (i) - alpha * q  (i)
      r_t(i) = r_t(i) - alpha * q_t(i)
    end do

# ifdef DEBUG
    print '(a,i3,a,1pe14.7)',  &
      "iter: ", iter, " r_new/r_ini = ", sqrt(rho_new/rho_0)
# endif
    if(sqrt(rho_new) .lt. 1.0e-12) exit

    rho_old = rho_new

  end do

  ! Print final residual
# ifdef DEBUG
    print '(a,i3,a,1pe14.7)',  &
      "iter: ", iter, " r_new = ", sqrt(rho_new)
# endif

  !---------------------------------------------------------!
  !                                                         !
  !   Copy local b and x vectors back (Stueben's f and u)   !
  !                                                         !
  !---------------------------------------------------------!
  do ig = amg % imin(level), amg % imax(level)
    i = ig - amg % imin(level) + 1  ! local unkown number
    f(ig) = b(i)
    u(ig) = x(i)
  end do


  call amg % timer_stop(13)

  ia(amg % imax(level)+1) = iaux

  end subroutine
