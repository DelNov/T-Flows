! #define DEBUG  1

!==============================================================================!
  subroutine Cg_On_Level(Amg, level, max_iter,  &
                         a, u, f, ia, ja,       &  ! defines system
                         iw, icg)
!------------------------------------------------------------------------------!
!   Conjugate gradient on the coarsest level level, as described here:
!   w3.pppl.gov/~hammett/comp/numerical_tricks/templates.pdf
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: level, max_iter
  real            :: a(:), u(:), f(:)
  integer         :: ia(:), ja(:)
  integer         :: iw(:), icg(:)
!-----------------------------------[locals]-----------------------------------!
  real    :: s
  integer :: ig, iaux, jg
!---------------------------------[new locals]---------------------------------!
  integer              :: i, iter, j, k, n, nnz
  real,    allocatable :: m(:), b(:), x(:)
  real,    allocatable :: p(:),   q(:),   r(:),   z(:)
! real,    allocatable :: p_t(:), q_t(:), r_t(:), z_t(:)
  real,    allocatable :: a_val(:)
  integer, allocatable :: row_ptr(:), col_idx(:)
! real,    allocatable :: a_t_val(:)
! integer, allocatable :: row_t_ptr(:), col_t_idx(:)
! integer, allocatable :: counter(:)
  real                 :: alpha, beta
  real                 :: pq, rho_0, rho_old, rho_new
!==============================================================================!

  call Amg % timer_start()

  ! See comment in source "Coarsening.f90" at line 180
  iaux = ia(Amg % imax(level)+1)
  ia(Amg % imax(level)+1) = iw(Amg % iminw(level))


  !---------------------------------------------!
  !                                             !
  !                                             !
  !   Create local linear system of equations   !
  !                                             !
  !                                             !
  !---------------------------------------------!

  ! Number of unknowns and non-zeros on this level
  n = Amg % imax(level) - Amg % imin(level) + 1
  nnz = ia(Amg % imax(level)+1) - ia(Amg % imin(level))

# ifdef AMG_VERBOSE
    write(*,'(a,i4,a,i9,a)', advance = 'no')  &
      ' # CG on level ', level, ' with ', n, ' unknowns; '
# endif

  !-------------------------------------------------------!
  !                                                       !
  !   Create local CRS matrix from the global ia, ja, a   !
  !                                                       !
  !-------------------------------------------------------!
  allocate(a_val(nnz));   a_val(:)   = 0.0
  allocate(col_idx(nnz)); col_idx(:) = 0
  allocate(row_ptr(n+1)); row_ptr(:) = 0
  do ig = Amg % imin(level), Amg % imax(level)
    i = ig - Amg % imin(level) + 1                   ! local unknown number
    row_ptr(i) = ia(ig) - ia(Amg % imin(level)) + 1  ! local row pointer
    do jg = ia(ig), ia(ig+1) - 1              ! browse through row ig
      j          = jg - ia(Amg % imin(level)) + 1    ! local nonzero index
      col_idx(j) = ja(jg) - Amg % imin(level) + 1    ! local column number
      a_val(j)   = a(jg)                      ! local matrix value
    end do
  end do
  row_ptr(n+1) = ia(Amg % imax(level)+1) - ia(Amg % imin(level)) + 1  ! final row pointer

  !------------------------------------------------------!
  !                                                      !
  !   Create local b and x vectors (Stueben's f and u)   !
  !                                                      !
  !------------------------------------------------------!
  allocate(b(n))
  allocate(x(n))
  do ig = Amg % imin(level), Amg % imax(level)
    i = ig - Amg % imin(level) + 1  ! local unkown number
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
  m(:) = 0.0
  p(:) = 0.0
  q(:) = 0.0
  r(:) = 0.0
  z(:) = 0.0

  !---------------------------------!
  !   Form preconditioning matrix   !
  !---------------------------------!
  do i = 1, n
    m(i) = 1.0 / a_val(row_ptr(i))  ! diagonal is at a_val(row_ptr(i))
  end do

  !--------------------------------------------------------------!
  !   Compute r^(0) = b - A x^(0) for some initial guess x^(0)   !
  !--------------------------------------------------------------!
  rho_0 = 0.0
  do i = 1, n
    s = b(i)
    do k = row_ptr(i), row_ptr(i+1) - 1
      j = col_idx(k)
      s = s - a_val(k) * x(j)
    end do
    r(i) = s
    rho_0 = rho_0 + r(i) * r(i)
  end do

# ifdef DEBUG
    print *, "r_ini = ", sqrt(rho_0)
# endif

# ifdef AMG_VERBOSE
    write(*,'(a,1es12.3,a)', advance = 'no')  &
      ' res_ini = ', sqrt(rho_0), '; '
# endif


  !--------------------!
  !                    !
  !   Iteration loop   !
  !                    !
  !--------------------!
  if(sqrt(rho_0) .gt. Amg % eps) then
  do iter = 1, max_iter

    !-------------------------------!
    !   Solve M z^(i-1) = r^(i-1)   !
    !-------------------------------!
    do i = 1, n
      z(i) = m(i) * r(i)
    end do

    !-----------------------------------!
    !   rho_(i-1) = r^(i-1)^T z^(i-1)   !
    !-----------------------------------!
    rho_new = 0.0
    do i = 1, n
      rho_new = rho_new + r(i) * z(i)
    end do

    if(iter .eq. 1) then
      !-------------------!
      !   p^(1) = z^(0)   !
      !-------------------!
      do i = 1, n
        p(i) = z(i)
      end do
    else
      !--------------------------------------------!
      !   beta_(i-1) = rho_(i-1) / rho_(i-2)       !
      !                                            !
      !   p^(i) = z^(i-1) + beta_(i-1) * p^(i-1)   !
      !--------------------------------------------!
      beta = rho_new / rho_old
      do i = 1, n
        p(i) = z(i) + beta * p(i)
      end do
    end if

    !-----------------------!
    !   q^(i) = A * p^(i)   !
    !-----------------------!
    do i = 1, n
      s = 0.0
      do k = row_ptr(i), row_ptr(i+1) - 1
        j = col_idx(k)
        s = s + a_val(k) * p(j)
      end do
      q(i) = s
    end do

    !-------------------------------------!
    !   alpha = rho_(i-1)/p^(i)^T q^(i)   !
    !-------------------------------------!
    pq = 0.0
    do i = 1, n
      pq = pq + p(i) * q(i)
    end do
    alpha = rho_new / pq

    !-------------------------------------!
    !   x^(i) = x^(i-1) + alpha * p^(i)   !
    !   r^(i) = r^(i-1) - alpha * q^(i)   !
    !-------------------------------------!
    do i = 1, n
      x(i) = x(i) + alpha * p(i)
      r(i) = r(i) - alpha * q(i)
    end do

# ifdef DEBUG
    print '(a,i3,a,1pe14.7)',  &
      "iter: ", iter, " r_new/r_ini = ", sqrt(rho_new/rho_0)
# endif
    if(sqrt(rho_new) .lt. Amg % eps) exit

    rho_old = rho_new

  end do
  end if

  ! Print final residual
# ifdef DEBUG
    print '(a,i3,a,1pe14.7)',  &
      "iter: ", iter, " r_new = ", sqrt(rho_new)
# endif
# ifdef AMG_VERBOSE
    write(*, '(a, 1es12.3)')  ' res_fin = ', sqrt(rho_new)
# endif

  !---------------------------------------------------------!
  !                                                         !
  !   Copy local b and x vectors back (Stueben's f and u)   !
  !                                                         !
  !---------------------------------------------------------!
  do ig = Amg % imin(level), Amg % imax(level)
    i = ig - Amg % imin(level) + 1  ! local unkown number
    f(ig) = b(i)
    u(ig) = x(i)
  end do


  call Amg % timer_stop(13)

  ia(Amg % imax(level)+1) = iaux

  end subroutine
