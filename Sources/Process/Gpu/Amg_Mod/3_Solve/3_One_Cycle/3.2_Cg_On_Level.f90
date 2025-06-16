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
  class(Amg_Type), target :: Amg
  integer                 :: level, max_iter
  real                    :: a(:), u(:), f(:)
  integer                 :: ia(:), ja(:)
  integer                 :: iw(:), icg(:)
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
  real                 :: alpha, beta, pq
  real                 :: res_ini, res_cur  ! don't mix rhos and res's ...
  real                 :: rho_old, rho_new  ! ... they are not the same thing

  real,    contiguous, pointer :: lev_a(:), lev_u(:), lev_f(:)
  integer, contiguous, pointer :: lev_ia(:), lev_ja(:)
!------------------------------------[save]------------------------------------!
  save  ! this is really needed for local allocatable arrays
!==============================================================================!

  call Amg % timer_start()

  ! Initialize residuals
  res_ini = 0.0
  res_cur = 0.0
  rho_old = 0.0
  rho_new = 0.0

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

!old:  ! Number of unknowns and non-zeros on this level
!old:  n = Amg % imax(level) - Amg % imin(level) + 1
!old:  nnz = ia(Amg % imax(level)+1) - ia(Amg % imin(level))
!old:
!old:  if(Amg % iout .gt. 3) then
!old:    write(*,'(a,i4,a,i9,a)', advance = 'no')  &
!old:      ' # CG on level ', level, ' with ', n, ' unknowns; '
!old:  end if

  n      =  Amg % lev(level) % n
  nnz    =  Amg % lev(level) % nnz
  lev_a  => Amg % lev(level) % a
  lev_u  => Amg % lev(level) % u
  lev_f  => Amg % lev(level) % f
  lev_ia => Amg % lev(level) % ia
  lev_ja => Amg % lev(level) % ja

  !-------------------------------------------------------!
  !                                                       !
  !   Create local CRS matrix from the global ia, ja, a   !
  !                                                       !
  !-------------------------------------------------------!
  call Amg % Enlarge_Real(a_val,   nnz);  a_val(:)   = 0.0
  call Amg % Enlarge_Int (col_idx, nnz);  col_idx(:) = 0
  call Amg % Enlarge_Int (row_ptr, n+1);  row_ptr(:) = 0

!old:  do ig = Amg % imin(level), Amg % imax(level)
!old:    i = ig - Amg % imin(level) + 1                   ! local unknown number
!old:    row_ptr(i) = ia(ig) - ia(Amg % imin(level)) + 1  ! local row pointer
!old:    do jg = ia(ig), ia(ig+1) - 1                     ! browse through row ig
!old:      j          = jg - ia(Amg % imin(level)) + 1    ! local nonzero index
!old:      col_idx(j) = ja(jg) - Amg % imin(level) + 1    ! local column number
!old:      a_val(j)   = a(jg)                             ! local matrix value
!old:    end do
!old:  end do
!old:  ! Final row pointer
!old:  row_ptr(n+1) = ia(Amg % imax(level)+1) - ia(Amg % imin(level)) + 1

  do i = 1, n + 1
    row_ptr(i) = lev_ia(i)
  end do

  do k = 1, nnz
    col_idx(k) = lev_ja(k)
    a_val(k)   = lev_a(k)
  end do

  !------------------------------------------------------!
  !                                                      !
  !   Create local b and x vectors (Stueben's f and u)   !
  !                                                      !
  !------------------------------------------------------!
  call Amg % Enlarge_Real(b, n)
  call Amg % Enlarge_Real(x, n)
  call Amg % Update_U_And_F_At_Level(level, vec_u=u, vec_f=f)

!old  do ig = Amg % imin(level), Amg % imax(level)
!old    i = ig - Amg % imin(level) + 1  ! local unkown number
!old    b(i) = f(ig)
!old    x(i) = u(ig)
!old  end do

  do i = 1, n
    b(i) = lev_f(i)
    x(i) = lev_u(i)
  end do

  !---------------------------------------------------------!
  !                                                         !
  !   Allocate memory for local vectors for the algorithm   !
  !                                                         !
  !---------------------------------------------------------!
  call Amg % Enlarge_Real(m, n)
  call Amg % Enlarge_Real(p, n)
  call Amg % Enlarge_Real(q, n)
  call Amg % Enlarge_Real(r, n)
  call Amg % Enlarge_Real(z, n)
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
  res_ini = 0.0
  do i = 1, n
    s = b(i)
    do k = row_ptr(i), row_ptr(i+1) - 1
      j = col_idx(k)
      s = s - a_val(k) * x(j)
    end do
    r(i) = s
    res_ini = res_ini + r(i) * r(i)
  end do
  res_cur = res_ini

# ifdef DEBUG
    print *, "res_ini = ", sqrt(res_ini)
# endif

  if(Amg % iout .gt. 3) then
    write(*,'(a,1es12.3,a)', advance = 'no')  &
      ' res_ini = ', sqrt(res_ini), '; '
  end if


  !--------------------!
  !                    !
  !   Iteration loop   !
  !                    !
  !--------------------!
  if(sqrt(res_ini) .gt. Amg % eps) then
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

# ifdef DEBUG
    print *, "rho_new = ", sqrt(rho_new)
# endif

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

    ! Compute current residual
    res_cur = 0.0
    do i = 1, n
      s = b(i)
      do k = row_ptr(i), row_ptr(i+1) - 1
        j = col_idx(k)
        s = s - a_val(k) * x(j)
      end do
      res_cur = res_cur + s * s
    end do

# ifdef DEBUG
    print '(a,i3,a,1pe14.7)',  &
      "iter: ", iter, " res_new/res_ini = ", sqrt(res_cur/res_ini)
# endif
    if(sqrt(res_cur) .lt. Amg % eps) exit

    rho_old = rho_new

  end do
  end if

  ! Print final residual
# ifdef DEBUG
    print '(a,i3,a,1pe14.7)',  &
      "iter: ", iter, " res_fin = ", sqrt(res_cur)
# endif
  if(Amg % iout .gt. 3) then
    write(*, '(a, 1es12.3)')  ' res_fin = ', sqrt(res_cur)
  end if

  !---------------------------------------------------------!
  !                                                         !
  !   Copy local b and x vectors back (Stueben's f and u)   !
  !                                                         !
  !---------------------------------------------------------!
!old:  do ig = Amg % imin(level), Amg % imax(level)
!old:    i = ig - Amg % imin(level) + 1  ! local unkown number
!old:    f(ig) = b(i)
!old:    u(ig) = x(i)
!old:  end do

  do i = 1, n
    lev_f(i) = b(i)
    lev_u(i) = x(i)
  end do
  call Amg % Update_U_And_F_Globally(level, vec_u=u, vec_f=f)

  call Amg % timer_stop(13)

  ia(Amg % imax(level)+1) = iaux

  end subroutine
