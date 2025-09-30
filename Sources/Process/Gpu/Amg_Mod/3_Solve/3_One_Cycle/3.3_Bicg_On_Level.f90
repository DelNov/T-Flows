!==============================================================================!
  subroutine Bicg_On_Level(Amg, level, max_iter)
!------------------------------------------------------------------------------!
!   Bi-conjugate gradient on the coarsest level level
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level, max_iter
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  real    :: s
  integer :: ig, jg
!---------------------------------[New locals]---------------------------------!
  integer                      :: i, iter, j, ij, n, nnz
  real,    allocatable         :: m(:)
  real,    allocatable         :: p(:),   q(:),   r(:),   z(:)
  real,    allocatable         :: p_t(:), q_t(:), r_t(:), z_t(:)
  real,    contiguous, pointer :: a(:), u(:), f(:)
  integer, contiguous, pointer :: ia(:), ja(:)
  real,    allocatable         :: a_t(:)
  integer, allocatable         :: ia_t(:), ja_t(:)
  integer, allocatable         :: counter(:)
  real                         :: alpha, beta, pq
  real                         :: res_ini, res_cur  ! don't compare rho and res
  real                         :: rho_old, rho_new  ! they're not the same thing
!------------------------------------[Save]------------------------------------!
  save  ! this is really needed for local allocatable arrays
!==============================================================================!

  ! Initialize residuals
  res_ini = 0.0
  res_cur = 0.0
  rho_old = 0.0
  rho_new = 0.0

  !---------------------------------------------!
  !                                             !
  !                                             !
  !   Create local linear system of equations   !
  !                                             !
  !                                             !
  !---------------------------------------------!
  n   =  Amg % lev(level) % n
  nnz =  Amg % lev(level) % nnz
  a   => Amg % lev(level) % a
  u   => Amg % lev(level) % u
  f   => Amg % lev(level) % f
  ia  => Amg % lev(level) % ia
  ja  => Amg % lev(level) % ja

  !---------------------------------------------------------!
  !                                                         !
  !   Allocate memory for local vectors for the algorithm   !
  !                                                         !
  !---------------------------------------------------------!
  call Amg % Enlarge_Real(m, n)
  call Amg % Enlarge_Real(p, n);  call Amg % Enlarge_Real(p_t, n)
  call Amg % Enlarge_Real(q, n);  call Amg % Enlarge_Real(q_t, n)
  call Amg % Enlarge_Real(r, n);  call Amg % Enlarge_Real(r_t, n)
  call Amg % Enlarge_Real(z, n);  call Amg % Enlarge_Real(z_t, n)

  !----------------------------------------------!
  !                                              !
  !   Create transpose of the local CRS matrix   !
  !                                              !
  !----------------------------------------------!
  call Amg % Enlarge_Real(a_t,  nnz)  ! ;  a_t(:)  = 0.0
  call Amg % Enlarge_Int (ja_t, nnz)  ! ;  ja_t(:) = 0
  call Amg % Enlarge_Int (ia_t, n+1)  ! ;  ia_t(:) = 0

  !----------------------------------------------------!
  !   Phase 1 - count the columns in original matrix   !
  !----------------------------------------------------!
  call Amg % Enlarge_Int(counter, nnz)
  counter(:) = 0  ! acts as column counter

  do i = 1, n                      ! rows in the original matrix
    do ij = ia(i), ia(i+1) - 1
      j = ja(ij)                   ! columns in the original matrix
      counter(j) = counter(j) + 1  ! increase the column count
    end do
  end do

  if(DEBUG) then
    print *, "Size of each column:"
    print '(64i4)', counter(1:n)
  end if

  !---------------------------------------------------!
  !   Phase 2 - form the ia_t based on counter   !
  !---------------------------------------------------!
  ia_t(1) = 1  ! must start at 1
  do i = 2, n + 1   ! ends at n+1 by convention
    ia_t(i) = ia_t(i-1) + counter(i-1)
  end do

  if(DEBUG) then
    print *, "Row pointer for transposed matrix:"
    print '(64i4)', ia_t(1:n+1)
  end if

  !-------------------------------------------!
  !   Phase 3 - fill the rest of the matrix   !
  !-------------------------------------------!
  counter(:) = 0      ! re-initialize the column/row counter

  ! Phase 3.1 - fill the diagonals up
  do i = 1, n       ! rows in the original matrix
    ij = ia(i)
    j = ja(ij)      ! column in the original matrix
                    ! but row in the transposed
    a_t  (ia_t(j)) = a(ij)
    ja_t(ia_t(j)) = i
    counter(j) = 1  ! first entry is the diagonal
  end do

  ! Phase 3.2 - populate the off-diagonal terms
  do i = 1, n         ! rows in the original matrix
    do ij = ia(i) + 1, ia(i+1) - 1
      j = ja(ij)      ! columns in the original matrix
                      ! but rows in the transposed
      a_t  (ia_t(j) + counter(j)) = a(ij)
      ja_t(ia_t(j) + counter(j)) = i
      counter(j) = counter(j) + 1  ! update the counter
    end do
  end do

!#  !--------------------------------------!
!#  !   Phase 4 - validate the transpose   !
!#  !--------------------------------------!
!#
!#  ! Fill p with something non-trivial
!#  do i = 1, n
!#    p(i) = real(i) * 0.1
!#  end do
!#
!#  ! Compute q = A*p
!#  do i = 1, n
!#    do ij = ia(i), ia(i+1)-1
!#      j = ja(ij)
!#      q(i) = q(i) + a(ij) * p(j)
!#    end do
!#  end do
!#
!#  ! Compute r = A^T*p
!#  do i = 1, n
!#    do ij = ia_t(i), ia_t(i+1)-1
!#      j = ja_t(ij)
!#      r(i) = r(i) + a_t(ij) * p(j)
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
    m(i) = 1.0 / a(ia(i))  ! diagonal is at a(ia(i))
  end do

  !--------------------------------------------------------------!
  !   Compute r^(0) = b - A x^(0) for some initial guess x^(0)   !
  !           _                  _                               !
  !   Choose  r^(0) (for example r^(0) = r^(0))                  !
  !--------------------------------------------------------------!
  res_ini = 0.0
  do i = 1, n
    s = f(i)
    do ij = ia(i), ia(i+1) - 1
      j = ja(ij)
      s = s - a(ij) * u(j)
    end do
    r  (i) = s
    r_t(i) = r(i)
    res_ini = res_ini + r(i) * r(i)
  end do
  res_cur = res_ini

  if(DEBUG) then
    print *, "res_ini = ", sqrt(res_ini)
  end if

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
    rho_new = 0.0
    do i = 1, n
      rho_new = rho_new + r_t(i) * z(i)
    end do

  if(DEBUG) then
    print *, "rho_new = ", sqrt(rho_new)
  end if

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
      s = 0.0
      do ij = ia(i), ia(i+1) - 1
        j = ja(ij)
        s = s + a(ij) * p(j)
      end do
      q(i) = s
    end do

    !---_-------------_-------!
    !   q^(i) = A^T * p^(i)   !
    !-------------------------!
    do i = 1, n
      s = 0.0
      do ij = ia_t(i), ia_t(i+1) - 1
        j = ja_t(ij)
        s = s + a_t(ij) * p_t(j)
      end do
      q_t(i) = s
    end do

    !---------------------_---------------!
    !   alpha = rho_(i-1)/p^(i)^T q^(i)   !
    !-------------------------------------!
    pq = 0.0
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
      u  (i) = u  (i) + alpha * p  (i)
      r  (i) = r  (i) - alpha * q  (i)
      r_t(i) = r_t(i) - alpha * q_t(i)
    end do

    ! Compute current residual
    res_cur = 0.0
    do i = 1, n
      s = f(i)
      do ij = ia(i), ia(i+1) - 1
        j = ja(ij)
        s = s - a(ij) * u(j)
      end do
      res_cur = res_cur + s * s
    end do

  if(DEBUG) then
    print '(a,i3,a,1pe14.7)',  &
      "iter: ", iter, " res_new/res_ini = ", sqrt(res_cur/res_ini)
  end if
    if(sqrt(res_cur) .lt. Amg % eps) exit

    rho_old = rho_new

  end do
  end if

  ! Print final residual
  if(DEBUG) then
    print '(a,i3,a,1pe14.7)',  &
      "iter: ", iter, " res_fin = ", sqrt(res_cur)
  end if
  if(Amg % iout .gt. 3) then
    write(*, '(a, 1es12.3)')  ' res_fin = ', sqrt(res_cur)
  end if

  end subroutine
