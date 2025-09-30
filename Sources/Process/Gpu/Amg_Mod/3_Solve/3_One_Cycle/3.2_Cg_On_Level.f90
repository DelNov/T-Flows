!==============================================================================!
  subroutine Cg_On_Level(Amg, level, max_iter)
!------------------------------------------------------------------------------!
!   Conjugate gradient on the coarsest level level, as described here:
!   w3.pppl.gov/~hammett/comp/numerical_tricks/templates.pdf
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
! real,    allocatable         :: p_t(:), q_t(:), r_t(:), z_t(:)
  real,    contiguous, pointer :: a(:), u(:), f(:)
  integer, contiguous, pointer :: ia(:), ja(:)
! real,    allocatable         :: a_t(:)
! integer, allocatable         :: ia_t(:), ja_t(:)
! integer, allocatable         :: counter(:)
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
  call Amg % Enlarge_Real(p, n)
  call Amg % Enlarge_Real(q, n)
  call Amg % Enlarge_Real(r, n)
  call Amg % Enlarge_Real(z, n)

  !---------------------------------!
  !   Form preconditioning matrix   !
  !---------------------------------!
  do i = 1, n
    m(i) = 1.0 / a(ia(i))  ! diagonal is at a(ia(i))
  end do

  !--------------------------------------------------------------!
  !   Compute r^(0) = b - A x^(0) for some initial guess x^(0)   !
  !--------------------------------------------------------------!
  res_ini = 0.0
  do i = 1, n
    s = f(i)
    do ij = ia(i), ia(i+1) - 1
      j = ja(ij)
      s = s - a(ij) * u(j)
    end do
    r(i) = s
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

  if(DEBUG) then
    print *, "rho_new = ", sqrt(rho_new)
  end if

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
      do ij = ia(i), ia(i+1) - 1
        j = ja(ij)
        s = s + a(ij) * p(j)
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
      u(i) = u(i) + alpha * p(i)
      r(i) = r(i) - alpha * q(i)
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
