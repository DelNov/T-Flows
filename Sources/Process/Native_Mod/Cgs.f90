!==============================================================================!
  subroutine Cgs(Nat, A, x, b, prec, miter, niter, tol, fin_res, norm)
!------------------------------------------------------------------------------!
!   Solves the linear systems of equations by a preconditioned CGS method      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type), target :: Nat
  type(Matrix_Type)          :: A
  real                       :: x(-Nat % pnt_grid % n_bnd_cells :  &
                                   Nat % pnt_grid % n_cells)
  real                       :: b( Nat % pnt_grid % n_cells)
  character(SL)              :: prec     ! preconditioner
  integer                    :: miter    ! maximum and actual ...
  integer                    :: niter    ! ... number of iterations
  real                       :: tol      ! tolerance
  real                       :: fin_res  ! final residual
  real, optional             :: norm     ! normalization
!-----------------------------------[Locals]-----------------------------------!
  type(Matrix_Type), pointer :: D
  integer                    :: nt, ni, nb
  real                       :: alfa, beta, rho, rho_old, bnrm2, res
  integer                    :: i, j, k, iter
  real                       :: sum_a, fn
  integer                    :: sum_n
  real, contiguous, pointer  :: p1(:), q1(:), r1(:), p2(:), q2(:), r2(:),  &
                                u1(:), u2(:), v2(:), u1_plus_q1(:)
!==============================================================================!

  call Work % Connect_Real_Cell(p1, q1, r1, p2, q2, r2,  &
                                u1, u2, v2, u1_plus_q1)

  ! Take some aliases
  D => Nat % D
  nt = A % pnt_grid % n_cells
  ni = A % pnt_grid % n_cells - A % pnt_grid % comm % n_buff_cells
  nb = A % pnt_grid % n_bnd_cells

  res = 0.0

  !--------------------------!
  !   Normalize the system   !
  !--------------------------!
  sum_a = 0.0
  sum_n = 0
  do i = 1, ni
    sum_a = sum_a + A % val(A % dia(i))
    sum_n = sum_n + 1
  end do
  call Comm_Mod_Global_Sum_Real(sum_a)
  call Comm_Mod_Global_Sum_Int (sum_n)  ! this is stored somewhere, check
  sum_a = sum_a / sum_n
  fn = 1.0 / sum_a
  do i = 1, nt
    do j = A % row(i), A % row(i+1)-1
      A % val(j) = A % val(j) * fn
    end do
    b(i) = b(i) * fn
  end do

  !---------------------!
  !   Preconditioning   !
  !---------------------!
  call Nat % Prec_Form(ni, A, D, prec)

  !-----------------------------------!
  !    This is quite tricky point.    !
  !   What if bnrm2 is very small ?   !
  !-----------------------------------!
  if(.not. present(norm)) then
    bnrm2 = Nat % Normalized_Root_Mean_Square(ni, b(1:nt), A, x(1:nt))
  else
    bnrm2 = Nat % Normalized_Root_Mean_Square(ni, b(1:nt), A, x(1:nt), norm)
  end if

  if(bnrm2 < tol) then
    iter=0
    goto 1
  end if

  !-----------------!
  !   r1 = b - Ax   !
  !-----------------!
  call Nat % Residual_Vector(ni, r1(1:nt), b(1:nt), A, x(1:nt))

  !--------------------------------!
  !   Calculate initial residual   !
  !--------------------------------!
  res = Nat % Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt))

  if(res < tol) then
    iter = 0
    goto 1
  end if

  !-------------!
  !   r2 = r1   !
  !-------------!
  do i = 1, ni
    r2(i) = r1(i)
  end do

  !---------------!
  !               !
  !   Main loop   !
  !               !
  !---------------!
  do iter = 1, miter

    !-------------------!
    !   rho = (r2,z1)   !
    !-------------------!
    rho = dot_product(r1(1:ni), r2(1:ni))
    call Comm_Mod_Global_Sum_Real(rho)

    if(iter .eq. 1) then
      do i = 1, ni
        u1(i) = r1(i)
        u2(i) = u1(i)
      end do
    else
      beta = rho / rho_old
      do i = 1, ni
        u1(i) = r1(i) + beta *  q1(i)
        u2(i) = u1(i) + beta * (q1(i) + beta*u2(i))
      end do
    end if

    !---------------------!
    !   Solve M p2 = u2   !
    !---------------------!
    call Nat % Prec_Solve(ni, A, D, p2(1:nt), u2(1:nt), prec)

    !--------------!
    !   v2 = Ap2   !
    !--------------!
    call A % pnt_grid % Exchange_Cells_Real(p2(-nb:ni))
    do i = 1, ni
      v2(i) = 0.0
      do j = A % row(i), A % row(i+1)-1
        k = A % col(j)
        v2(i) = v2(i) + A % val(j) * p2(k)
      end do
      alfa = alfa + r2(i) * v2(i)
    end do

    !------------------------!
    !   alfa = rho/(r2,v2)   !
    !------------------------!
    alfa = dot_product(r2(1:ni), v2(1:ni))
    call Comm_Mod_Global_Sum_Real(alfa)
    alfa=rho / alfa

    !-------------------------!
    !   q1 = u1 - alfa * v2   !
    !-------------------------!
    do i = 1, ni
      q1(i) = u1(i) - alfa * v2(i)
    end do

    !-------------------------!
    !   solve Mp1 = u1 + q1   !
    !-------------------------!
    do i = 1, ni
      u1_plus_q1(i) = u1(i) + q1(i)
    end do
    call Nat % Prec_Solve(ni, A, D, p1(1:nt), u1_plus_q1(1:nt), prec)

    !---------------------!
    !   x = x + alfa p1   !
    !---------------------!
    do i = 1, ni
      x(i) = x(i) + alfa * p1(i)
    end do

    !---------------!
    !   q2 = A p1   !
    !---------------!
    call A % pnt_grid % Exchange_Cells_Real(p1(-nb:ni))
    do i = 1, ni
      q2(i) = 0.0
      do j = A % row(i), A % row(i+1)-1
        k = A % col(j)
        q2(i) = q2(i) + A % val(j) * p1(k)
      end do
    end do

    !---------------------!
    !   r = r - alfa q2   !
    !---------------------!
    do i = 1, ni
      r1(i) = r1(i) - alfa * q2(i)
    end do

    !-----------------------!
    !   Check convergence   !
    !-----------------------!
    if(.not. present(norm)) then
      res = Nat % Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt))
    else
      res = Nat % Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt), norm)
    end if

    if(res < tol) goto 1

    rho_old = rho

  end do  ! iter

  ! Correct the number of executed iterations, because
  ! Fortran goes one number up do loop's upper limit
  iter = iter - 1

  !----------------------------------!
  !                                  !
  !   Convergence has been reached   !
  !                                  !
  !----------------------------------!
1 continue

  !-------------------------------------------!
  !   Refresh the solution vector's buffers   !
  !-------------------------------------------!
  call A % pnt_grid % Exchange_Cells_Real(x(-nb:ni))

  !-----------------------------!
  !   De-normalize the system   !
  !-----------------------------!
  do i = 1, nt
    do j = A % row(i), A % row(i+1)-1
      A % val(j) = A % val(j) / fn
    end do
    b(i) = b(i) / fn
  end do

  fin_res = res
  niter   = iter

  call Work % Disconnect_Real_Cell(p1, q1, r1, p2, q2, r2,  &
                                   u1, u2, v2, u1_plus_q1)

  end subroutine
