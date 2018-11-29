!==============================================================================!
  subroutine Cgs(sol, x, b, prec, niter, tol, ini_res, fin_res, norm)
!------------------------------------------------------------------------------!
!   Solves the linear systems of equations by a precond. CGS Method.           !
!------------------------------------------------------------------------------!
!   Allows preconditioning of the system by:                                   !
!     1. Diagonal preconditioning                                              !
!     2. Incomplete Cholesky preconditioning                                   !
!                                                                              !
!   The type of precondtioning is chosen by setting the variable prec to 0     !
!   (for no preconditioning), 1 (for diagonal preconditioning) or 2 (for       !
!   incomplete Cholesky preconditioning)                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Matrix_Mod
  use Work_Mod, only: p1         => r_cell_01,  &
                      p2         => r_cell_02,  &
                      q1         => r_cell_03,  &
                      q2         => r_cell_04,  &
                      r1         => r_cell_06,  &
                      r2         => r_cell_07,  &
                      u1         => r_cell_08,  &
                      u2         => r_cell_09,  &
                      v2         => r_cell_10,  &
                      u1_plus_q1 => r_cell_11
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Solver_Type), target :: sol
  real              :: x(-sol % pnt_grid % n_bnd_cells :  &
                          sol % pnt_grid % n_cells)
  real              :: b( sol % pnt_grid % n_cells)      ! [A]{x}={b}
  character(len=80) :: prec                              ! preconditioner
  integer           :: niter                             ! number of iterations
  real              :: tol                               ! tolerance
  real              :: ini_res, fin_res                  ! residual
  real, optional    :: norm                              ! normalization
!-----------------------------------[Locals]-----------------------------------!
  type(Matrix_Type), pointer :: a
  integer :: nt, ni, nb
  real    :: alfa, beta, rho, rho_old, bnrm2, error
  integer :: i, j, k, iter, sub
!==============================================================================!

  ! Take some aliases
  a => sol % a
  nt = a % pnt_grid % n_cells
  ni = a % pnt_grid % n_cells - a % pnt_grid % comm % n_buff_cells
  nb = a % pnt_grid % n_bnd_cells

  error = 0.0
  r1(:) = 0
  r1(1:nt) = b(1:nt)

  !---------------------!
  !   Preconditioning   !
  !---------------------!
  call Prec_Form(sol, prec)

  !-----------------------------------!
  !    This is quite tricky point.    !
  !   What if bnrm2 is very small ?   !
  !-----------------------------------!
  if(.not. present(norm)) then
    bnrm2 = Normalized_Root_Mean_Square(ni, r1(1:ni), a, x(1:nt))
  else
    bnrm2 = Normalized_Root_Mean_Square(ni, r1(1:ni), a, x(1:nt), norm)
  end if

  if(bnrm2 < tol) then 
    iter=0
    goto 1
  end if  

  !-----------------!
  !   r1 = b - Ax   !
  !-----------------!
  call Residual_Vector(ni, r1(1:ni), r1(1:ni), a, x(1:nt))

  !-------------!
  !   r2 = r1   !
  !-------------!
  do i=1,ni
    r2(i)=r1(i) 
  end do

  !--------------------------------!
  !   Calculate initial residual   !
  !--------------------------------!
  error = Normalized_Root_Mean_Square(ni, r1(1:ni), a, x(1:nt))

  !---------------------------------------------------------------!
  !   Residual after the correction and before the new solution   !
  !---------------------------------------------------------------!
  ini_res = error  

  if(error < tol) then
    iter=0
    goto 1
  end if  

  !---------------!
  !               !
  !   Main loop   !
  !               !
  !---------------!
  do iter = 1, niter 

    !-------------------!
    !   rho = (r2,z1)   !
    !-------------------!
    rho=0.0
    do i=1,ni
      rho=rho+r1(i)*r2(i)
    end do
    call Comm_Mod_Global_Sum_Real(rho)

    if(iter .eq. 1) then
      do i=1,ni
        u1(i) = r1(i)
        u2(i) = u1(i)
      end do        
    else
      beta=rho/rho_old
      do i=1,ni
        u1(i) = r1(i) + beta*q1(i) 
        u2(i) = u1(i) + beta*(q1(i) + beta*u2(i)) 
      end do
    end if
                   
    !---------------------!
    !   Solve M p2 = u2   !
    !---------------------!
    call Prec_Solve(sol, p2, u2(1), prec) 

    !--------------!
    !   v2 = Ap2   !  
    !--------------!
    call Comm_Mod_Exchange_Real(a % pnt_grid, p2)
    do i=1,ni
      v2(i) = 0.0
      do j=a % row(i), a % row(i+1)-1
        k=a % col(j)
        v2(i) = v2(i) + a % val(j) * p2(k)
      end do
      alfa=alfa+r2(i)*v2(i)
    end do

    !------------------------!
    !   alfa = rho/(r2,v2)   !
    !------------------------!
    alfa=0.0
    do i=1,ni
      alfa=alfa+r2(i)*v2(i)
    end do
    call Comm_Mod_Global_Sum_Real(alfa) 
    alfa=rho/alfa

    !-------------------------!
    !   q1 = u1 - alfa * v2   !
    !-------------------------!
    do i=1,ni
      q1(i) = u1(i) - alfa*v2(i)
    end do

    !-------------------------------!
    !   solve Mp1 = u1(i) + q1(i)   !
    !-------------------------------!
    do i=1,ni
      u1_plus_q1(i) = u1(i) + q1(i)
    end do
    call Prec_Solve(sol, p1, u1_plus_q1(1), prec) 

    !---------------------!
    !   x = x + alfa p1   !
    !---------------------!
    do i=1,ni
      x(i)=x(i) + alfa*p1(i)
    end do

    !---------------!
    !   q2 = A p1   !
    !---------------!
    call Comm_Mod_Exchange_Real(a % pnt_grid, p1)
    do i=1,ni
      q2(i) = 0.0
      do j=a % row(i), a % row(i+1)-1
        k=a % col(j)
        q2(i) = q2(i) + a % val(j) * p1(k)
      end do
    end do

    !---------------------!
    !   r = r - alfa q2   !
    !---------------------!
    do i=1,ni
      r1(i)=r1(i) - alfa*q2(i)
    end do

    !-----------------------!
    !   Check convergence   !
    !-----------------------!
    if(.not. present(norm)) then
      error = Normalized_Root_Mean_Square(ni, r1(1:ni), a, x(1:nt))
    else
      error = Normalized_Root_Mean_Square(ni, r1(1:ni), a, x(1:nt), norm)
    end if

    if(error < tol) goto 1

    rho_old=rho

  end do                ! iter 

1 fin_res = error
  niter  = iter

  end subroutine
