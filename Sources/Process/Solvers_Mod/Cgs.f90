!==============================================================================!
  subroutine Cgs(mat_a, x, r1, prec, niter, tol, ini_res, fin_res)
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
                      r2         => r_cell_05,  &
                      u1         => r_cell_06,  &
                      u2         => r_cell_07,  &
                      v2         => r_cell_08,  &   
                      u1_plus_q1 => r_cell_09
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Matrix_Type) :: mat_a           
  real    :: x(-mat_a % pnt_grid % n_bnd_cells : mat_a % pnt_grid % n_cells)
  real    :: r1(mat_a % pnt_grid % n_cells)    !  [A]{x}={r1}
  integer :: niter              ! number of iterations
  real    :: tol                ! tolerance
  real    :: ini_res, fin_res   ! residual
  character(len=80) :: prec     ! preconditioner
!-----------------------------------[Locals]-----------------------------------!
  integer :: n, nb
  real    :: alfa, beta, rho, rho_old, bnrm2, error
  integer :: i, j, k, iter, sub
!==============================================================================!
           
  n  = mat_a % pnt_grid % n_cells
  nb = mat_a % pnt_grid % n_bnd_cells

  !---------------------!
  !   Preconditioning   !
  !---------------------!
  call Prec_Form(mat_a, prec)    

  !???????????????????????????????????!
  !    This is quite tricky point.    !
  !   What if bnrm2 is very small ?   !
  !???????????????????????????????????!
  bnrm2 = Normalized_Residual(n, nb, mat_a, x, r1)

  if(bnrm2 < tol) then 
    iter=0
    goto 1
  end if  

  !-----------------!
  !   r1 = b - Ax   !
  !-----------------!
  call Residual_Vector(n, nb, mat_a, x, r1) 

  !-------------!
  !   r2 = r1   !
  !-------------!
  do i=1,n
    r2(i)=r1(i) 
  end do

  !--------------------------------!
  !   Calculate initial residual   !
  !--------------------------------!
  error = Normalized_Residual(n, nb, mat_a, x, r1)

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
    do i=1,n
      rho=rho+r1(i)*r2(i)
    end do
    call Comm_Mod_Global_Sum_Real(rho)

    if(iter .eq. 1) then
      do i=1,n
        u1(i) = r1(i)
        u2(i) = u1(i)
      end do        
    else
      beta=rho/rho_old
      do i=1,n
        u1(i) = r1(i) + beta*q1(i) 
        u2(i) = u1(i) + beta*(q1(i) + beta*u2(i)) 
      end do
    end if
                   
    !---------------------!
    !   Solve M p2 = u2   !
    !---------------------!
    call Prec_Solve(mat_a, p2, u2(1), prec) 

    !--------------!
    !   v2 = Ap2   !  
    !--------------!
    do i=1,n
      v2(i) = 0.0                    
      do j=mat_a % row(i), mat_a % row(i+1)-1   
        k=mat_a % col(j)                  
        v2(i) = v2(i) + mat_a % val(j) * p2(k)   
      end do
      alfa=alfa+r2(i)*v2(i)
    end do
    call Comm_Mod_Exchange(mat_a % pnt_grid, p2)
    do sub=1,n_proc
      if(nbb_e(sub)  <=  nbb_s(sub)) then
        do k=nbb_s(sub),nbb_e(sub),-1
          i=buffer_index(k)
          v2(i) = v2(i) + mat_a % bou(k)*p2(k)
        end do
      end if
    end do

    !------------------------!
    !   alfa = rho/(r2,v2)   !
    !------------------------!
    alfa=0.0
    do i=1,n
      alfa=alfa+r2(i)*v2(i)
    end do
    call Comm_Mod_Global_Sum_Real(alfa) 
    alfa=rho/alfa

    !-------------------------!
    !   q1 = u1 - alfa * v2   !
    !-------------------------!
    do i=1,n
      q1(i) = u1(i) - alfa*v2(i)
    end do
           
    !-------------------------------!
    !   solve Mp1 = u1(i) + q1(i)   !
    !-------------------------------!
    do i=1,n
      u1_plus_q1(i) = u1(i) + q1(i)
    end do
    call Prec_Solve(mat_a, p1, u1_plus_q1(1), prec) 

    !---------------------!
    !   x = x + alfa p1   !
    !---------------------!
    do i=1,n
      x(i)=x(i) + alfa*p1(i)
    end do

    !---------------!
    !   q2 = A p1   !     
    !---------------!
    do i=1,n
      q2(i) = 0.0
      do j=mat_a % row(i), mat_a % row(i+1)-1
        k=mat_a % col(j)
        q2(i) = q2(i) + mat_a % val(j) * p1(k)
      end do
    end do
    call Comm_Mod_Exchange(mat_a % pnt_grid, p1)
    do sub=1,n_proc
      if(nbb_e(sub)  <=  nbb_s(sub)) then
        do k=nbb_s(sub),nbb_e(sub),-1
          i=buffer_index(k)
          q2(i) = q2(i) + mat_a % bou(k)*p1(k)
        end do
      end if
    end do

    !---------------------!
    !   r = r - alfa q2   !
    !---------------------!
    do i=1,n
      r1(i)=r1(i) - alfa*q2(i)
    end do

    !???????????????????????!
    !   Check convergence   !
    !???????????????????????!
    error = Normalized_Residual(n, nb, mat_a, x, r1)

    if(error < tol) goto 1

    rho_old=rho

  end do                ! iter 

1 fin_res = error
  Niter  = iter

  end subroutine
