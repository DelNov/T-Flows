!==============================================================================!
  subroutine Cg_Level(lev, a, d, x, r1, niter, tol, ini_res, fin_res)
!------------------------------------------------------------------------------!
!   Conjugate gradient method for one level of the multigrid.                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Matrix_Mod
  use Work_Mod, only: p1 => r_cell_01,  &
                      q1 => r_cell_02
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: lev
  type(Matrix_Type) :: a
  type(Matrix_Type) :: d
  real              :: x (a % pnt_grid % level(lev) % n_cells)
! real              :: x(-a % pnt_grid % n_bnd_cells :  &
!                         a % pnt_grid % n_cells)
  real              :: r1(a % pnt_grid % level(lev) % n_cells)  ! [A]{x}={r1}
  integer           :: niter                             ! number of iterations
  real              :: tol                               ! tolerance
  real              :: ini_res                           ! initial residual
  real              :: fin_res                           ! final residual
!-----------------------------------[Locals]-----------------------------------!
  integer                    :: nt, ni
  real                       :: alfa, beta, rho, rho_old, bnrm2, res
  integer                    :: i, j, k, iter, sub
!==============================================================================!

  ! Take some aliases
  nt = a % pnt_grid % level(lev) % n_cells
  ni = a % pnt_grid % level(lev) % n_cells ! - a % pnt_grid % comm % n_buff_cells

  res = 0.0

  !------------------------------!
  !   Diagonal preconditioning   !
  !------------------------------!
  do i = 1, ni
    d % val(d % dia(i)) = a % val(a % dia(i))
  end do

  !-----------------------------------!
  !    This is quite tricky point.    !
  !   What if bnrm2 is very small ?   !
  !-----------------------------------!
  bnrm2 = Normalized_Root_Mean_Square(ni, r1(1:ni), a, x(1:nt))

  if(bnrm2 < tol) then
    iter = 0
    goto 1
  end if

  !----------------!
  !   r = b - Ax   !
  !----------------!
  call Residual_Vector(ni, r1(1:ni), r1(1:ni), a, x(1:nt))

  !--------------------------------!
  !   Calculate initial residual   !
  !--------------------------------!
  res = Normalized_Root_Mean_Square(ni, r1(1:ni), a, x(1:nt))

  if(res < tol) then
    iter = 0
    goto 1
  end if

  !-----------!
  !   p = r   !
  !-----------!
  do i = 1, ni
    p1(i) = r1(i)
  end do

  !---------------!
  !               !
  !   Main loop   !
  !               !
  !---------------!
  do iter = 1, niter

    !----------------------!
    !     solve Mz = r     !
    !   (q instead of z)   !
    !----------------------!
    do i = 1, ni
      q1(i) = r1(i) / d % val(d % dia(i))
    end do

    !-----------------!
    !   rho = (r,z)   !
    !-----------------!
    rho = 0.
    do i = 1, ni
      rho = rho + r1(i)*q1(i)
    end do
    call Comm_Mod_Global_Sum_Real(rho)

    if(iter .eq. 1) then
      do i = 1, ni
        p1(i) = q1(i)
      end do
    else
      beta = rho/rho_old
      do i = 1, ni
        p1(i) = q1(i) + beta*p1(i)
      end do
    end if

    !------------!
    !   q = Ap   !
    !------------!
    call Comm_Mod_Exchange_Real(a % pnt_grid, p1)
    do i = 1, ni
      q1(i) = 0.
      do j = a % row(i), a % row(i+1)-1
        k = a % col(j)
        q1(i) = q1(i) + a % val(j) * p1(k)
      end do
    end do

    !------------------------!
    !   alfa = (r,z)/(p,q)   !
    !------------------------!
    alfa = 0.
    do i = 1, ni
      alfa = alfa + p1(i)*q1(i)
    end do
    call Comm_Mod_Global_Sum_Real(alfa)
    alfa = rho/alfa

    !---------------------!
    !   x = x + alfa p    !
    !   r = r - alfa Ap   !
    !---------------------!
    do i = 1, ni
      x(i)  = x(i)  + alfa*p1(i)
      r1(i) = r1(i) - alfa*q1(i)
    end do

    !-----------------------!
    !   Check convergence   !
    !-----------------------!
    res = Normalized_Root_Mean_Square(ni, r1(1:ni), a, x(1:nt))
    if(iter .eq. 1) then
      ini_res = res
    end if
PRINT *, 'ERROR = ', res

    if(res         < tol) goto 1
    if(res/ini_res < 0.1) goto 1

    rho_old = rho

  end do ! iter

1 continue
  fin_res = res
  niter = iter

  end subroutine
