!==============================================================================!
  subroutine Cg_Level(lev, a, d, x, b, niter, tol, ini_res, fin_res)
!------------------------------------------------------------------------------!
!   Conjugate gradient method for one level of the multigrid.                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Matrix_Mod
  use Work_Mod, only: p1 => r_cell_01,  &
                      q1 => r_cell_02,  &
                      r1 => r_cell_03
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: lev
  type(Matrix_Type) :: a
  type(Matrix_Type) :: d
  real              :: x (a % pnt_grid % level(lev) % n_cells)
  real              :: b (a % pnt_grid % level(lev) % n_cells)  ! [A]{x}={b}
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
  d % val(d % dia(1:ni)) = a % val(a % dia(1:ni))

  !-----------------------------------!
  !    This is quite tricky point.    !
  !   What if bnrm2 is very small ?   !
  !-----------------------------------!
  bnrm2 = Normalized_Root_Mean_Square(ni, b(1:nt), a, x(1:nt))
PRINT *, ' INITIAL BNRM2 = ', bnrm2

  if(bnrm2 < tol) then
    iter = 0
    goto 1
  end if

  !----------------!
  !   r = b - Ax   !
  !----------------!
  call Residual_Vector(ni, r1(1:nt), b(1:nt), a, x(1:nt))

  !--------------------------------!
  !   Calculate initial residual   !
  !--------------------------------!
  res = Normalized_Root_Mean_Square(ni, r1(1:nt), a, x(1:nt))
PRINT *, ' INITIAL ERROR = ', res

  if(res < tol) then
    iter = 0
    goto 1
  end if

  !-----------!
  !   p = r   !
  !-----------!
  p1(1:ni) = r1(1:ni)

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
    q1(1:ni) = r1(1:ni) / d % val(d % dia(1:ni))

    !-----------------!
    !   rho = (r,z)   !
    !-----------------!
    rho = dot_product(r1(1:ni), q1(1:ni))
    call Comm_Mod_Global_Sum_Real(rho)

    if(iter .eq. 1) then
      p1(1:ni) = q1(1:ni)
    else
      beta = rho / rho_old
      p1(1:ni) = q1(1:ni) + beta * p1(1:ni)
    end if

    !------------!
    !   q = Ap   !
    !------------!
    call Comm_Mod_Exchange_Real(a % pnt_grid, p1)
    do i = 1, ni
      q1(i) = 0.0
      do j = a % row(i), a % row(i+1)-1
        k = a % col(j)
        q1(i) = q1(i) + a % val(j) * p1(k)
      end do
    end do

    !------------------------!
    !   alfa = (r,z)/(p,q)   !
    !------------------------!
    alfa = dot_product(p1(1:ni), q1(1:ni))
    call Comm_Mod_Global_Sum_Real(alfa)
    alfa = rho/alfa

    !---------------------!
    !   x = x + alfa p    !
    !   r = r - alfa Ap   !
    !---------------------!
    x (1:ni) = x (1:ni) + alfa * p1(1:ni)
    r1(1:ni) = r1(1:ni) - alfa * q1(1:ni)

    !-----------------------!
    !   Check convergence   !
    !-----------------------!
    res = Normalized_Root_Mean_Square(ni, r1(1:nt), a, x(1:nt))
    if(iter .eq. 1) then
      ini_res = res
    end if
PRINT *, ' CURRENT ERROR = ', res

    if(res         < tol) goto 1
    if(res/ini_res < 0.1) goto 1

    rho_old = rho

  end do ! iter

1 continue
  fin_res = res
  niter = iter

  end subroutine
