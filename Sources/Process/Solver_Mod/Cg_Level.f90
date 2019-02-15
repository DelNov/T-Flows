!==============================================================================!
  subroutine Cg_Level(lev, a, d, x, b, p, r, prec,    &
                           cyc, max_iter, act_iter,         &
                           tol, res_rat, fin_res)
!------------------------------------------------------------------------------!
!   Conjugate gradient method for one level of the multigrid.                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: z  => r_cell_02,  &
                      r1 => r_cell_03
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: lev
  type(Matrix_Type) :: a
  type(Matrix_Type) :: d
  real              :: x(a % pnt_grid % level(lev) % n_cells)
  real              :: b(a % pnt_grid % level(lev) % n_cells)  ! [A]{x}={b}
  real              :: p(a % pnt_grid % level(lev) % n_cells)
  real              :: r(a % pnt_grid % level(lev) % n_cells)  ! {r}={b}-[A]{x}
  character(len=80) :: prec                     ! preconditioner
  integer           :: cyc                      ! current cycle
  integer           :: max_iter                 ! maximum number of iterations
  integer           :: act_iter                 ! actual number of iterations
  real              :: tol                      ! tolerance
  real              :: res_rat                  ! residual ratio
  real              :: fin_res                  ! final residual
!-----------------------------------[Locals]-----------------------------------!
  integer :: nt, ni, i, j, k, iter
  real    :: alfa, beta, rho, rho_old, bnrm2, ini_res, res
!==============================================================================!

  ! Take some aliases
  nt = a % pnt_grid % level(lev) % n_cells
  ni = a % pnt_grid % level(lev) % n_cells ! - a % pnt_grid % comm % n_buff_cells

  res = 0.0

  !---------------------!
  !   Preconditioning   !
  !---------------------!
  call Prec_Form(ni, a, d, prec)

  !-----------------------------------!
  !    This is quite tricky point.    !
  !   What if bnrm2 is very small ?   !
  !-----------------------------------!
  bnrm2 = Root_Mean_Square(ni, b(1:nt))
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
  res = Root_Mean_Square(ni, r1(1:nt))
  ini_res = res
  PRINT '(a,i2.2,a,i3.3,a,es12.3)', ' LEVEL_',    lev,   &
                                    '; ITER_',    0,     &
                                    '; ERROR = ', res

  if(res < tol) then
    iter = 0
    goto 1
  end if

  !---------------!
  !               !
  !   Main loop   !
  !               !
  !---------------!
  do iter = 1, max_iter

    !----------------------!
    !     solve Mz = r     !
    !   (q instead of z)   !
    !----------------------!
    call Prec_Solve(ni, a, d, z(1:nt), r1(1:nt), prec)

    !-----------------!
    !   rho = (r,z)   !
    !-----------------!
    rho = dot_product(r1(1:ni), z(1:ni))
    call Comm_Mod_Global_Sum_Real(rho)

    if(iter .eq. 1 .and. cyc .eq. 1) then
      p(1:ni) = z(1:ni)
    else
      beta = rho / rho_old
      p(1:ni) = z(1:ni) + beta * p(1:ni)
    end if

    !------------!
    !   q = Ap   !
    !------------!
    call Comm_Mod_Exchange_Real(a % pnt_grid, p)
    do i = 1, ni
      z(i) = 0.0
      do j = a % row(i), a % row(i+1)-1
        k = a % col(j)
        z(i) = z(i) + a % val(j) * p(k)
      end do
    end do

    !------------------------!
    !   alfa = (r,z)/(p,q)   !
    !------------------------!
    alfa = dot_product(p(1:ni), z(1:ni))
    call Comm_Mod_Global_Sum_Real(alfa)
    alfa = rho/alfa

    !---------------------!
    !   x = x + alfa p    !
    !   r = r - alfa Ap   !
    !---------------------!
    x (1:ni) = x (1:ni) + alfa * p(1:ni)
    r1(1:ni) = r1(1:ni) - alfa * z(1:ni)

    !-----------------------!
    !   Check convergence   !
    !-----------------------!
    res = Root_Mean_Square(ni, r1(1:nt))
    PRINT '(a,i2.2,a,i3.3,a,es12.3)', ' LEVEL_',    lev,   &
                                      '; ITER_',     iter,  &
                                      '; ERROR = ', res

    if(res         < tol)     goto 1
    if(res/ini_res < res_rat) goto 1

    rho_old = rho

  end do ! iter

  !----------------------------------!
  !                                  !
  !   Convergence has been reached   !
  !                                  !
  !----------------------------------!
1 continue
  act_iter = iter
  fin_res  = res

  !---------------------------------!
  !   Compute: {r} = {b} - [A]{x}   !
  !---------------------------------!
  call Residual_Vector(ni, r(1:nt), b(1:nt), a, x(1:nt))

  end subroutine
