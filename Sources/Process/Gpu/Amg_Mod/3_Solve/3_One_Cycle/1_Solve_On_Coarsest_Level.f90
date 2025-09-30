!==============================================================================!
  subroutine Solve_On_Coarsest_Level(Amg, m, ifac)
!------------------------------------------------------------------------------!
!   Solves on coarsest grid, either with Gauss-Seidel relaxation, Conjugate
!   Gradient or Bi-Cojugate Gradient solvers.
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Arguments]----------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: m, ifac
!-----------------------------------[Locals]-----------------------------------!
  real    :: fmax, resnew, resold
  integer :: n,aaux,esp,flag,i,ihi,ii,ilo,is, &
             iter,j,jhi,jj,jlo,jpos,js,np,npoint,nsp,path
  real, contiguous, pointer :: f(:)

  !------------------------------------------------------------------------!
  !   conv: if coarse grid solution is done with gs-relaxation and         !
  !   n_relax_coarse=0, as many gs-sweeps are performed as are necessary   !
  !   to reduce the residual by the factor conv                            !
  !------------------------------------------------------------------------!
  real, parameter :: conv =1.0e-2
!------------------------------------[Save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  !-------------------------------------------!
  !                                           !
  !   Solution with gauss-seidel relaxation   !
  !                                           !
  !-------------------------------------------!
  if(Amg % coarse_solver .eq. AMG_SOLVER_GS) then

    !---------------------------------------------!
    !   Just perform given number of iterations   !
    !---------------------------------------------!
    if(Amg % n_relax_coarse .ne. 0) then
      do iter = 1, Amg % n_relax_coarse
        call Amg % Gauss_Seidel_Sweep(m, 2)
      end do

    !-----------------------------------------------------!
    !   Reduce residual on coarsest grid by factor conv   !
    !   if not yet in the range of the truncation error   !
    !-----------------------------------------------------!
    else

     ! Calculate supremum norm of right hand side
      fmax = 0.0
      n =  Amg % lev(m) % n
      f => Amg % lev(m) % f
      do i = 1, n
        fmax = max(fmax, abs(f(i)))
      end do
      call Amg % Calculate_Residual(m, resold)

      resold = max(resold*conv,fmax*1.0e-12)
      do i = 1, 10
        do j = 1, 10
          call Amg % Gauss_Seidel_Sweep(m, 2)
        end do
        call Amg % Calculate_Residual(m, resnew)

        if(resnew .le. resold) return
      end do
    end if

  !-----------------------------!
  !                             !
  !   Solution with CG solver   !
  !                             !
  !-----------------------------!
  else if(Amg % coarse_solver .eq. AMG_SOLVER_CG) then

    call Amg % Cg_On_Level(m, 2)

  !-------------------------------!
  !                               !
  !   Solution with BiCG solver   !
  !                               !
  !-------------------------------!
  else if(Amg % coarse_solver .eq. AMG_SOLVER_BICG) then

    call Amg % Bicg_On_Level(m, 2)

  end if

  end subroutine
