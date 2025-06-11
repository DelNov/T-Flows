!==============================================================================!
  subroutine Solve_On_Coarsest_Level(Amg, m, ifac,     &
                                     a, u, f, ia, ja,  &
                                     iw, icg)
!------------------------------------------------------------------------------!
!   Solves on coarsest grid, either with Gauss-Seidel relaxation, Conjugate
!   Gradient or Bi-Cojugate Gradient solvers.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: nifac
  real            :: a(:), u(:), f(:)
  integer         :: ia(:), ja(:)
  integer         :: iw(:), icg(:)
!-----------------------------------[locals]-----------------------------------!
  real    :: fmax, resnew, resold
  integer :: m,ifac,aaux,esp,flag,i,iaux,ihi,ii,ilo,is, &
             iter,j,jhi,jj,jlo,jpos,js,np,npoint,nsp,path

  !------------------------------------------------------------------------!
  !   conv: if coarse grid solution is done with gs-relaxation and         !
  !   n_relax_coarse=0, as many gs-sweeps are performed as are necessary   !
  !   to reduce the residual by the factor conv                            !
  !------------------------------------------------------------------------!
  real, parameter :: conv =1.0e-2
!------------------------------------[save]------------------------------------!
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
        call Amg % Gauss_Seidel_Sweep(m, 2,             &
                                      a, u, f, ia, ja,  &
                                      iw, icg)
      end do

    !-----------------------------------------------------!
    !   Reduce residual on coarsest grid by factor conv   !
    !   if not yet in the range of the truncation error   !
    !-----------------------------------------------------!
    else

      call Amg % timer_start()

     ! Calculate supremum norm of right hand side
      fmax = 0.0
      do i = Amg % imin(m), Amg % imax(m)
        fmax = max(fmax,abs(f(i)))
      end do
      call Amg % Calculate_Residual(m, resold,  &
                                    a, u, f, ia, ja,  &
                                    iw)

      call Amg % timer_stop(15)
      resold = max(resold*conv,fmax*1.0e-12)
      do i = 1, 10
        do j = 1, 10
          call Amg % Gauss_Seidel_Sweep(m, 2, a, u, f, ia, ja,  &
                                        iw, icg)
        end do
        call Amg % timer_start()
        call Amg % Calculate_Residual(m, resnew, a, u, f, ia, ja,  &
                                      iw)

        call Amg % timer_stop(15)
        if(resnew .le. resold) return
      end do
    end if

  !-----------------------------!
  !                             !
  !   Solution with CG solver   !
  !                             !
  !-----------------------------!
  else if(Amg % coarse_solver .eq. AMG_SOLVER_CG) then

    call Amg % Cg_On_Level(m, 2,             &
                           a, u, f, ia, ja,  &
                           iw, icg)

  !-------------------------------!
  !                               !
  !   Solution with BiCG solver   !
  !                               !
  !-------------------------------!
  else if(Amg % coarse_solver .eq. AMG_SOLVER_BICG) then

    call Amg % Bicg_On_Level(m, 2,             &
                             a, u, f, ia, ja,  &
                             iw,icg)

  end if

  end subroutine
