!==============================================================================!
  subroutine Grad_Least_Variable(Flow, var)
!------------------------------------------------------------------------------!
!   Calculates gradient of a variable from field Flow with least squares       !
!                                                                              !
!   With OpenMP, this procedure got a speedup of 3.2 on 1M mesh and 4 threads. !
!   Yet, in the tests, it was called only to initialize Gaussian methods and   !
!   the measuring times where therefore short.  The numbers may be innacurate. !
!   It should be tested on a case with least squares only.                     !
!                                                                              !
!   On a case with least squares only, speed-up was 3.3, still not whopping.   !
!                                                                              !
!   But, the new version, the one which call all three components at once,     !
!   showed a speed-up of 3.7 on 1M mesh and 4 threads, which is good.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
  type(Var_Type)    :: var
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
!==============================================================================!

  call Profiler % Start('Grad_Least_Variable')

  ! Take alias
  Grid => Flow % pnt_grid

  ! Refresh buffers for variable
  call Grid % Exchange_Cells_Real(var % n)

  ! Old: Compute individual gradients without refreshing buffers
  ! call Flow % Grad_Component_No_Refresh(var % n, 1, var % x)  ! dp/dx
  ! call Flow % Grad_Component_No_Refresh(var % n, 2, var % y)  ! dp/dy
  ! call Flow % Grad_Component_No_Refresh(var % n, 3, var % z)  ! dp/dz

  ! New: Compute all gradient components without refreshing buffers
  call Flow % Grad_Three_Components_No_Refresh(var % n,  &
                                               var % x, var % y, var % z)

  ! Refresh buffers for gradient components
  call Grid % Exchange_Cells_Real(var % x)
  call Grid % Exchange_Cells_Real(var % y)
  call Grid % Exchange_Cells_Real(var % z)

  call Profiler % Stop('Grad_Least_Variable')

  end subroutine
