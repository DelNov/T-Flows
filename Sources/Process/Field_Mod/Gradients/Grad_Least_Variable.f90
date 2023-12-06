!==============================================================================!
  subroutine Grad_Least_Variable(Flow, var)
!------------------------------------------------------------------------------!
!>  This subroutine calculates the gradient of a generic variable using the
!>  least squares method.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Buffer refresh: Initiates by refreshing the buffers for the target       !
!     variable (var % n) to ensure up-to-date data across all processors.      !
!   * Gradient calculation: Uses a least squares method to compute the         !
!     gradient components (var % x, var % y, and var % z) of the variable.     !
!     This method is executed without refreshing buffers between component     !
!     calculations, enhancing computational efficiency.                        !
!   * Buffer refresh for gradients: Post gradient computation, the buffers     !
!     for each gradient component are refreshed to maintain data consistency   !
!     across the computational domain.                                         !
!------------------------------------------------------------------------------!
!   Notes                                                                      !
!                                                                              !
!   * The procedure has shown a moderate speedup with OpenMP on a 1M mesh and  !
!     4 threads. Although initial tests indicated a speedup of 3.2, further    !
!     tests with least squares only demonstrated a speedup of 3.3. A refined   !
!     version executing all three components at once achieved a speedup of     !
!     3.7, highlighting its efficiency in a OMP parallel environment.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow  !! parent flow object
  type(Var_Type)    :: var   !! variable whose gradients are calculated
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
!==============================================================================!

  call Profiler % Start('Grad_Least_Variable')

  ! Take alias
  Grid => Flow % pnt_grid

  ! Refresh buffers for variable
  call Grid % Exchange_Cells_Real(var % n)

  ! New: Compute all gradient components without refreshing buffers
  call Flow % Grad_Three_Components_No_Refresh(Grid, var % n,  &
                                               var % x, var % y, var % z)

  ! Refresh buffers for gradient components
  call Grid % Exchange_Cells_Real(var % x)
  call Grid % Exchange_Cells_Real(var % y)
  call Grid % Exchange_Cells_Real(var % z)

  call Profiler % Stop('Grad_Least_Variable')

  end subroutine
