!==============================================================================!
  subroutine Grad_Variable(Flow, Grid, phi)
!------------------------------------------------------------------------------!
!   Essentially the same as Grad_Pressure, but without iterations              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), intent(inout) :: Flow  !! parent flow object
  type(Grid_Type),   intent(in)    :: Grid  !! grid object
  type(Var_Type)                   :: phi   !! pressure or pressure correction
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, iter
  real    :: dx, dy, dz
!==============================================================================!

  call Profiler % Start('Grad_Variable')

  !----------------------------------!
  !   Nullify arrays on the device   !
  !----------------------------------!

  !$acc parallel loop independent
  do c = -Grid % n_bnd_cells, Grid % n_cells
    phi % x(c) = 0.0
    phi % y(c) = 0.0
    phi % z(c) = 0.0
  end do
  !$acc end parallel

  !---------------------------------------------------------------!
  !   Compute pressure gradients again with extrapolated values   !
  !---------------------------------------------------------------!
  call Flow % Grad_Component(Grid, phi % n, 1, phi % x)
  call Flow % Grad_Component(Grid, phi % n, 2, phi % y)
  call Flow % Grad_Component(Grid, phi % n, 3, phi % z)

  call Profiler % Stop('Grad_Variable')

  end subroutine
