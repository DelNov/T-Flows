!==============================================================================!
  subroutine Grad_Pressure(Flow, Grid, phi)
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

  call Profiler % Start('Grad_Pressure')

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

  !------------------------------------!
  !                                    !
  !   Iterativelly improve gradients   !
  !                                    !
  !------------------------------------!
  do iter = 1, Flow % least_miter

    !--------------------------------------!
    !   Extrapolate values to boundaries   !
    !--------------------------------------!

    !$acc parallel loop independent
    do c2 = -Grid % n_bnd_cells, -1
      c1 = Grid % cells_c(1,c2)
      dx = Grid % xc(c2) - Grid % xc(c1)
      dy = Grid % yc(c2) - Grid % yc(c1)
      dz = Grid % zc(c2) - Grid % zc(c1)
      phi % n(c2) = phi % n(c1) + phi % x(c1) * dx  &
                                + phi % y(c1) * dy  &
                                + phi % z(c1) * dz
    end do
    !$acc end parallel

    !---------------------------------------------------------------!
    !   Compute pressure gradients again with extrapolated values   !
    !---------------------------------------------------------------!
    call Flow % Grad_Component(Grid, phi % n, 1, phi % x)
    call Flow % Grad_Component(Grid, phi % n, 2, phi % y)
    call Flow % Grad_Component(Grid, phi % n, 3, phi % z)

  end do

  call Profiler % Stop('Grad_Pressure')

  end subroutine
