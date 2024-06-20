!==============================================================================!
  subroutine Grad_Pressure(Flow, Grid, phi)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target, intent(inout) :: Flow  !! parent flow object
  type(Grid_Type),           intent(in)    :: Grid  !! grid object
  type(Var_Type)                           :: phi   !! pressure (correction)
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: phi_x(:), phi_y(:), phi_z(:)
  integer                   :: c, c1, c2, iter
  real                      :: dx, dy, dz
!==============================================================================!

  call Profiler % Start('Grad_Pressure')

  ! Store the name of the variable whose gradients you are computing
  Flow % stores_gradients_of = phi % name

  phi_x => Flow % phi_x
  phi_y => Flow % phi_y
  phi_z => Flow % phi_z

  !----------------------------------!
  !   Nullify arrays on the device   !
  !----------------------------------!

  !$acc parallel loop independent
  do c = Cells_At_Boundaries_In_Domain_And_Buffers()
    phi_x(c) = 0.0
    phi_y(c) = 0.0
    phi_z(c) = 0.0
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

    !$acc parallel loop independent           &
    !$acc present(grid_cells_c,               &
    !$acc         grid_xc, grid_yc, grid_zc,  &
    !$acc         phi_x,   phi_y,   phi_z)
    do c2 = Cells_At_Boundaries()
      c1 = grid_cells_c(1,c2)
      dx = grid_xc(c2) - grid_xc(c1)
      dy = grid_yc(c2) - grid_yc(c1)
      dz = grid_zc(c2) - grid_zc(c1)
      phi % n(c2) = phi % n(c1) + phi_x(c1) * dx  &
                                + phi_y(c1) * dy  &
                                + phi_z(c1) * dz
    end do
    !$acc end parallel

    !---------------------------------------------------------------!
    !   Compute pressure gradients again with extrapolated values   !
    !---------------------------------------------------------------!
    call Flow % Grad_Component(Grid, phi % n, 1, phi_x)
    call Flow % Grad_Component(Grid, phi % n, 2, phi_y)
    call Flow % Grad_Component(Grid, phi % n, 3, phi_z)

  end do

  call Profiler % Stop('Grad_Pressure')

  end subroutine
