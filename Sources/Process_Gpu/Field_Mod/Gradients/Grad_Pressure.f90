!==============================================================================!
  subroutine Grad_Pressure(Flow, Grid, phi, phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), intent(inout) :: Flow  !! parent flow object
  type(Grid_Type),   intent(in)    :: Grid  !! grid object
  type(Var_Type)                   :: phi
  real,    optional, intent(out)   :: phi_x(-Grid % n_bnd_cells:Grid % n_cells)
  real,    optional, intent(out)   :: phi_y(-Grid % n_bnd_cells:Grid % n_cells)
  real,    optional, intent(out)   :: phi_z(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, iter
  real    :: dx, dy, dz
!==============================================================================!

  call Profiler % Start('Grad_Pressure')

  ! All or nothing must be present
  Assert(present(phi_x) .eqv. present(phi_y))
  Assert(present(phi_y) .eqv. present(phi_z))

  !----------------------------------!
  !   Nullify arrays on the device   !
  !----------------------------------!

  ! External arrays are present, use them
  if(present(phi_x)) then
    !$acc parallel loop independent
    do c = -Grid % n_bnd_cells, Grid % n_cells
      phi_x(c) = 0.0
      phi_y(c) = 0.0
      phi_z(c) = 0.0
    end do
    !$acc end parallel

  ! External arrays are not present, use internal gradient components
  else
    !$acc parallel loop independent
    do c = -Grid % n_bnd_cells, Grid % n_cells
      phi % x(c) = 0.0
      phi % y(c) = 0.0
      phi % z(c) = 0.0
    end do
    !$acc end parallel

  end if

  !------------------------------------!
  !                                    !
  !   Iterativelly improve gradients   !
  !                                    !
  !------------------------------------!
  do iter = 1, 4

    !--------------------------------------!
    !   Extrapolate values to boundaries   !
    !--------------------------------------!

    ! External arrays are present, use them
    if(present(phi_x)) then
      !$acc parallel loop independent
      do c2 = -Grid % n_bnd_cells, -1
        c1 = Grid % cells_c(1,c2)
        dx = Grid % xc(c2) - Grid % xc(c1)
        dy = Grid % yc(c2) - Grid % yc(c1)
        dz = Grid % zc(c2) - Grid % zc(c1)
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

    ! External arrays are not present, use internal gradient components
    else
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
    end if

  end do

  call Profiler % Stop('Grad_Pressure')

  end subroutine
