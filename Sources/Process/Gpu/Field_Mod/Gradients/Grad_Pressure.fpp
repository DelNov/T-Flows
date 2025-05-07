!==============================================================================!
  subroutine Grad_Pressure(Flow, Grid, phi)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target, intent(inout) :: Flow  !! parent flow object
  type(Grid_Type),           intent(in)    :: Grid  !! grid object
  type(Var_Type),    target                :: phi   !! pressure (correction)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, iter
  real    :: dx, dy, dz
!==============================================================================!

  call Profiler % Start('Grad_Pressure')

  ! Store the name of the variable whose gradients you are computing
  Flow % stores_gradients_of = phi % name

  !----------------------------------!
  !   Nullify arrays on the device   !
  !----------------------------------!

  !$tf-acc loop begin
  do c = Cells_At_Boundaries_In_Domain_And_Buffers()  ! all present
    Flow % phi_x(c) = 0.0
    Flow % phi_y(c) = 0.0
    Flow % phi_z(c) = 0.0
  end do
  !$tf-acc loop end

  !------------------------------------!
  !                                    !
  !   Iterativelly improve gradients   !
  !                                    !
  !------------------------------------!
  do iter = 1, Flow % least_miter

    !--------------------------------------!
    !   Extrapolate values to boundaries   !
    !--------------------------------------!

    !$tf-acc loop begin
    do c2 = Cells_At_Boundaries()  ! all present
      c1 = Grid % cells_c(1,c2)
      dx = Grid % xc(c2) - Grid % xc(c1)
      dy = Grid % yc(c2) - Grid % yc(c1)
      dz = Grid % zc(c2) - Grid % zc(c1)
      phi % n(c2) = phi % n(c1) + Flow % phi_x(c1) * dx  &
                                + Flow % phi_y(c1) * dy  &
                                + Flow % phi_z(c1) * dz
    end do
    !$tf-acc loop end

    !---------------------------------------------------------------!
    !   Compute pressure gradients again with extrapolated values   !
    !---------------------------------------------------------------!
    call Flow % Grad_Component(Grid, phi % n, 1,  &
                               Flow % phi_x, boundary_updated=.true.)
    call Flow % Grad_Component(Grid, phi % n, 2,  &
                               Flow % phi_y, boundary_updated=.true.)
    call Flow % Grad_Component(Grid, phi % n, 3,  &
                               Flow % phi_z, boundary_updated=.true.)

  end do

  call Profiler % Stop('Grad_Pressure')

  end subroutine
