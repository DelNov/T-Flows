!==============================================================================!
  subroutine Grad_Pressure(Flow, Grid, phi)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target, intent(inout) :: Flow  !! parent flow object
  type(Grid_Type),           intent(in)    :: Grid  !! grid object
  type(Var_Type),    target                :: phi   !! pressure (correction)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, iter, s, reg
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
    do reg = Boundary_Regions()

      if(Grid % region % type(reg) .ne. PRESSURE) then

        !$tf-acc loop begin
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          phi % n(c2) = phi % n(c1) + Flow % phi_x(c1) * Grid % dx(s)  &
                                    + Flow % phi_y(c1) * Grid % dy(s)  &
                                    + Flow % phi_z(c1) * Grid % dz(s)
        end do
        !$tf-acc loop end

      else

        !$tf-acc loop begin
        do s = Faces_In_Region(reg)
          c2 = Grid % faces_c(2,s)
          phi % n(c2) = 0.0
        end do
        !$tf-acc loop end

      end if  ! pressure or not

    end do    ! regions

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
