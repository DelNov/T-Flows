!==============================================================================!
  subroutine Grad_Variable(Flow, Grid, phi)
!------------------------------------------------------------------------------!
!   Essentially the same as Grad_Pressure, but without iterations              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target, intent(inout) :: Flow  !! parent flow object
  type(Grid_Type),   target, intent(in)    :: Grid  !! grid object
  type(Var_Type),    target, intent(in)    :: phi   !! some variable
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, reg, s
!==============================================================================!

  call Profiler % Start('Grad_Variable')

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

  !----------------------------------------!
  !   Copy values to symmetry boundaries   !
  !   (Probably not the most consistent)   !
  !----------------------------------------!
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. SYMMETRY) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)  ! all present
        c1 = Grid % faces_c(1,s)   ! inside cell
        c2 = Grid % faces_c(2,s)   ! boundary cell
        phi % n(c2) = phi % n(c1)
      end do
      !$tf-acc loop end

    end if
  end do

  !---------------------------------------------------------------!
  !   Compute variable gradients again with extrapolated values   !
  !---------------------------------------------------------------!
  call Flow % Grad_Component(Grid, phi % n, 1,  &
                             Flow % phi_x, boundary_updated=.true.)
  call Flow % Grad_Component(Grid, phi % n, 2,  &
                             Flow % phi_y, boundary_updated=.true.)
  call Flow % Grad_Component(Grid, phi % n, 3,  &
                             Flow % phi_z, boundary_updated=.true.)

  call Profiler % Stop('Grad_Variable')

  end subroutine
