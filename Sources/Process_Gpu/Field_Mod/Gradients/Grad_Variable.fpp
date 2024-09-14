!==============================================================================!
  subroutine Grad_Variable(Flow, Grid, phi)
!------------------------------------------------------------------------------!
!   Essentially the same as Grad_Pressure, but without iterations              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target, intent(inout) :: Flow  !! parent flow object
  type(Grid_Type),   target, intent(in)    :: Grid  !! grid object
  type(Var_Type)                           :: phi   !! some variable
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: phi_x(:), phi_y(:), phi_z(:)
  integer                   :: c
!==============================================================================!

  call Profiler % Start('Grad_Variable')

  ! Store the name of the variable whose gradients you are computing
  Flow % stores_gradients_of = phi % name

  phi_x => Flow % phi_x
  phi_y => Flow % phi_y
  phi_z => Flow % phi_z

  !----------------------------------!
  !   Nullify arrays on the device   !
  !----------------------------------!

  !$tf-acc loop begin
  do c = Cells_At_Boundaries_In_Domain_And_Buffers()  ! all present
    phi_x(c) = 0.0
    phi_y(c) = 0.0
    phi_z(c) = 0.0
  end do
  !$tf-acc loop end

  !---------------------------------------------------------------!
  !   Compute pressure gradients again with extrapolated values   !
  !---------------------------------------------------------------!
  call Flow % Grad_Component(Grid, phi % n, 1, phi_x)
  call Flow % Grad_Component(Grid, phi % n, 2, phi_y)
  call Flow % Grad_Component(Grid, phi % n, 3, phi_z)

  call Profiler % Stop('Grad_Variable')

  end subroutine
