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
  real, contiguous, pointer :: phi_x(:), phi_y(:), phi_z(:)
  integer                   :: c, c1, c2, reg, s
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

  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   phi_x,  &
  !$acc   phi_y,  &
  !$acc   phi_z   &
  !$acc )
  do c = grid_region_f_cell(1), grid_region_l_cell(grid_n_regions+1)  ! all present
    phi_x(c) = 0.0
    phi_y(c) = 0.0
    phi_z(c) = 0.0
  end do
  !$acc end parallel

  !----------------------------------------!
  !   Copy values to symmetry boundaries   !
  !   (Probably not the most consistent)   !
  !----------------------------------------!
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. SYMMETRY) then

      phi_n => phi % n
      !$acc parallel loop  &
      !$acc present(  &
      !$acc   grid_region_f_face,  &
      !$acc   grid_region_l_face,  &
      !$acc   grid_faces_c,  &
      !$acc   phi_n   &
      !$acc )
      do s = grid_region_f_face(reg), grid_region_l_face(reg)  ! all present
        c1 = grid_faces_c(1,s)   ! inside cell
        c2 = grid_faces_c(2,s)   ! boundary cell
        phi_n(c2) = phi_n(c1)
      end do
      !$acc end parallel

    end if
  end do

  !---------------------------------------------------------------!
  !   Compute pressure gradients again with extrapolated values   !
  !---------------------------------------------------------------!
  call Flow % Grad_Component(Grid, phi % n, 1, phi_x, boundary_updated=.true.)
  call Flow % Grad_Component(Grid, phi % n, 2, phi_y, boundary_updated=.true.)
  call Flow % Grad_Component(Grid, phi % n, 3, phi_z, boundary_updated=.true.)

  call Profiler % Stop('Grad_Variable')

  end subroutine
