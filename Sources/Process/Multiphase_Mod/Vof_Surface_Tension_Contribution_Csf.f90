!==============================================================================!
  subroutine Multiphase_Mod_Vof_Surface_Tension_Contribution_Csf(mult)
!------------------------------------------------------------------------------!
!   Computes the Curvature based on Brackbill's CSF approach                   !
!   It actually smooths vof or distance function and computes its gradients    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: vof
  type(Var_Type),   pointer :: smooth
  integer                   :: c, nb, nc
!==============================================================================!

  ! First take aliases
  flow   => mult % pnt_flow
  grid   => mult % pnt_grid
  vof    => mult % vof
  smooth => mult % smooth

  nb = grid % n_bnd_cells
  nc = grid % n_cells

  !---------------------------------------------------!
  !   Curvature smoothing was engaged                 !
  !---------------------------------------------------!
  if (mult % n_conv_curv > 0) then

    ! Calculate smooth variable from vof ...
    call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, vof % n,   &
                                   smooth % n(-nb:nc), mult % n_conv_curv)

    ! ... and find its gradients as well
    call Field_Mod_Grad_Variable(flow, smooth)

  !-------------------------------------------------------!
  !   Curvature smoothing was not engaged                 !
  !-------------------------------------------------------!
  else

    ! Take smooth variable to be the same as VOF itself ...
    do c = -grid % n_bnd_cells, grid % n_cells
      smooth % n(c) = vof % n(c)
    end do

    ! ... and find its gradients as well
    call Field_Mod_Grad_Variable(flow, smooth)
  end if

  !-------------------------------------------------------------!
  !   Smoothing of normal is engaged                            !
  !-------------------------------------------------------------!
  if (mult % n_conv_norm > 0) then
    call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, vof % n,   &
                                   smooth % n(-nb:nc), mult % n_conv_norm)
    call Field_Mod_Grad_Variable(flow, smooth)
  end if

  !----------------------------------------------------------!
  !                                                          !
  !   No matter what you did above, compute curvature next   !
  !                                                          !
  !----------------------------------------------------------!
  call Multiphase_Mod_Vof_Curvature_Csf(mult)

  end subroutine
