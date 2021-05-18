!==============================================================================!
  subroutine Smooth_For_Curvature_Csf(Vof)
!------------------------------------------------------------------------------!
!   Smooths vof function or distance function and computes its gradients       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: fun
  type(Var_Type),   pointer :: smooth
  integer                   :: c, nb, nc
!==============================================================================!

  ! First take aliases
  Flow   => Vof % pnt_flow
  Grid   => Vof % pnt_grid
  fun    => Vof % fun
  smooth => Vof % smooth

  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

  !---------------------------------------------------!
  !   Curvature smoothing was engaged                 !
  !---------------------------------------------------!
  if(Vof % n_conv_curv > 0) then

    ! Calculate smooth variable from vof function ...
    call Vof % Smooth_Scalar(Grid, fun % n,   &
                             smooth % n(-nb:nc), Vof % n_conv_curv)

    ! ... and find its gradients as well
    call Flow % Grad_Variable(smooth)

  !-------------------------------------------------------!
  !   Curvature smoothing was not engaged                 !
  !-------------------------------------------------------!
  else

    ! Take smooth variable to be the same as VOF itself ...
    do c = -Grid % n_bnd_cells, Grid % n_cells
      smooth % n(c) = fun % n(c)
    end do

    ! ... and find its gradients as well
    call Flow % Grad_Variable(smooth)
  end if

  !-------------------------------------------------------------!
  !   Smoothing of normal is engaged                            !
  !-------------------------------------------------------------!
  if(Vof % n_conv_norm > 0) then
    call Vof % Smooth_Scalar(Grid, fun % n,   &
                             smooth % n(-nb:nc), Vof % n_conv_norm)
    call Flow % Grad_Variable(smooth)
  end if

  end subroutine
