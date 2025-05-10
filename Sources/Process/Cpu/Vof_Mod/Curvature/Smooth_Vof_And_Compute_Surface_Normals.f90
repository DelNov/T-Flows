!==============================================================================!
  subroutine Smooth_Vof_And_Compute_Surface_Normals(Vof)
!------------------------------------------------------------------------------!
!   Smooths VOF, computes its gradients and eventually surface normals         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  integer                   :: c, nb, nc
  real                      :: norm_grad
!==============================================================================!

  ! First take aliases
  Flow => Vof % pnt_flow
  Grid => Vof % pnt_grid

  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

  !---------------------------------------------------!
  !   Curvature smoothing was engaged                 !
  !---------------------------------------------------!
  if(Vof % n_smooth_for_curv_csf > 0) then

    ! Calculate smooth variable from vof function
    call Vof % Smooth_Scalar(Grid,                      &
                             Vof % fun    % n,          &
                             Vof % smooth % n(-nb:nc),  &
                             Vof % n_smooth_for_curv_csf)

  !-------------------------------------------------------!
  !   Curvature smoothing was not engaged                 !
  !-------------------------------------------------------!
  else

    ! Take smooth variable to be the same as VOF itself
    do c = Cells_At_Boundaries_In_Domain_And_Buffers()
      Vof % smooth % n(c) = Vof % fun % n(c)
    end do

  end if

  !----------------------------------------------------------!
  !   Compute gradients of the smooth function in any case   !
  !----------------------------------------------------------!
  call Flow % Grad_Variable(Vof % smooth)

  !-------------------------------------------!
  !   Once you have gradients of the smooth   !
  !    compute normals to the surface too     !
  !-------------------------------------------!
  do c = Cells_In_Domain()
    norm_grad = sqrt(  Vof % smooth % x(c) ** 2  &
                     + Vof % smooth % y(c) ** 2  &
                     + Vof % smooth % z(c) ** 2)
    if(norm_grad >= FEMTO) then
      Vof % nx(c) = Vof % smooth % x(c) / norm_grad
      Vof % ny(c) = Vof % smooth % y(c) / norm_grad
      Vof % nz(c) = Vof % smooth % z(c) / norm_grad
    else
      Vof % nx(c) = 0.0
      Vof % ny(c) = 0.0
      Vof % nz(c) = 0.0
    end if
  end do
  call Grid % Exchange_Cells_Real(Vof % nx(-nb:nc))
  call Grid % Exchange_Cells_Real(Vof % ny(-nb:nc))
  call Grid % Exchange_Cells_Real(Vof % nz(-nb:nc))

  end subroutine
