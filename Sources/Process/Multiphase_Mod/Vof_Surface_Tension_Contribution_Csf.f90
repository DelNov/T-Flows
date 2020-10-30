!==============================================================================!
  subroutine Multiphase_Mod_Vof_Surface_Tension_Contribution_Csf(mult)
!------------------------------------------------------------------------------!
!   Computes the Curvature based on Brackbill's CSF approach                   !
!   It actually smooths vof or distance function and computes its gradients    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: grad_kx    => r_cell_06,  & ! grad on x of vof for curvat
                      grad_ky    => r_cell_07,  & ! grad on y of vof for curvat
                      grad_kz    => r_cell_08
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: vof
  type(Var_Type),   pointer :: smooth
  integer                   :: s, c, c1, c2, c_iter, nb, nc
  real                      :: vol_face, grad_face(3), grad_control(3)
  real                      :: dotprod, sxyz_mod, sxyz_control, fs, epsloc
  real                      :: d_n(3)     ! normal pointing to the wall
  real                      :: norm_grad  ! normal of a gradient
!==============================================================================!

  epsloc = epsilon(epsloc)

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
    call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, vof % n,   &
                                   smooth % n(-nb:nc), mult % n_conv_curv)
    call Field_Mod_Grad_Variable(flow, smooth)

  !-------------------------------------------------------!
  !   Curvature smoothing was not engaged                 !
  !-------------------------------------------------------!
  else
    smooth % n(:) = vof % n(:)
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

  !------------------------------------------------------!
  !                                                      !
  !   No matter what you did above, compute curvature    !
  !   using the smoothed variable and gradients stored   !
  !   in vof % x, vof % y and vof % z                    !
  !                                                      !
  !------------------------------------------------------!
  call Multiphase_Mod_Vof_Curvature_Csf(mult)

  end subroutine
