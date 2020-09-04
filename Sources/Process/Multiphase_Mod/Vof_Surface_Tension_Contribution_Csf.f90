!==============================================================================!
  subroutine Multiphase_Mod_Vof_Surface_Tension_Contribution_Csf(mult)
!------------------------------------------------------------------------------!
!    Computes the Curvature based on old-fashioned Brackbill's CSF approach    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: grad_nx    => r_cell_03,  & ! grad on x of vof for normal
                      grad_ny    => r_cell_04,  & ! grad on y of vof for normal
                      grad_nz    => r_cell_05,  & ! grad on z of vof for normal
                      grad_kx    => r_cell_06,  & ! grad on x of vof for curvat
                      grad_ky    => r_cell_07,  & ! grad on y of vof for curvat
                      grad_kz    => r_cell_08,  & ! grad on z of vof for curvat
                      smooth_var => r_cell_09
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),pointer :: flow
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: vof
  integer                  :: s, c, c1, c2, c_iter, nb, nc
  real                     :: vol_face, grad_face(3), grad_control(3)
  real                     :: dotprod, sxyz_mod, sxyz_control, fs, epsloc
  real                     :: d_n(3)     ! normal pointing to the wall
  real                     :: norm_grad  ! normal of a gradient
  real                     :: cpc
!==============================================================================!

  epsloc = epsilon(epsloc)

  ! First take aliases
  flow => mult % pnt_flow
  grid => mult % pnt_grid
  vof  => mult % vof

  nb = grid % n_bnd_cells
  nc = grid % n_cells

  cpc = 0.5
  if(mult % d_func) then  ! using distance function

    if (mult % n_conv_curv > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, mult % dist_func % oo, &
                                     smooth_var(-nb:nc), mult % n_conv_curv)

      call Field_Mod_Grad(flow, smooth_var(-nb:nc), grad_kx(-nb:nc),  &
                                                    grad_ky(-nb:nc),  &
                                                    grad_kz(-nb:nc))
    else
      call Field_Mod_Grad(flow, mult % dist_func % oo, grad_kx(-nb:nc),  &
                                                       grad_ky(-nb:nc),  &
                                                       grad_kz(-nb:nc))
    end if

    if (mult % n_conv_norm > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, mult % dist_func % oo, &
                                     smooth_var(-nb:nc), mult % n_conv_norm)

      call Field_Mod_Grad(flow, smooth_var(-nb:nc), grad_nx(-nb:nc),  &
                                                    grad_ny(-nb:nc),  &
                                                    grad_nz(-nb:nc))
      vof % x(-nb:nc) = grad_nx(-nb:nc)
      vof % y(-nb:nc) = grad_ny(-nb:nc)
      vof % z(-nb:nc) = grad_nz(-nb:nc)
    else
      call Field_Mod_Grad(flow, mult % dist_func % oo, vof % x,  &
                                                       vof % y,  &
                                                       vof % z)
    end if

  else ! using VOF

    if (mult % n_conv_curv > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, vof % n,   &
                                     smooth_var(-nb:nc), mult % n_conv_curv)

      call Field_Mod_Grad(flow, smooth_var(-nb:nc), grad_kx(-nb:nc),  &
                                                    grad_ky(-nb:nc),  &
                                                    grad_kz(-nb:nc))
    else
      grad_kx(-nb:nc) = vof % x(-nb:nc)
      grad_ky(-nb:nc) = vof % y(-nb:nc)
      grad_kz(-nb:nc) = vof % z(-nb:nc)
    end if

    if (mult % n_conv_norm > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, vof % n,   &
                                     smooth_var(-nb:nc), mult % n_conv_norm)

      call Field_Mod_Grad(flow, smooth_var(-nb:nc), grad_nx(-nb:nc),  &
                                                    grad_ny(-nb:nc),  &
                                                    grad_nz(-nb:nc))

      vof % x(-nb:nc) = grad_nx(-nb:nc)
      vof % y(-nb:nc) = grad_ny(-nb:nc)
      vof % z(-nb:nc) = grad_nz(-nb:nc)
    end if
  end if

  call Multiphase_Mod_Vof_Curvature_Csf(grid, mult,                  &
                grad_kx(-nb:nc), grad_ky(-nb:nc), grad_kz(-nb:nc),   &
                smooth_var(-nb:nc))

  end subroutine
