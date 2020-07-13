!==============================================================================!
  subroutine Multiphase_Mod_Vof_Surface_Tension_Contribution(mult)
!------------------------------------------------------------------------------!
!    Computes the Curvature based on old-fashioned Brackbill's CSF approach    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: grad_nx   => r_cell_03,  & !grad on x of vof for normal
                      grad_ny   => r_cell_04,  & !grad on y of vof for normal
                      grad_nz   => r_cell_05,  & !grad on z of vof for normal
                      grad_kx   => r_cell_06,  & !grad on x of vof for curvat
                      grad_ky   => r_cell_07,  & !grad on y of vof for curvat
                      grad_kz   => r_cell_08,  & !grad on z of vof for curvat
                      tmp_curv  => r_cell_30,  & !curvature holder
                      smooth_var=> r_cell_09
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),pointer :: flow
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: vof
  integer                  :: s, c, c1, c2, c_iter
  real                     :: vol_face, grad_face(3), grad_control(3)
  real                     :: dotprod, sxyz_mod, sxyz_control, fs, epsloc
  real                     :: d_n(3) !normal pointing to the wall
  real                     :: norm_grad !normal of a gradient
  real                     :: cpc
!==============================================================================!

  epsloc = epsilon(epsloc)

  ! First take aliases
  flow => mult % pnt_flow
  grid => mult % pnt_grid
  vof  => mult % vof

  cpc = 0.5
  if(mult % d_func) then  ! using distance function

    if (mult % n_conv_curv > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, mult % dist_func % oo, &
                                            smooth_var, mult % n_conv_curv)

      call Field_Mod_Grad_Component(flow, smooth_var, 1, grad_kx)
      call Field_Mod_Grad_Component(flow, smooth_var, 2, grad_ky)
      call Field_Mod_Grad_Component(flow, smooth_var, 3, grad_kz)
    else
      call Field_Mod_Grad_Component(flow, mult % dist_func % oo, 1, grad_kx)
      call Field_Mod_Grad_Component(flow, mult % dist_func % oo, 2, grad_ky)
      call Field_Mod_Grad_Component(flow, mult % dist_func % oo, 3, grad_kz)
    end if

    if (mult % n_conv_norm > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, mult % dist_func % oo, &
                                            smooth_var, mult % n_conv_norm)

      call Field_Mod_Grad_Component(flow, smooth_var, 1, grad_nx)
      call Field_Mod_Grad_Component(flow, smooth_var, 2, grad_ny)
      call Field_Mod_Grad_Component(flow, smooth_var, 3, grad_nz)
      vof % x = grad_nx
      vof % y = grad_ny
      vof % z = grad_nz
    else
      call Field_Mod_Grad_Component(flow, mult % dist_func % oo, 1, vof % x)
      call Field_Mod_Grad_Component(flow, mult % dist_func % oo, 2, vof % y)
      call Field_Mod_Grad_Component(flow, mult % dist_func % oo, 3, vof % z)
    end if

  else ! using VOF

    if (mult % n_conv_curv > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, vof % n,   &
                                           smooth_var, mult % n_conv_curv)

      call Field_Mod_Grad_Component(flow, smooth_var, 1, grad_kx)
      call Field_Mod_Grad_Component(flow, smooth_var, 2, grad_ky)
      call Field_Mod_Grad_Component(flow, smooth_var, 3, grad_kz)
    else
      grad_kx = vof % x
      grad_ky = vof % y
      grad_kz = vof % z
    end if

    if (mult % n_conv_norm > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, vof % n,   &
                                            smooth_var, mult % n_conv_norm)

      call Field_Mod_Grad_Component(flow, smooth_var, 1, grad_nx)
      call Field_Mod_Grad_Component(flow, smooth_var, 2, grad_ny)
      call Field_Mod_Grad_Component(flow, smooth_var, 3, grad_nz)

      vof % x = grad_nx
      vof % y = grad_ny
      vof % z = grad_nz
    end if
  end if

  call Multiphase_Mod_Vof_Curvature_Csf(grid, mult,                  &
                                        grad_kx, grad_ky, grad_kz,   &
                                        smooth_var)

  vof % oo = mult % curv

  end subroutine
