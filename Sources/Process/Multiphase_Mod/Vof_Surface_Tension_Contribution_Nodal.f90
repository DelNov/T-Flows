!==============================================================================!
  subroutine Multiphase_Mod_Vof_Surface_Tension_Contribution_Nodal(mult)
!------------------------------------------------------------------------------!
!    Computes the Curvature based on old-fashioned Brackbill's CSF approach    !
!    but using nodal values                                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: smooth_k   => r_cell_08,    &
                      smooth_n   => r_cell_09,    &
                      tmp_curv   => r_cell_10,    &
                      var_node_k => r_node_08,    &
                      var_node_n => r_node_09
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
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, mult % dist_func % n, &
                                            smooth_k, mult % n_conv_curv)
    else
      smooth_k = mult % dist_func % n
    end if
    call Multiphase_Mod_Vof_Interpolate_Cells_Nodes(grid,                     &
                                                    smooth_k, var_node_k)

    if (mult % n_conv_norm > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, mult % dist_func % n, &
                                            smooth_n, mult % n_conv_norm)
    else
      smooth_n = mult % dist_func % oo
    end if
    call Multiphase_Mod_Vof_Interpolate_Cells_Nodes(grid,                     &
                                                    smooth_n, var_node_n)

  else ! using VOF
    if (mult % n_conv_curv > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, vof % n, &
                                            smooth_k, mult % n_conv_curv)
    else
      smooth_k = vof % n
    end if
    call Multiphase_Mod_Vof_Interpolate_Cells_Nodes(grid,                     &
                                                    smooth_k, var_node_k)

    if (mult % n_conv_norm > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, vof % n, &
                                            smooth_n, mult % n_conv_norm)
    else
      smooth_n = vof % n
    end if
    call Multiphase_Mod_Vof_Interpolate_Cells_Nodes(grid,                     &
                                                    smooth_n, var_node_n)

  end if

  call Multiphase_Mod_Vof_Curvature_Nodal(grid, mult, smooth_k, var_node_k)

  vof % oo = mult % curv

  end subroutine
