!==============================================================================!
  subroutine Multiphase_Mod_Vof_Surface_Tension_Contribution_Nodal(mult)
!------------------------------------------------------------------------------!
!    Computes the Curvature based on old-fashioned Brackbill's CSF approach    !
!    but using nodal values                                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: smooth_k   => r_cell_08,    &
                      smooth_n   => r_cell_09,    &
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
  integer                  :: s, c, c1, c2, c_iter, nb, nc, nn
  real                     :: vol_face, grad_face(3), grad_control(3)
  real                     :: dotprod, sxyz_mod, sxyz_control, fs, epsloc
  real                     :: d_n(3)     ! normal pointing to the wall
  real                     :: norm_grad  ! normal of a gradient
!==============================================================================!

  epsloc = epsilon(epsloc)

  ! First take aliases
  flow => mult % pnt_flow
  grid => mult % pnt_grid
  vof  => mult % vof

  nb = grid % n_bnd_cells
  nc = grid % n_cells
  nn = grid % n_nodes

  !-------------------------------!
  !                               !
  !   Distance function is used   !
  !                               !
  !-------------------------------!
  if(mult % d_func) then  ! using distance function

    !----------------------------------------------------------------------!
    !   Curvature smoothing was engaged                                    !
    !                                                                      !
    !   Smooth distance function (why n?) and store it in smooth_k array   !
    !----------------------------------------------------------------------!
    if (mult % n_conv_curv > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, mult % dist_func % n, &
                                     smooth_k(-nb:nc), mult % n_conv_curv)

    !----------------------------------------------------------!
    !   Curvature smoothing was not engaged                    !
    !                                                          !
    !   Store distance function's (why n?) in smooth_k array   !
    !----------------------------------------------------------!
    else
      smooth_k(-nb:nc) = mult % dist_func % n
    end if
    call Field_Mod_Interpolate_Cells_To_Nodes(flow,                     &
                                        smooth_k(-nb:nc), var_node_k(1:nn))

    !----------------------------------------------------------------------!
    !   Smoothing of normal is engaged                                     !
    !                                                                      !
    !   Smooth distance function (why n?) and store it in smooth_n array   !
    !----------------------------------------------------------------------!
    if (mult % n_conv_norm > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, mult % dist_func % n, &
                                     smooth_n(-nb:nc), mult % n_conv_norm)

    !-----------------------------------------------------------!
    !   Smoothing of normal is not engaged                      !
    !                                                           !
    !   Store distance function's (why oo?) in smooth_n array   !
    !-----------------------------------------------------------!
    else
      smooth_n(-nb:nc) = mult % dist_func % oo(-nb:nc)
    end if
    call Field_Mod_Interpolate_Cells_To_Nodes(flow,                     &
                                        smooth_n(-nb:nc), var_node_n(1:nn))

  !-----------------------------!
  !                             !
  !   Volume of fluid is used   !
  !                             !
  !-----------------------------!
  else ! using VOF

    !--------------------------------------------------------------------!
    !   Curvature smoothing was engaged                                  !
    !                                                                    !
    !   Smooth volume of fluid function and store it in smooth_k array   !
    !--------------------------------------------------------------------!
    if (mult % n_conv_curv > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, vof % n, &
                                     smooth_k(-nb:nc), mult % n_conv_curv)

    !------------------------------------------------------!
    !   Curvature smoothing was not engaged                !
    !                                                      !
    !   Store volume of fluid function in smooth_k array   !
    !------------------------------------------------------!
    else
      smooth_k(-nb:nc) = vof % n(-nb:nc)
    end if
    call Field_Mod_Interpolate_Cells_To_Nodes(flow,                     &
                                        smooth_k(-nb:nc), var_node_k(1:nn))

    !--------------------------------------------------------------------!
    !   Smoothing of normal is engaged                                   !
    !                                                                    !
    !   Smooth volume of fluid function and store it in smooth_n array   !
    !--------------------------------------------------------------------!
    if (mult % n_conv_norm > 0) then
      call Multiphase_Mod_Vof_Smooth_Scalar(grid, mult, vof % n, &
                                     smooth_n(-nb:nc), mult % n_conv_norm)

    !------------------------------------------------------!
    !   Smoothing of normal is not engaged                 !
    !                                                      !
    !   Store volume of fluid function in smooth_n array   !
    !------------------------------------------------------!
    else
      smooth_n(-nb:nc) = vof % n(-nb:nc)
    end if
    call Field_Mod_Interpolate_Cells_To_Nodes(flow,                     &
                                        smooth_n(-nb:nc), var_node_n(1:nn))

  end if

  !----------------------------------------------------------------------!
  !                                                                      !
  !   No matter what you did above, compute curvature from nodes using   !
  !   the smoothed variable in smooth_k and nodal values in var_node_k   !
  !                                                                      !
  !----------------------------------------------------------------------!
  call Multiphase_Mod_Vof_Curvature_Nodal(grid, mult,  &
                                          smooth_k(-nb:nc), var_node_k(1:nn))

  end subroutine
