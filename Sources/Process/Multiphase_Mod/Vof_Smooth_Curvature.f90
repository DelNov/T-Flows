!==============================================================================!
  subroutine Multiphase_Mod_Vof_Smooth_Curvature(mult)
!------------------------------------------------------------------------------!
!   Smoothes curvature in two steps: first a smoothing curvature around the    !
!   Interface and second in the direction of the normal. This technique can    !
!   be found at https://spiral.imperial.ac.uk/handle/10044/1/28101
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: k_star       => r_cell_14,    &
                      gradk_x      => r_cell_15,    &
                      gradk_y      => r_cell_16,    &
                      gradk_z      => r_cell_17,    &
                      sum_k_weight => r_cell_18,    &
                      sum_weight   => r_cell_19
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: col   ! colour; could be vof or smooth
  integer                  :: s, c, c1, c2, nb, nc
  real                     :: fs, w_v1, w_v2, w_m1, w_m2
  real                     :: weight_s, weight_n
  real                     :: curvf
!==============================================================================!

  ! Take aliases
  grid => mult % pnt_grid

  col => mult % smooth
  ! col => mult % vof

  nb = grid % n_bnd_cells
  nc = grid % n_cells

  sum_k_weight(-nb:nc) = 0.0
  sum_weight  (-nb:nc) = 0.0
  weight_s             = 8.0
  weight_n             = 8.0

  !-------------------------!
  !   Smoothing curvature   !
  !-------------------------!

  ! What is the curvature at boundaries??? For now zero gradient
  ! Preliminary results using wetting (with symmetry in some boundaries) show
  ! it is better not to take into aacount the boundaries

  !  gradient of curvature, only interior
  gradk_x(-nb:nc) = 0.0
  gradk_y(-nb:nc) = 0.0
  gradk_z(-nb:nc) = 0.0

  do s = grid % n_bnd_faces + 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    curvf = fs * mult % curv(c1) + (1.0 - fs) * mult % curv(c2)

    gradk_x(c1) = gradk_x(c1) + curvf * grid % sx(s) / grid % vol(c1)
    gradk_y(c1) = gradk_y(c1) + curvf * grid % sy(s) / grid % vol(c1)
    gradk_z(c1) = gradk_z(c1) + curvf * grid % sz(s) / grid % vol(c1)

    gradk_x(c2) = gradk_x(c2) - curvf * grid % sx(s) / grid % vol(c2)
    gradk_y(c2) = gradk_y(c2) - curvf * grid % sy(s) / grid % vol(c2)
    gradk_z(c2) = gradk_z(c2) - curvf * grid % sz(s) / grid % vol(c2)

  end do

  call Grid_Mod_Exchange_Cells_Real(grid, gradk_x(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, gradk_y(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, gradk_z(-nb:nc))

  ! Interior faces
  do s = grid % n_bnd_faces + 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    w_v1 = (1.0 - 2.0 * abs(0.5 - col % n(c1))) ** weight_s
    w_v2 = (1.0 - 2.0 * abs(0.5 - col % n(c2))) ** weight_s
    sum_k_weight(c1) = sum_k_weight(c1) + mult % curv(c2) * w_v2
    sum_weight(c1) = sum_weight(c1) + w_v2

    sum_k_weight(c2) = sum_k_weight(c2) + mult % curv(c1) * w_v1
    sum_weight(c2) = sum_weight(c2) + w_v1
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, sum_k_weight(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, sum_weight  (-nb:nc))

  do c = 1, grid % n_cells
    w_v1 = (1.0 - 2.0 * abs(0.5 - col % n(c))) ** weight_s
    k_star(c) = (w_v1 * mult % curv(c) + sum_k_weight(c))    &
              / (w_v1 + sum_weight(c) + FEMTO)
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, k_star(-nb:nc))

  !-------------------------------------------------------------------------!
  !   Smoothing curvature in the direction of the normal to the interface   !
  !-------------------------------------------------------------------------!

  sum_k_weight(-nb:nc) = 0.0
  sum_weight  (-nb:nc) = 0.0

  !  gradient of curvature, only interior
  gradk_x(-nb:nc) = 0.0
  gradk_y(-nb:nc) = 0.0
  gradk_z(-nb:nc) = 0.0

  do s = grid % n_bnd_faces + 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    curvf = fs * k_star(c1) + (1.0 - fs) * k_star(c2)

    gradk_x(c1) = gradk_x(c1) + curvf * grid % sx(s) / grid % vol(c1)
    gradk_y(c1) = gradk_y(c1) + curvf * grid % sy(s) / grid % vol(c1)
    gradk_z(c1) = gradk_z(c1) + curvf * grid % sz(s) / grid % vol(c1)

    gradk_x(c2) = gradk_x(c2) - curvf * grid % sx(s) / grid % vol(c2)
    gradk_y(c2) = gradk_y(c2) - curvf * grid % sy(s) / grid % vol(c2)
    gradk_z(c2) = gradk_z(c2) - curvf * grid % sz(s) / grid % vol(c2)

  end do

  call Grid_Mod_Exchange_Cells_Real(grid, gradk_x(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, gradk_y(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, gradk_z(-nb:nc))

  ! Interior faces
  do s = grid % n_bnd_faces + 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    w_v1 = (1.0 - 2.0 * abs(0.5 - col % n(c1))) ** weight_n
    w_v2 = (1.0 - 2.0 * abs(0.5 - col % n(c2))) ** weight_n

    w_m1 = abs(dot_product((/mult % nx(c1), mult % ny(c1), mult % nz(c1)/),  &
                           (/grid % dx(s),  grid % dy(s),  grid % dz(s)/)    &
                           / grid % d(s))) ** weight_n

    w_m2 = abs(dot_product((/mult % nx(c2), mult % ny(c2), mult % nz(c2)/),  &
                           (/grid % dx(s),  grid % dy(s),  grid % dz(s)/)    &
                           / (-grid % d(s)))) ** weight_n

    sum_k_weight(c1) = sum_k_weight(c1) + k_star(c2) * w_v2 * w_m2
    sum_weight(c1) = sum_weight(c1) + w_v2 * w_m2

    sum_k_weight(c2) = sum_k_weight(c2) + k_star(c1) * w_v1 * w_m1
    sum_weight(c2) = sum_weight(c2) + w_v1 * w_m1
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, sum_k_weight(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, sum_weight  (-nb:nc))

  do c = 1, grid % n_cells
    w_v1 = (1.0 - 2.0 * abs(0.5 - col % n(c))) ** weight_n
    mult % curv(c) = (w_v1 * k_star(c) + sum_k_weight(c))    &
                   / (w_v1 + sum_weight(c) + FEMTO)
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, mult % curv)

  end subroutine
