!==============================================================================!
  subroutine Smooth_Curvature(Vof)
!------------------------------------------------------------------------------!
!   Smoothes curvature in two steps: first a smoothing curvature around the    !
!   Interface and second in the direction of the normal. This technique can    !
!   be found at https://spiral.imperial.ac.uk/handle/10044/1/28101
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer  :: Grid
  type(Var_Type),  pointer  :: col   ! colour; could be vof or smooth
  integer                   :: s, c, c1, c2, nb, nc
  real                      :: fs, w_v1, w_v2, w_m1, w_m2
  real                      :: weight_s, weight_n
  real                      :: curvf
  real, contiguous, pointer :: gradk_x(:), gradk_y(:), gradk_z(:)
  real, contiguous, pointer :: k_star(:), sum_k_weight(:), sum_weight(:)
!==============================================================================!

  call Work % Connect_Real_Cell(gradk_x, gradk_y, gradk_z,  &
                                k_star, sum_k_weight, sum_weight)

  ! Take aliases
  Grid => Vof % pnt_grid

  col => Vof % smooth
  ! col => Vof % vof

  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

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

  ! Gradient of curvature, only interior
  gradk_x(-nb:nc) = 0.0
  gradk_y(-nb:nc) = 0.0
  gradk_z(-nb:nc) = 0.0

  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    fs = Grid % f(s)
    curvf = fs * Vof % curv(c1) + (1.0 - fs) * Vof % curv(c2)
    gradk_x(c1) = gradk_x(c1) + curvf * Grid % sx(s) / Grid % vol(c1)
    gradk_y(c1) = gradk_y(c1) + curvf * Grid % sy(s) / Grid % vol(c1)
    gradk_z(c1) = gradk_z(c1) + curvf * Grid % sz(s) / Grid % vol(c1)
    gradk_x(c2) = gradk_x(c2) - curvf * Grid % sx(s) / Grid % vol(c2)
    gradk_y(c2) = gradk_y(c2) - curvf * Grid % sy(s) / Grid % vol(c2)
    gradk_z(c2) = gradk_z(c2) - curvf * Grid % sz(s) / Grid % vol(c2)
  end do  ! faces

  call Grid % Exchange_Cells_Real(gradk_x(-nb:nc))
  call Grid % Exchange_Cells_Real(gradk_y(-nb:nc))
  call Grid % Exchange_Cells_Real(gradk_z(-nb:nc))

  ! Interior faces
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    w_v1 = (1.0 - 2.0 * abs(0.5 - col % n(c1))) ** weight_s
    w_v2 = (1.0 - 2.0 * abs(0.5 - col % n(c2))) ** weight_s
    sum_k_weight(c1) = sum_k_weight(c1) + Vof % curv(c2) * w_v2
    sum_weight(c1)   = sum_weight(c1)   + w_v2

    sum_k_weight(c2) = sum_k_weight(c2) + Vof % curv(c1) * w_v1
    sum_weight(c2)   = sum_weight(c2)   + w_v1
  end do  ! faces
  call Grid % Exchange_Cells_Real(sum_k_weight(-nb:nc))
  call Grid % Exchange_Cells_Real(sum_weight  (-nb:nc))

  do c = Cells_In_Domain()
    w_v1 = (1.0 - 2.0 * abs(0.5 - col % n(c))) ** weight_s
    k_star(c) = (w_v1 * Vof % curv(c) + sum_k_weight(c))    &
              / (w_v1 + sum_weight(c) + FEMTO)
  end do
  call Grid % Exchange_Cells_Real(k_star(-nb:nc))

  !-------------------------------------------------------------------------!
  !   Smoothing curvature in the direction of the normal to the interface   !
  !-------------------------------------------------------------------------!

  sum_k_weight(-nb:nc) = 0.0
  sum_weight  (-nb:nc) = 0.0

  ! Gradient of curvature, only interior
  gradk_x(-nb:nc) = 0.0
  gradk_y(-nb:nc) = 0.0
  gradk_z(-nb:nc) = 0.0

  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    fs = Grid % f(s)
    curvf = fs * k_star(c1) + (1.0 - fs) * k_star(c2)
    gradk_x(c1) = gradk_x(c1) + curvf * Grid % sx(s) / Grid % vol(c1)
    gradk_y(c1) = gradk_y(c1) + curvf * Grid % sy(s) / Grid % vol(c1)
    gradk_z(c1) = gradk_z(c1) + curvf * Grid % sz(s) / Grid % vol(c1)
    gradk_x(c2) = gradk_x(c2) - curvf * Grid % sx(s) / Grid % vol(c2)
    gradk_y(c2) = gradk_y(c2) - curvf * Grid % sy(s) / Grid % vol(c2)
    gradk_z(c2) = gradk_z(c2) - curvf * Grid % sz(s) / Grid % vol(c2)
  end do
  call Grid % Exchange_Cells_Real(gradk_x(-nb:nc))
  call Grid % Exchange_Cells_Real(gradk_y(-nb:nc))
  call Grid % Exchange_Cells_Real(gradk_z(-nb:nc))

  ! Interior faces
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    w_v1 = (1.0 - 2.0 * abs(0.5 - col % n(c1))) ** weight_n
    w_v2 = (1.0 - 2.0 * abs(0.5 - col % n(c2))) ** weight_n

    w_m1 = abs(dot_product((/Vof % nx(c1), Vof % ny(c1), Vof % nz(c1)/),  &
                           (/Grid % dx(s), Grid % dy(s), Grid % dz(s)/)   &
                           / Grid % d(s))) ** weight_n

    w_m2 = abs(dot_product((/Vof % nx(c2), Vof % ny(c2), Vof % nz(c2)/),  &
                           (/Grid % dx(s), Grid % dy(s), Grid % dz(s)/)   &
                           / (-Grid % d(s)))) ** weight_n

    sum_k_weight(c1) = sum_k_weight(c1) + k_star(c2) * w_v2 * w_m2
    sum_weight(c1)   = sum_weight(c1)   + w_v2 * w_m2

    sum_k_weight(c2) = sum_k_weight(c2) + k_star(c1) * w_v1 * w_m1
    sum_weight(c2)   = sum_weight(c2)   + w_v1 * w_m1
  end do
  call Grid % Exchange_Cells_Real(sum_k_weight(-nb:nc))
  call Grid % Exchange_Cells_Real(sum_weight  (-nb:nc))

  do c = Cells_In_Domain()
    w_v1 = (1.0 - 2.0 * abs(0.5 - col % n(c))) ** weight_n
    Vof % curv(c) = (w_v1 * k_star(c) + sum_k_weight(c))    &
                   / (w_v1 + sum_weight(c) + FEMTO)
  end do

  call Grid % Exchange_Cells_Real(Vof % curv)

  call Work % Disconnect_Real_Cell(gradk_x, gradk_y, gradk_z,  &
                                   k_star, sum_k_weight, sum_weight)

  end subroutine
