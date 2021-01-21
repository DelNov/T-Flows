!==============================================================================!
  subroutine Interpolate_From_Nodes(grid, mult, var_node, s, i_probe, i_s)
!------------------------------------------------------------------------------!
!   This function interpolates nodal values to probe on face s                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type),       target :: grid
  type(Multiphase_Type), target :: mult
  real                          :: var_node(1 : grid % n_nodes)
  integer                       :: s
  integer                       :: i_probe, i_s
!-----------------------------------[Locals]-----------------------------------!
  integer             :: n, i_nod   ! counters
  real,   allocatable :: weights(:)
  real,   pointer     :: s_coor(:)
  real                :: rx, ry, rz, lambda_x, lambda_y, lambda_z
  real                :: ixx, iyy, izz, ixz, iyz, ixy, d
  real                :: a11, a12, a13, a21, a22, a23, a31, a32, a33
  real                :: corr_x, corr_y, corr_z
  real                :: xf, yf, zf
  real                :: epsloc, sum1, sum2
!==============================================================================!

  ! Take alias
  s_coor => probes(i_probe) % s_coor(i_s,:)

  epsloc = epsilon(epsloc)

  !---------------------------------------!
  !   Browse through all nodes to form    !
  !              node weights             !
  !---------------------------------------!
  allocate(weights(grid % faces_n_nodes(s)))
  weights = 0.0

  rx = 0.0; ry = 0.0; rz = 0.0
  ixx = 0.0; iyy = 0.0; izz = 0.0; ixz = 0.0; iyz = 0.0; ixy = 0.0

  ! Loop on nodes
  do i_nod = 1, grid % faces_n_nodes(s)
    n = grid % faces_n(i_nod, s)
    rx = rx + (grid % xn(n) - s_coor(1))
    ry = ry + (grid % yn(n) - s_coor(2))
    rz = rz + (grid % zn(n) - s_coor(3))
    ixx = ixx + (grid % xn(n) - s_coor(1)) ** 2
    iyy = iyy + (grid % yn(n) - s_coor(2)) ** 2
    izz = izz + (grid % zn(n) - s_coor(3)) ** 2
    ixy = ixy + (grid % xn(n) - s_coor(1)) * (grid % yn(n) - s_coor(2))
    ixz = ixz + (grid % xn(n) - s_coor(1)) * (grid % zn(n) - s_coor(3))
    iyz = iyz + (grid % yn(n) - s_coor(2)) * (grid % zn(n) - s_coor(3))
  end do

  a11 = iyy * izz - iyz ** 2
  a12 = ixz * iyz - ixy * izz
  a13 = ixy * iyz - ixz * iyy
  a21 = ixz * iyz - ixy * izz
  a22 = ixx * izz - ixz ** 2
  a23 = ixy * ixz - ixx * iyz
  a31 = ixy * iyz - ixz * iyy
  a32 = ixz * ixy - ixx * iyz
  a33 = ixx * iyy - ixy ** 2

  d = ixx * iyy * izz - ixx * iyz **2 - iyy * ixz ** 2 - izz * ixy ** 2   &
    + 2.0 * ixy * ixz * iyz

  lambda_x = (rx * a11 + ry * a12 + rz * a13) / (d + epsloc)
  lambda_y = (rx * a21 + ry * a22 + rz * a23) / (d + epsloc)
  lambda_z = (rx * a31 + ry * a32 + rz * a33) / (d + epsloc)

  do i_nod = 1, grid % faces_n_nodes(s)
    n = grid % faces_n(i_nod, s)
    weights(i_nod) = 1.0                                     &
                   + lambda_x * (grid % xn(n) - s_coor(1))   &
                   + lambda_y * (grid % yn(n) - s_coor(2))   &
                   + lambda_z * (grid % zn(n) - s_coor(3))
  end do

  ! Interpolate
  sum1 = 0.0; sum2 = 0.0

  do i_nod = 1, grid % faces_n_nodes(s)
    n = grid % faces_n(i_nod, s)
    sum1 = sum1 + weights(i_nod) * var_node(n)
    sum2 = sum2 + weights(i_nod)
  end do
  probes(i_probe) % s_vof(i_s) = sum1 / sum2

  end subroutine
