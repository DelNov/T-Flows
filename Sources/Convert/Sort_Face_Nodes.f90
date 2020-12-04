!==============================================================================!
  subroutine Sort_Face_Nodes(grid, s)
!------------------------------------------------------------------------------!
!   Sort nodes of a given face in rotational fashion                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: s
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, k, n, nn
  real    :: a(3,MAX_FACES_N_NODES), b(MAX_FACES_N_NODES)
  real    :: a_p(3,3), b_p(3), det
  real    :: normal_p(3), center_p(3), x_p(3), y_p(3)
  real    :: rp_3d(3,MAX_FACES_N_NODES)
  real    :: rp_2d(2,MAX_FACES_N_NODES)
  real    :: sorting(MAX_FACES_N_NODES)  ! sorting criterion
  integer :: order(MAX_FACES_N_NODES)    ! carry-on array with indices
!==============================================================================!

  ! Take alias
  nn = grid % faces_n_nodes(s)

  !----------------------------------!
  !   Find the planes center point   !
  !----------------------------------!
  center_p(:) = 0
  do i = 1, nn
    n = grid % faces_n(i, s)
    center_p(1) = center_p(1) + grid % xn(n)
    center_p(2) = center_p(2) + grid % yn(n)
    center_p(3) = center_p(3) + grid % zn(n)
  end do
  center_p(1:3) = center_p(1:3) / real(nn)

  !---------------------------!
  !   Find relative vectors   !
  !---------------------------!
  do i = 1, nn  ! use "i", not "i_nod" here
    n = grid % faces_n(i, s)
    rp_3d(1, i) = grid % xn(n) - center_p(1)
    rp_3d(2, i) = grid % yn(n) - center_p(2)
    rp_3d(3, i) = grid % zn(n) - center_p(3)
  end do

  !--------------------------------!
  !   Find a normal of the plane   !
  !--------------------------------!

  ! Fill up the original matrices
  do i = 1, nn  ! use "i", not "i_nod" here
    n = grid % faces_n(i, s)
    a(1,i) = rp_3d(1,i) + MICRO * rand() ! add some noise on coords ...
    a(2,i) = rp_3d(2,i) + MICRO * rand() ! ... to avoid over-dermined ...
    a(3,i) = rp_3d(3,i) + MICRO * rand() ! ... system of equations
    b(i)   = 1.0
  end do

  ! Form the system for least squares method
  a_p(:,:) = 0.0
  do j = 1, 3
    do k = 1, 3
      do i = 1, nn
        a_p(j,k) = a_p(j,k) + a(j,i) * a(k,i)
      end do
    end do
  end do

  b_p(:) = 0.0
  do j = 1, 3
    do i = 1, nn
      b_p(j) = b_p(j) + a(j,i) * b(i)
    end do
  end do

  ! Solve the system
  call Math_Mod_Invert_3x3(a_p, det)
  normal_p(:) = 0.0
  do j = 1, 3
    do i = 1, 3
      normal_p(j) = normal_p(j) + a_p(j,i) * b_p(i)
    end do
  end do

  normal_p(1:3) = normal_p(1:3) / norm2(normal_p(1:3))

  ! Define x in the plane
  x_p(1:3) = Math_Mod_Cross_Product(normal_p(1:3), rp_3d(1:3, 1))
  x_p(1:3) = -x_p(1:3) / norm2(x_p(1:3))

  ! Define y in the plane
  y_p(1:3) = Math_Mod_Cross_Product(normal_p(1:3), x_p(1:3))
  y_p(1:3) = y_p(1:3) / norm2(y_p(1:3))

  ! Project relative points to the plane's coordinate system
  do i = 1, nn  ! use "i", not "i_nod" here
    n = grid % faces_n(i, s)
    rp_2d(1, i) = dot_product(x_p(1:3), rp_3d(1:3, i))
    rp_2d(2, i) = dot_product(y_p(1:3), rp_3d(1:3, i))
  end do

  do i = 1, nn  ! use "i", not "i_nod" here
    sorting(i) = atan2(rp_2d(1,i), rp_2d(2,i)) * 57.2957795131
    order(i) = i  ! like old number
  end do

  call Sort_Mod_Real_Carry_2_Int(sorting(1:nn),           &
                                 grid % faces_n(1:nn,s),  &
                                 order(1:nn))

  ! Rotate for better visualisation with Paraview.  Paraview constructs
  ! a polygon from triangles which all originate in the first node.  If
  ! the first node is defined exactly at the sharpest edge, which is also
  ! the first one to be inserted in polynomial face in the calling func.
  ! All in all: try to keep the first node first for better visuals.
  do while(order(1) .ne. 1)
    order(1:nn)             = cshift(order(1:nn), 1)
    grid % faces_n(1:nn, s) = cshift(grid % faces_n(1:nn, s), 1)
  end do

  end subroutine
