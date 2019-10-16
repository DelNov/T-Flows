!==============================================================================!
  subroutine Compute_Gradient_Matrix(grid)
!------------------------------------------------------------------------------!
!   Calculates gradient matrix.                                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s
  real    :: dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2
  real    :: jac, g_inv(6)
!==============================================================================!

  do c = 1, grid % n_cells
    g(1,c) = 0.0
    g(2,c) = 0.0
    g(3,c) = 0.0
    g(4,c) = 0.0
    g(5,c) = 0.0
    g(6,c) = 0.0
  end do

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s) 

    dx_c1 = grid % dx(s)
    dy_c1 = grid % dy(s)
    dz_c1 = grid % dz(s)
    dx_c2 = grid % dx(s)
    dy_c2 = grid % dy(s)
    dz_c2 = grid % dz(s)

    g(1,c1)=g(1,c1) + dx_c1*dx_c1    ! 1,1
    g(2,c1)=g(2,c1) + dy_c1*dy_c1    ! 2,2
    g(3,c1)=g(3,c1) + dz_c1*dz_c1    ! 3,3
    g(4,c1)=g(4,c1) + dx_c1*dy_c1    ! 1,2  &  2,1
    g(5,c1)=g(5,c1) + dx_c1*dz_c1    ! 1,3  &  3,1
    g(6,c1)=g(6,c1) + dy_c1*dz_c1    ! 2,3  &  3,2
    if(c2 > 0) then  ! this is enough even for parallel
      g(1,c2)=g(1,c2) + dx_c2*dx_c2  ! 1,1
      g(2,c2)=g(2,c2) + dy_c2*dy_c2  ! 2,2
      g(3,c2)=g(3,c2) + dz_c2*dz_c2  ! 3,3
      g(4,c2)=g(4,c2) + dx_c2*dy_c2  ! 1,2  &  2,1
      g(5,c2)=g(5,c2) + dx_c2*dz_c2  ! 1,3  &  3,1
      g(6,c2)=g(6,c2) + dy_c2*dz_c2  ! 2,3  &  3,2
    end if

  end do

  !----------------------------------!
  !   Find the inverse of matrix g   !
  !----------------------------------!
  do c = 1, grid % n_cells
    jac  =         g(1,c) * g(2,c) * g(3,c)  &
           -       g(1,c) * g(6,c) * g(6,c)  &
           -       g(4,c) * g(4,c) * g(3,c)  &
           + 2.0 * g(4,c) * g(5,c) * g(6,c)  &
           -       g(5,c) * g(5,c) * g(2,c)

    g_inv(1) = +( g(2,c)*g(3,c) - g(6,c)*g(6,c) ) / (jac+TINY)
    g_inv(2) = +( g(1,c)*g(3,c) - g(5,c)*g(5,c) ) / (jac+TINY)
    g_inv(3) = +( g(1,c)*g(2,c) - g(4,c)*g(4,c) ) / (jac+TINY)
    g_inv(4) = -( g(4,c)*g(3,c) - g(5,c)*g(6,c) ) / (jac+TINY)
    g_inv(5) = +( g(4,c)*g(6,c) - g(5,c)*g(2,c) ) / (jac+TINY)
    g_inv(6) = -( g(1,c)*g(6,c) - g(4,c)*g(5,c) ) / (jac+TINY)

    g(1,c) = g_inv(1) 
    g(2,c) = g_inv(2)
    g(3,c) = g_inv(3)
    g(4,c) = g_inv(4)
    g(5,c) = g_inv(5)
    g(6,c) = g_inv(6)
  end do 

  end subroutine
