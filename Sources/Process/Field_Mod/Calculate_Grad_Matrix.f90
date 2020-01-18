!==============================================================================!
  subroutine Field_Mod_Calculate_Grad_Matrix(flow)
!------------------------------------------------------------------------------!
!   Calculates gradient matrix.                                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c, c1, c2, s
  real                     :: dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2
  real                     :: jac, g_inv(6)
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  do c = 1, grid % n_cells
    flow % grad(1,c) = 0.0
    flow % grad(2,c) = 0.0
    flow % grad(3,c) = 0.0
    flow % grad(4,c) = 0.0
    flow % grad(5,c) = 0.0
    flow % grad(6,c) = 0.0
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

    flow % grad(1,c1)=flow % grad(1,c1) + dx_c1*dx_c1    ! 1,1
    flow % grad(2,c1)=flow % grad(2,c1) + dy_c1*dy_c1    ! 2,2
    flow % grad(3,c1)=flow % grad(3,c1) + dz_c1*dz_c1    ! 3,3
    flow % grad(4,c1)=flow % grad(4,c1) + dx_c1*dy_c1    ! 1,2  &  2,1
    flow % grad(5,c1)=flow % grad(5,c1) + dx_c1*dz_c1    ! 1,3  &  3,1
    flow % grad(6,c1)=flow % grad(6,c1) + dy_c1*dz_c1    ! 2,3  &  3,2
    if(c2 > 0) then  ! this is enough even for parallel
      flow % grad(1,c2)=flow % grad(1,c2) + dx_c2*dx_c2  ! 1,1
      flow % grad(2,c2)=flow % grad(2,c2) + dy_c2*dy_c2  ! 2,2
      flow % grad(3,c2)=flow % grad(3,c2) + dz_c2*dz_c2  ! 3,3
      flow % grad(4,c2)=flow % grad(4,c2) + dx_c2*dy_c2  ! 1,2  &  2,1
      flow % grad(5,c2)=flow % grad(5,c2) + dx_c2*dz_c2  ! 1,3  &  3,1
      flow % grad(6,c2)=flow % grad(6,c2) + dy_c2*dz_c2  ! 2,3  &  3,2
    end if

  end do

  !--------------------------------!
  !   Find the inverse of matrix   !
  !--------------------------------!
  do c = 1, grid % n_cells
    jac  =       flow % grad(1,c) * flow % grad(2,c) * flow % grad(3,c)  &
         -       flow % grad(1,c) * flow % grad(6,c) * flow % grad(6,c)  &
         -       flow % grad(4,c) * flow % grad(4,c) * flow % grad(3,c)  &
         + 2.0 * flow % grad(4,c) * flow % grad(5,c) * flow % grad(6,c)  &
         -       flow % grad(5,c) * flow % grad(5,c) * flow % grad(2,c)

    g_inv(1) = +(  flow % grad(2,c) * flow % grad(3,c)  &
                 - flow % grad(6,c) * flow % grad(6,c) ) / (jac+TINY)
    g_inv(2) = +(  flow % grad(1,c) * flow % grad(3,c)  &
                 - flow % grad(5,c) * flow % grad(5,c) ) / (jac+TINY)
    g_inv(3) = +(  flow % grad(1,c) * flow % grad(2,c)  &
                 - flow % grad(4,c) * flow % grad(4,c) ) / (jac+TINY)
    g_inv(4) = -(  flow % grad(4,c) * flow % grad(3,c)  &
                 - flow % grad(5,c) * flow % grad(6,c) ) / (jac+TINY)
    g_inv(5) = +(  flow % grad(4,c) * flow % grad(6,c)  &
                 - flow % grad(5,c) * flow % grad(2,c) ) / (jac+TINY)
    g_inv(6) = -(  flow % grad(1,c) * flow % grad(6,c)  &
                 - flow % grad(4,c) * flow % grad(5,c) ) / (jac+TINY)

    flow % grad(1,c) = g_inv(1)
    flow % grad(2,c) = g_inv(2)
    flow % grad(3,c) = g_inv(3)
    flow % grad(4,c) = g_inv(4)
    flow % grad(5,c) = g_inv(5)
    flow % grad(6,c) = g_inv(6)
  end do

  end subroutine
