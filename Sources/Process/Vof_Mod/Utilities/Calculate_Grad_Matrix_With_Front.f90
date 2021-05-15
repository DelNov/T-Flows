!==============================================================================!
  subroutine Calculate_Grad_Matrix_With_Front(Vof)
!------------------------------------------------------------------------------!
!   Calculates gradient matrix for cases with front tracking (boiling)         !
!                                                                              !
!   (Closely related to Field_Mod_Calculate_Grad_Matrix)                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type) :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  integer                   :: c, c1, c2, s
  real                      :: dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2
  real                      :: jac, g_inv(6)
!==============================================================================!

  ! Take alias
  grid => Vof % pnt_grid
  flow => Vof % pnt_flow

  !--------------------------------------------!
  !   Initialize gradient matrices for cells   !
  !--------------------------------------------!
  do c = 1, grid % n_cells
    flow % grad_c2c(1,c) = 0.0
    flow % grad_c2c(2,c) = 0.0
    flow % grad_c2c(3,c) = 0.0
    flow % grad_c2c(4,c) = 0.0
    flow % grad_c2c(5,c) = 0.0
    flow % grad_c2c(6,c) = 0.0
  end do

  !----------------------------------------------------------------------!
  !   Compute the gradient matrix for all cells browsing through faces   !
  !----------------------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    dx_c1 = grid % dx(s)
    dy_c1 = grid % dy(s)
    dz_c1 = grid % dz(s)
    dx_c2 = grid % dx(s)
    dy_c2 = grid % dy(s)
    dz_c2 = grid % dz(s)

    ! If face is at the front, reduce the extents of the stencil
    if(any(Vof % Front % face_at_elem(1:2,s) .ne. 0)) then
      dx_c1 = grid % xs(s) - grid % xc(c1)
      dy_c1 = grid % ys(s) - grid % yc(c1)
      dz_c1 = grid % zs(s) - grid % zc(c1)
      dx_c2 = grid % xc(c2) - grid % xs(s)
      dy_c2 = grid % yc(c2) - grid % ys(s)
      dz_c2 = grid % zc(c2) - grid % zs(s)
    end if

    flow % grad_c2c(1,c1)=flow % grad_c2c(1,c1) + dx_c1*dx_c1    ! 1,1
    flow % grad_c2c(2,c1)=flow % grad_c2c(2,c1) + dy_c1*dy_c1    ! 2,2
    flow % grad_c2c(3,c1)=flow % grad_c2c(3,c1) + dz_c1*dz_c1    ! 3,3
    flow % grad_c2c(4,c1)=flow % grad_c2c(4,c1) + dx_c1*dy_c1    ! 1,2  &  2,1
    flow % grad_c2c(5,c1)=flow % grad_c2c(5,c1) + dx_c1*dz_c1    ! 1,3  &  3,1
    flow % grad_c2c(6,c1)=flow % grad_c2c(6,c1) + dy_c1*dz_c1    ! 2,3  &  3,2
    if(c2 > 0) then  ! this is enough even for parallel
      flow % grad_c2c(1,c2)=flow % grad_c2c(1,c2) + dx_c2*dx_c2  ! 1,1
      flow % grad_c2c(2,c2)=flow % grad_c2c(2,c2) + dy_c2*dy_c2  ! 2,2
      flow % grad_c2c(3,c2)=flow % grad_c2c(3,c2) + dz_c2*dz_c2  ! 3,3
      flow % grad_c2c(4,c2)=flow % grad_c2c(4,c2) + dx_c2*dy_c2  ! 1,2  &  2,1
      flow % grad_c2c(5,c2)=flow % grad_c2c(5,c2) + dx_c2*dz_c2  ! 1,3  &  3,1
      flow % grad_c2c(6,c2)=flow % grad_c2c(6,c2) + dy_c2*dz_c2  ! 2,3  &  3,2
    end if

  end do

  !--------------------------------!
  !   Find the inverse of matrix   !
  !--------------------------------!
  do c = 1, grid % n_cells
    jac = flow % grad_c2c(1,c) * flow % grad_c2c(2,c) * flow % grad_c2c(3,c)  &
        - flow % grad_c2c(1,c) * flow % grad_c2c(6,c) * flow % grad_c2c(6,c)  &
        - flow % grad_c2c(4,c) * flow % grad_c2c(4,c) * flow % grad_c2c(3,c)  &
        + flow % grad_c2c(4,c) * flow % grad_c2c(5,c) * flow % grad_c2c(6,c)  &
        + flow % grad_c2c(4,c) * flow % grad_c2c(5,c) * flow % grad_c2c(6,c)  &
        - flow % grad_c2c(5,c) * flow % grad_c2c(5,c) * flow % grad_c2c(2,c)

    g_inv(1) = +(  flow % grad_c2c(2,c) * flow % grad_c2c(3,c)  &
                 - flow % grad_c2c(6,c) * flow % grad_c2c(6,c) ) / (jac+TINY)
    g_inv(2) = +(  flow % grad_c2c(1,c) * flow % grad_c2c(3,c)  &
                 - flow % grad_c2c(5,c) * flow % grad_c2c(5,c) ) / (jac+TINY)
    g_inv(3) = +(  flow % grad_c2c(1,c) * flow % grad_c2c(2,c)  &
                 - flow % grad_c2c(4,c) * flow % grad_c2c(4,c) ) / (jac+TINY)
    g_inv(4) = -(  flow % grad_c2c(4,c) * flow % grad_c2c(3,c)  &
                 - flow % grad_c2c(5,c) * flow % grad_c2c(6,c) ) / (jac+TINY)
    g_inv(5) = +(  flow % grad_c2c(4,c) * flow % grad_c2c(6,c)  &
                 - flow % grad_c2c(5,c) * flow % grad_c2c(2,c) ) / (jac+TINY)
    g_inv(6) = -(  flow % grad_c2c(1,c) * flow % grad_c2c(6,c)  &
                 - flow % grad_c2c(4,c) * flow % grad_c2c(5,c) ) / (jac+TINY)

    flow % grad_c2c(1,c) = g_inv(1)
    flow % grad_c2c(2,c) = g_inv(2)
    flow % grad_c2c(3,c) = g_inv(3)
    flow % grad_c2c(4,c) = g_inv(4)
    flow % grad_c2c(5,c) = g_inv(5)
    flow % grad_c2c(6,c) = g_inv(6)
  end do

  end subroutine
