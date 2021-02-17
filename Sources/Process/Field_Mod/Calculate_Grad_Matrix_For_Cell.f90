!==============================================================================!
  subroutine Field_Mod_Calculate_Grad_Matrix_For_Cell(flow, c)
!------------------------------------------------------------------------------!
!   Calculates gradient matrix.                                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  integer          :: c
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, i_fac
  real                     :: dx_c, dy_c, dz_c
  real                     :: jac, g_inv(6)
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  !-------------------------------------------!
  !   Initialize gradient matrices for cell   !
  !-------------------------------------------!
  flow % grad_c2c(1:6,c) = 0.0

  !--------------------------------------!
  !   Form gradient matrices for cells   !
  !--------------------------------------!
  do i_fac = 1, grid % cells_n_faces(c)
    s = grid % cells_f(i_fac, c)

    dx_c = grid % dx(s)
    dy_c = grid % dy(s)
    dz_c = grid % dz(s)

    flow % grad_c2c(1,c)=flow % grad_c2c(1,c) + dx_c**2    ! 1,1
    flow % grad_c2c(2,c)=flow % grad_c2c(2,c) + dy_c**2    ! 2,2
    flow % grad_c2c(3,c)=flow % grad_c2c(3,c) + dz_c**2    ! 3,3
    flow % grad_c2c(4,c)=flow % grad_c2c(4,c) + dx_c*dy_c  ! 1,2  &  2,1
    flow % grad_c2c(5,c)=flow % grad_c2c(5,c) + dx_c*dz_c  ! 1,3  &  3,1
    flow % grad_c2c(6,c)=flow % grad_c2c(6,c) + dy_c*dz_c  ! 2,3  &  3,2

  end do  ! i_fac

  !--------------------------------!
  !   Find the inverse of matrix   !
  !--------------------------------!
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

  end subroutine
