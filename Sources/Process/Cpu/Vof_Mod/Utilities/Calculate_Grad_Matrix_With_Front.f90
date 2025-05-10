!==============================================================================!
  subroutine Calculate_Grad_Matrix_With_Front(Vof)
!------------------------------------------------------------------------------!
!   Calculates gradient matrix for cases with front tracking (boiling)         !
!                                                                              !
!   (Closely related to Field_Mod/Gradients/Calculate_Grad_Matrix)             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Front_Type), pointer :: Front
  integer                   :: c, c1, c2, s
  real                      :: dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2
  real                      :: jac, g_inv(6)
  real, contiguous, pointer :: g1(:), g2(:), g3(:), g4(:), g5(:), g6(:)
!==============================================================================!

  if(DEBUG) then
    call Work % Connect_Real_Cell(g1, g2, g3, g4, g5, g6)
  end if

  ! Take alias
  Grid  => Vof % pnt_grid
  Flow  => Vof % pnt_flow
  Front => Vof % Front

  !--------------------------------------------!
  !   Initialize gradient matrices for cells   !
  !--------------------------------------------!
  do c = Cells_In_Domain_And_Buffers()
    Flow % grad_c2c(1,c) = 0.0
    Flow % grad_c2c(2,c) = 0.0
    Flow % grad_c2c(3,c) = 0.0
    Flow % grad_c2c(4,c) = 0.0
    Flow % grad_c2c(5,c) = 0.0
    Flow % grad_c2c(6,c) = 0.0
  end do

  !----------------------------------------------------------------------!
  !   Compute the gradient matrix for all cells browsing through faces   !
  !----------------------------------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    dx_c1 = Grid % dx(s)
    dy_c1 = Grid % dy(s)
    dz_c1 = Grid % dz(s)
    dx_c2 = Grid % dx(s)
    dy_c2 = Grid % dy(s)
    dz_c2 = Grid % dz(s)

    ! If face is at the front, reduce the extents of the stencil
    if(Front % intersects_face(s)) then
      dx_c1 = Front % xs(s) - Grid % xc(c1)
      dy_c1 = Front % ys(s) - Grid % yc(c1)
      dz_c1 = Front % zs(s) - Grid % zc(c1)
      dx_c2 = Grid % xc(c2) - Front % xs(s)
      dy_c2 = Grid % yc(c2) - Front % ys(s)
      dz_c2 = Grid % zc(c2) - Front % zs(s)
    end if

    Flow % grad_c2c(1,c1)=Flow % grad_c2c(1,c1) + dx_c1*dx_c1    ! 1,1
    Flow % grad_c2c(2,c1)=Flow % grad_c2c(2,c1) + dy_c1*dy_c1    ! 2,2
    Flow % grad_c2c(3,c1)=Flow % grad_c2c(3,c1) + dz_c1*dz_c1    ! 3,3
    Flow % grad_c2c(4,c1)=Flow % grad_c2c(4,c1) + dx_c1*dy_c1    ! 1,2  &  2,1
    Flow % grad_c2c(5,c1)=Flow % grad_c2c(5,c1) + dx_c1*dz_c1    ! 1,3  &  3,1
    Flow % grad_c2c(6,c1)=Flow % grad_c2c(6,c1) + dy_c1*dz_c1    ! 2,3  &  3,2
    if(c2 > 0) then  ! this is enough even for parallel
      Flow % grad_c2c(1,c2)=Flow % grad_c2c(1,c2) + dx_c2*dx_c2  ! 1,1
      Flow % grad_c2c(2,c2)=Flow % grad_c2c(2,c2) + dy_c2*dy_c2  ! 2,2
      Flow % grad_c2c(3,c2)=Flow % grad_c2c(3,c2) + dz_c2*dz_c2  ! 3,3
      Flow % grad_c2c(4,c2)=Flow % grad_c2c(4,c2) + dx_c2*dy_c2  ! 1,2  &  2,1
      Flow % grad_c2c(5,c2)=Flow % grad_c2c(5,c2) + dx_c2*dz_c2  ! 1,3  &  3,1
      Flow % grad_c2c(6,c2)=Flow % grad_c2c(6,c2) + dy_c2*dz_c2  ! 2,3  &  3,2
    end if

  end do

  !--------------------------------!
  !   Find the inverse of matrix   !
  !--------------------------------!
  do c = Cells_In_Domain()
    jac = Flow % grad_c2c(1,c) * Flow % grad_c2c(2,c) * Flow % grad_c2c(3,c)  &
        - Flow % grad_c2c(1,c) * Flow % grad_c2c(6,c) * Flow % grad_c2c(6,c)  &
        - Flow % grad_c2c(4,c) * Flow % grad_c2c(4,c) * Flow % grad_c2c(3,c)  &
        + Flow % grad_c2c(4,c) * Flow % grad_c2c(5,c) * Flow % grad_c2c(6,c)  &
        + Flow % grad_c2c(4,c) * Flow % grad_c2c(5,c) * Flow % grad_c2c(6,c)  &
        - Flow % grad_c2c(5,c) * Flow % grad_c2c(5,c) * Flow % grad_c2c(2,c)
    Assert(jac > 0)

    g_inv(1) = +(  Flow % grad_c2c(2,c) * Flow % grad_c2c(3,c)  &
                 - Flow % grad_c2c(6,c) * Flow % grad_c2c(6,c) ) / jac
    g_inv(2) = +(  Flow % grad_c2c(1,c) * Flow % grad_c2c(3,c)  &
                 - Flow % grad_c2c(5,c) * Flow % grad_c2c(5,c) ) / jac
    g_inv(3) = +(  Flow % grad_c2c(1,c) * Flow % grad_c2c(2,c)  &
                 - Flow % grad_c2c(4,c) * Flow % grad_c2c(4,c) ) / jac
    g_inv(4) = -(  Flow % grad_c2c(4,c) * Flow % grad_c2c(3,c)  &
                 - Flow % grad_c2c(5,c) * Flow % grad_c2c(6,c) ) / jac
    g_inv(5) = +(  Flow % grad_c2c(4,c) * Flow % grad_c2c(6,c)  &
                 - Flow % grad_c2c(5,c) * Flow % grad_c2c(2,c) ) / jac
    g_inv(6) = -(  Flow % grad_c2c(1,c) * Flow % grad_c2c(6,c)  &
                 - Flow % grad_c2c(4,c) * Flow % grad_c2c(5,c) ) / jac

    Flow % grad_c2c(1,c) = g_inv(1)
    Flow % grad_c2c(2,c) = g_inv(2)
    Flow % grad_c2c(3,c) = g_inv(3)
    Flow % grad_c2c(4,c) = g_inv(4)
    Flow % grad_c2c(5,c) = g_inv(5)
    Flow % grad_c2c(6,c) = g_inv(6)
  end do

  if(DEBUG) then
    do c = Cells_In_Domain_And_Buffers()
      g1(c) = Flow % grad_c2c(1,c)
      g2(c) = Flow % grad_c2c(2,c)
      g3(c) = Flow % grad_c2c(3,c)
      g4(c) = Flow % grad_c2c(4,c)
      g5(c) = Flow % grad_c2c(5,c)
      g6(c) = Flow % grad_c2c(6,c)
    end do
    call Grid % Save_Debug_Vtu("grad_matrix_with_front",                  &
                               tensor_cell = (/g1, g2, g3, g4, g5, g6/),  &
                               tensor_comp = 6,                           &
                               tensor_name = "grad_matrix_with_front")

    call Work % Disconnect_Real_Cell(g1, g2, g3, g4, g5, g6)
  end if

  end subroutine
