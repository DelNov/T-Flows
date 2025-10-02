!==============================================================================!
  subroutine Calculate_Grad_Matrix(Flow)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to calculate the gradient matrix for all cells
!>  in a computational domain. This matrix is fundamental for calculation of
!>  gradients with the least squares method.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization: The subroutine begins by initializing the gradient       !
!     matrices for cells. This includes setting initial values to zero for     !
!     each cell in the domain.                                                 !
!   * Gradient matrix computation:                                             !
!     - The subroutine iteratively computes the gradient matrix by browsing    !
!       through the faces of each cell. It covers both boundary cells and      !
!       faces inside the domain.                                               !
!     - For each cell, it calculates the components of the gradient matrix     !
!       based on the cell's dimensions and position relative to its            !
!       neighboring cells.                                                     !
!     - This involves accumulating contributions from each face of the cell    !
!       to the respective components of the gradient matrix.                   !
!   * Matrix inversion:                                                        !
!     - Once the gradient matrix for each cell is computed, the subroutine     !
!       calculates the inverse of this matrix.                                 !
!     - The inversion is necessary for later computations in the fluid flow    !
!       simulation, particularly for calculating gradients of flow properties  !
!       like velocity and pressure.                                            !
!   * Debugging upport:                                                        !
!     - The subroutine includes a DEBUG flag that, when enabled, triggers the  !
!       recording and saving of the gradient matrices for further inspection.  !
!------------------------------------------------------------------------------!
!   Note                                                                       !
!                                                                              !
!   * It Has a sister in Vof_Mod/Utilitis/Calculate_Grad_Matrix_With_Front     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow  !! parent flow object
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer  :: Grid
  integer                   :: c, c1, c2, s, reg
  real                      :: dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2
  real                      :: jac, g_inv(6)
  real, contiguous, pointer :: g1(:), g2(:), g3(:), g4(:), g5(:), g6(:)
!==============================================================================!

  if(DEBUG) then
    call Work % Connect_Real_Cell(g1, g2, g3, g4, g5, g6)
  end if

  ! Take alias
  Grid => Flow % pnt_grid

  !--------------------------------------------!
  !   Initialize gradient matrices for cells   !
  !--------------------------------------------!
  do c = Cells_In_Domain()
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

  ! Boundary cells
  do reg = Boundary_Regions()
    do s = Faces_In_Region(reg)
      c1 = Grid % faces_c(1,s)

      dx_c1 = Grid % dx(s)
      dy_c1 = Grid % dy(s)
      dz_c1 = Grid % dz(s)

      Flow % grad_c2c(1,c1)=Flow % grad_c2c(1,c1) + dx_c1*dx_c1  ! 1,1
      Flow % grad_c2c(2,c1)=Flow % grad_c2c(2,c1) + dy_c1*dy_c1  ! 2,2
      Flow % grad_c2c(3,c1)=Flow % grad_c2c(3,c1) + dz_c1*dz_c1  ! 3,3
      Flow % grad_c2c(4,c1)=Flow % grad_c2c(4,c1) + dx_c1*dy_c1  ! 1,2  &  2,1
      Flow % grad_c2c(5,c1)=Flow % grad_c2c(5,c1) + dx_c1*dz_c1  ! 1,3  &  3,1
      Flow % grad_c2c(6,c1)=Flow % grad_c2c(6,c1) + dy_c1*dz_c1  ! 2,3  &  3,2
    end do
  end do

  ! Faces inside the domain
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    dx_c1 = Grid % dx(s)
    dy_c1 = Grid % dy(s)
    dz_c1 = Grid % dz(s)
    dx_c2 = Grid % dx(s)
    dy_c2 = Grid % dy(s)
    dz_c2 = Grid % dz(s)

    Flow % grad_c2c(1,c1)=Flow % grad_c2c(1,c1) + dx_c1*dx_c1  ! 1,1
    Flow % grad_c2c(2,c1)=Flow % grad_c2c(2,c1) + dy_c1*dy_c1  ! 2,2
    Flow % grad_c2c(3,c1)=Flow % grad_c2c(3,c1) + dz_c1*dz_c1  ! 3,3
    Flow % grad_c2c(4,c1)=Flow % grad_c2c(4,c1) + dx_c1*dy_c1  ! 1,2  &  2,1
    Flow % grad_c2c(5,c1)=Flow % grad_c2c(5,c1) + dx_c1*dz_c1  ! 1,3  &  3,1
    Flow % grad_c2c(6,c1)=Flow % grad_c2c(6,c1) + dy_c1*dz_c1  ! 2,3  &  3,2

    Flow % grad_c2c(1,c2)=Flow % grad_c2c(1,c2) + dx_c2*dx_c2  ! 1,1
    Flow % grad_c2c(2,c2)=Flow % grad_c2c(2,c2) + dy_c2*dy_c2  ! 2,2
    Flow % grad_c2c(3,c2)=Flow % grad_c2c(3,c2) + dz_c2*dz_c2  ! 3,3
    Flow % grad_c2c(4,c2)=Flow % grad_c2c(4,c2) + dx_c2*dy_c2  ! 1,2  &  2,1
    Flow % grad_c2c(5,c2)=Flow % grad_c2c(5,c2) + dx_c2*dz_c2  ! 1,3  &  3,1
    Flow % grad_c2c(6,c2)=Flow % grad_c2c(6,c2) + dy_c2*dz_c2  ! 2,3  &  3,2

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
    Assert(jac > 0.0)

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
    do c = Cells_In_Domain()
      g1(c) = Flow % grad_c2c(1,c)
      g2(c) = Flow % grad_c2c(2,c)
      g3(c) = Flow % grad_c2c(3,c)
      g4(c) = Flow % grad_c2c(4,c)
      g5(c) = Flow % grad_c2c(5,c)
      g6(c) = Flow % grad_c2c(6,c)
    end do
    call Grid % Save_Debug_Vtu("grad_matrix",                             &
                               tensor_cell = (/g1, g2, g3, g4, g5, g6/),  &
                               tensor_comp = 6,                           &
                               tensor_name = "grad_matrix")

    call Work % Disconnect_Real_Cell(g1, g2, g3, g4, g5, g6)
  end if

  end subroutine
