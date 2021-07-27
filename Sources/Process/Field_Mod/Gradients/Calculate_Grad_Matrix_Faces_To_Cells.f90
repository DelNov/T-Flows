!==============================================================================!
  subroutine Calculate_Grad_Matrix_Faces_To_Cells(Flow)
!------------------------------------------------------------------------------!
!   Calculates gradient matrix for computation of gradients in cells from      !
!   values in faces.  This function is very similar to more usual cell-based   !
!   version for calculation of gradient matrices, and was introduced as an     !
!   addition, or alternative, to Gauss based method of gradient computation    !
!   It is considered experimental, and I am not sure if it will ever be used   !
!   but could prove useful when implementing force-based face interpolations.  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c, c1, c2, s
  real                     :: dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2
  real                     :: jac, g_inv(6)
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  !--------------------------------------------!
  !   Initialize gradient matrices for cells   !
  !--------------------------------------------!
  do c = 1, Grid % n_cells
    Flow % grad_f2c(1,c) = 0.0
    Flow % grad_f2c(2,c) = 0.0
    Flow % grad_f2c(3,c) = 0.0
    Flow % grad_f2c(4,c) = 0.0
    Flow % grad_f2c(5,c) = 0.0
    Flow % grad_f2c(6,c) = 0.0
  end do

  !----------------------------------------------------------------------!
  !   Compute the gradient matrix for all cells browsing through faces   !
  !----------------------------------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    dx_c1 = Grid % xf(s) - Grid % xc(c1)
    dy_c1 = Grid % yf(s) - Grid % yc(c1)
    dz_c1 = Grid % zf(s) - Grid % zc(c1)
    dx_c2 = Grid % xf(s) - Grid % xc(c2)
    dy_c2 = Grid % yf(s) - Grid % yc(c2)
    dz_c2 = Grid % zf(s) - Grid % zc(c2)

    Flow % grad_f2c(1,c1)=Flow % grad_f2c(1,c1) + dx_c1*dx_c1    ! 1,1
    Flow % grad_f2c(2,c1)=Flow % grad_f2c(2,c1) + dy_c1*dy_c1    ! 2,2
    Flow % grad_f2c(3,c1)=Flow % grad_f2c(3,c1) + dz_c1*dz_c1    ! 3,3
    Flow % grad_f2c(4,c1)=Flow % grad_f2c(4,c1) + dx_c1*dy_c1    ! 1,2  &  2,1
    Flow % grad_f2c(5,c1)=Flow % grad_f2c(5,c1) + dx_c1*dz_c1    ! 1,3  &  3,1
    Flow % grad_f2c(6,c1)=Flow % grad_f2c(6,c1) + dy_c1*dz_c1    ! 2,3  &  3,2
    if(c2 > 0) then  ! this is enough even for parallel
      Flow % grad_f2c(1,c2)=Flow % grad_f2c(1,c2) + dx_c2*dx_c2  ! 1,1
      Flow % grad_f2c(2,c2)=Flow % grad_f2c(2,c2) + dy_c2*dy_c2  ! 2,2
      Flow % grad_f2c(3,c2)=Flow % grad_f2c(3,c2) + dz_c2*dz_c2  ! 3,3
      Flow % grad_f2c(4,c2)=Flow % grad_f2c(4,c2) + dx_c2*dy_c2  ! 1,2  &  2,1
      Flow % grad_f2c(5,c2)=Flow % grad_f2c(5,c2) + dx_c2*dz_c2  ! 1,3  &  3,1
      Flow % grad_f2c(6,c2)=Flow % grad_f2c(6,c2) + dy_c2*dz_c2  ! 2,3  &  3,2
    end if

  end do

  !--------------------------------!
  !   Find the inverse of matrix   !
  !--------------------------------!
  do c = 1, Grid % n_cells
    jac = Flow % grad_f2c(1,c) * Flow % grad_f2c(2,c) * Flow % grad_f2c(3,c)  &
        - Flow % grad_f2c(1,c) * Flow % grad_f2c(6,c) * Flow % grad_f2c(6,c)  &
        - Flow % grad_f2c(4,c) * Flow % grad_f2c(4,c) * Flow % grad_f2c(3,c)  &
        + Flow % grad_f2c(4,c) * Flow % grad_f2c(5,c) * Flow % grad_f2c(6,c)  &
        + Flow % grad_f2c(4,c) * Flow % grad_f2c(5,c) * Flow % grad_f2c(6,c)  &
        - Flow % grad_f2c(5,c) * Flow % grad_f2c(5,c) * Flow % grad_f2c(2,c)

    g_inv(1) = +(  Flow % grad_f2c(2,c) * Flow % grad_f2c(3,c)  &
                 - Flow % grad_f2c(6,c) * Flow % grad_f2c(6,c) ) / (jac+TINY)
    g_inv(2) = +(  Flow % grad_f2c(1,c) * Flow % grad_f2c(3,c)  &
                 - Flow % grad_f2c(5,c) * Flow % grad_f2c(5,c) ) / (jac+TINY)
    g_inv(3) = +(  Flow % grad_f2c(1,c) * Flow % grad_f2c(2,c)  &
                 - Flow % grad_f2c(4,c) * Flow % grad_f2c(4,c) ) / (jac+TINY)
    g_inv(4) = -(  Flow % grad_f2c(4,c) * Flow % grad_f2c(3,c)  &
                 - Flow % grad_f2c(5,c) * Flow % grad_f2c(6,c) ) / (jac+TINY)
    g_inv(5) = +(  Flow % grad_f2c(4,c) * Flow % grad_f2c(6,c)  &
                 - Flow % grad_f2c(5,c) * Flow % grad_f2c(2,c) ) / (jac+TINY)
    g_inv(6) = -(  Flow % grad_f2c(1,c) * Flow % grad_f2c(6,c)  &
                 - Flow % grad_f2c(4,c) * Flow % grad_f2c(5,c) ) / (jac+TINY)

    Flow % grad_f2c(1,c) = g_inv(1)
    Flow % grad_f2c(2,c) = g_inv(2)
    Flow % grad_f2c(3,c) = g_inv(3)
    Flow % grad_f2c(4,c) = g_inv(4)
    Flow % grad_f2c(5,c) = g_inv(5)
    Flow % grad_f2c(6,c) = g_inv(6)
  end do

  end subroutine
