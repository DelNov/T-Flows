!==============================================================================!
  subroutine Calculate_Grad_Matrix(Flow)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow  !! parent flow object
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer  :: Grid
  integer                   :: c, c1, c2, s
  real                      :: dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2
  real                      :: jac, g_inv(6)
!==============================================================================!

  call Profiler % Start('Calculate_Grad_Matrix')

  ! Take alias
  Grid => Flow % pnt_grid

  !--------------------------------------------!
  !   Initialize gradient matrices for cells   !
  !--------------------------------------------!
  do c = 1, Grid % n_cells
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
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 .lt. 0) then

      ! This is not correct on boundaries
      dx_c1 = Grid % dx(s)
      dy_c1 = Grid % dy(s)
      dz_c1 = Grid % dz(s)

      Flow % grad_c2c(1,c1)=Flow % grad_c2c(1,c1) + dx_c1*dx_c1  ! 1,1
      Flow % grad_c2c(2,c1)=Flow % grad_c2c(2,c1) + dy_c1*dy_c1  ! 2,2
      Flow % grad_c2c(3,c1)=Flow % grad_c2c(3,c1) + dz_c1*dz_c1  ! 3,3
      Flow % grad_c2c(4,c1)=Flow % grad_c2c(4,c1) + dx_c1*dy_c1  ! 1,2  &  2,1
      Flow % grad_c2c(5,c1)=Flow % grad_c2c(5,c1) + dx_c1*dz_c1  ! 1,3  &  3,1
      Flow % grad_c2c(6,c1)=Flow % grad_c2c(6,c1) + dy_c1*dz_c1  ! 2,3  &  3,2

    end if
  end do

  ! Faces inside the domain
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 .gt. 0) then

      ! This is correct only on uniform grids
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

    end if
  end do

  !--------------------------------!
  !   Find the inverse of matrix   !
  !--------------------------------!
  do c = 1, Grid % n_cells
    jac = Flow % grad_c2c(1,c) * Flow % grad_c2c(2,c) * Flow % grad_c2c(3,c)  &
        - Flow % grad_c2c(1,c) * Flow % grad_c2c(6,c) * Flow % grad_c2c(6,c)  &
        - Flow % grad_c2c(4,c) * Flow % grad_c2c(4,c) * Flow % grad_c2c(3,c)  &
        + Flow % grad_c2c(4,c) * Flow % grad_c2c(5,c) * Flow % grad_c2c(6,c)  &
        + Flow % grad_c2c(4,c) * Flow % grad_c2c(5,c) * Flow % grad_c2c(6,c)  &
        - Flow % grad_c2c(5,c) * Flow % grad_c2c(5,c) * Flow % grad_c2c(2,c)

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

  call Profiler % Stop('Calculate_Grad_Matrix')

  end subroutine
