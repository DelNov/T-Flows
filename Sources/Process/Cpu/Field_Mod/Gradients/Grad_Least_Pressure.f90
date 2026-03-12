!==============================================================================!
  subroutine Grad_Least_Pressure(Flow, p)
!------------------------------------------------------------------------------!
!>  This subroutine calculates gradient of pressure (or pressure correction)
!>  using the least squares method.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Buffer refresh: The subroutine starts by refreshing the buffers for the  !
!     pressure (or pressure correction) variable p % n.                        !
!   * Iterative calculation: It then enters an iterative loop (up to a         !
!     maximum number of iterations defined by Flow % least_miter):             !
!     - Boundary extrapolation: For each boundary region which is not of the   !
!       PRESSURE type, it extrapolates the pressure values to boundary cells   !
!       based on the neighboring cells' pressure and gradient values.          !
!     - Gradient calculation: Individual gradient components (p % x, p % y,    !
!       and p % z for dp/dx, dp/dy, dp/dz respectively) are calculated for     !
!       each cell in the domain without refreshing buffers.                    !
!     - Buffer refresh for gradient components: After calculating the gradient !
!       components, the subroutine refreshes the buffers for each component.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
  type(Var_Type)    :: p
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: s, c1, c2, iter, reg
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  ! Refresh buffers for variable
  call Grid % Exchange_Cells_Real(p % n)

  do iter = 1, Flow % least_miter

    ! Extrapolation to boundaries
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .ne. AMBIENT) then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          p % n(c2) = p % n(c1) + p % x(c1) * Grid % dx(s)  &
                                + p % y(c1) * Grid % dy(s)  &
                                + p % z(c1) * Grid % dz(s)
        end do  ! faces
      else
        ! This was used for checking
        ! do s = Faces_In_Region(reg)
        !   c2 = Grid % faces_c(2,s)
        ! end do
      end if    ! pressure, ambient or other
    end do      ! regions

    ! Compute individual gradients without refreshing buffers
    call Flow % Grad_Component_No_Refresh(Grid, p % n, 1, p % x)  ! dp/dx
    call Flow % Grad_Component_No_Refresh(Grid, p % n, 2, p % y)  ! dp/dy
    call Flow % Grad_Component_No_Refresh(Grid, p % n, 3, p % z)  ! dp/dz

    ! Refresh buffers for gradient components
    call Grid % Exchange_Cells_Real(p % x)
    call Grid % Exchange_Cells_Real(p % y)
    call Grid % Exchange_Cells_Real(p % z)

  end do

  end subroutine
