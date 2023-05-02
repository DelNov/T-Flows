!==============================================================================!
  subroutine Grad_Least_Pressure(Flow, p)
!------------------------------------------------------------------------------!
!   Calculates gradient of pressure with least squares method
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
      if(Grid % region % type(reg) .ne. PRESSURE) then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          p % n(c2) = p % n(c1) + p % x(c1) * Grid % dx(s)  &
                                + p % y(c1) * Grid % dy(s)  &
                                + p % z(c1) * Grid % dz(s)
        end do  ! faces
      end if    ! boundary not pressure
    end do      ! regions

    ! Compute individual gradients without refreshing buffers
    call Flow % Grad_Component_No_Refresh(p % n, 1, p % x)  ! dp/dx
    call Flow % Grad_Component_No_Refresh(p % n, 2, p % y)  ! dp/dy
    call Flow % Grad_Component_No_Refresh(p % n, 3, p % z)  ! dp/dz

    ! Refresh buffers for gradient components
    call Grid % Exchange_Cells_Real(p % x)
    call Grid % Exchange_Cells_Real(p % y)
    call Grid % Exchange_Cells_Real(p % z)

  end do

  end subroutine
