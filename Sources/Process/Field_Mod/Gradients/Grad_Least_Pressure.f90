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
  integer                  :: s, c1, c2, iter
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  ! Refresh buffers for variable
  call Grid % Exchange_Cells_Real(p % n)

  do iter = 1, Flow % least_miter

    ! Extrapolation to boundaries
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid % Bnd_Cond_Type(c2) .ne. PRESSURE) then
          p % n(c2) = p % n(c1) + p % x(c1) * Grid % dx(s)  &
                                + p % y(c1) * Grid % dy(s)  &
                                + p % z(c1) * Grid % dz(s)
        end if
      end if
    end do

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
