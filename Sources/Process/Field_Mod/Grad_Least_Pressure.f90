!==============================================================================!
  subroutine Field_Mod_Grad_Least_Pressure(flow, p)
!------------------------------------------------------------------------------!
!   Calculates gradient of pressure with least squares method
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  type(Var_Type)   :: p
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c1, c2, iter
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  ! Refresh buffers for variable
  call Grid_Mod_Exchange_Cells_Real(grid, p % n)

  do iter = 1, flow % least_miter

    ! Extrapolation to boundaries
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. PRESSURE) then
          p % n(c2) = p % n(c1) + p % x(c1) * grid % dx(s)  &
                                + p % y(c1) * grid % dy(s)  &
                                + p % z(c1) * grid % dz(s)
        end if
      end if
    end do

    ! Compute individual gradients without refreshing buffers
    call Field_Mod_Grad_Component_No_Refresh(flow, p % n, 1, p % x)  ! dp/dx
    call Field_Mod_Grad_Component_No_Refresh(flow, p % n, 2, p % y)  ! dp/dy
    call Field_Mod_Grad_Component_No_Refresh(flow, p % n, 3, p % z)  ! dp/dz

    ! Refresh buffers for gradient components
    call Grid_Mod_Exchange_Cells_Real(grid, p % x)
    call Grid_Mod_Exchange_Cells_Real(grid, p % y)
    call Grid_Mod_Exchange_Cells_Real(grid, p % z)

  end do

  end subroutine
