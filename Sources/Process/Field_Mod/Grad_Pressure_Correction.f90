!==============================================================================!
  subroutine Field_Mod_Grad_Pressure_Correction(flow, pp)
!------------------------------------------------------------------------------!
!   Calculates gradient of pressure correction.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  type(Var_Type)   :: pp
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c1, c2, iter
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  ! Refresh buffers for variable
  call Grid_Mod_Exchange_Cells_Real(grid, pp % n)

  do iter = 1, 4

    ! Extrapolation to boundaries
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. PRESSURE) then
          pp % n(c2) = pp % n(c1) + pp % x(c1) * grid % dx(s)  &
                                  + pp % y(c1) * grid % dy(s)  &
                                  + pp % z(c1) * grid % dz(s)
        end if
      end if
    end do

    ! Compute individual gradients without refreshing buffers
    call Field_Mod_Grad_Component_No_Refresh(flow, pp % n, 1, pp % x)  ! dp/dx
    call Field_Mod_Grad_Component_No_Refresh(flow, pp % n, 2, pp % y)  ! dp/dy
    call Field_Mod_Grad_Component_No_Refresh(flow, pp % n, 3, pp % z)  ! dp/dz

    ! Refresh buffers for gradient components
    call Grid_Mod_Exchange_Cells_Real(grid, pp % x)
    call Grid_Mod_Exchange_Cells_Real(grid, pp % y)
    call Grid_Mod_Exchange_Cells_Real(grid, pp % z)

  end do

  end subroutine
