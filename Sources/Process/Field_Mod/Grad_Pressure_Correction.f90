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
  integer                  :: s, c, c1, c2
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  call Comm_Mod_Exchange_Real(grid, pp % n)

  do c = 1, grid % n_cells
    pp % x(c) = 0.0
    pp % y(c) = 0.0
    pp % z(c) = 0.0
  end do

  !---------------------------------!
  !   No correction at boundaries   !
  !---------------------------------!

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. PRESSURE) then
        pp % n(c2) = 0.0
      end if
    end if
  end do

  call Field_Mod_Grad_Component(flow, pp % n, 1, pp % x)  ! dp/dx
  call Field_Mod_Grad_Component(flow, pp % n, 2, pp % y)  ! dp/dy
  call Field_Mod_Grad_Component(flow, pp % n, 3, pp % z)  ! dp/dz

  end subroutine
